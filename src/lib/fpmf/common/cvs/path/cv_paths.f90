!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

! PATHS - path-progress coordinate for PATH-based path

module cv_paths

use pmf_constants
use pmf_dat
use cv_common
use pmf_paths
use pmf_sizes

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePATHS

    character(PMF_MAX_CV_NAME)  :: srcpathname      ! path name
    integer                     :: srcpathindx      ! path index
    class(PathType),pointer     :: srcpath          ! path data

    real(PMFDP)                 :: alpha            ! unit less
    integer                     :: ioffset

    real(PMFDP),allocatable     :: dsc(:)           ! helper vector for derivatives

    contains
        procedure :: load_cv        => load_paths
        procedure :: join_cv2path   => join_cv2path_paths
        procedure :: calculate_cv   => calculate_paths

end type CVTypePATHS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_paths
!===============================================================================

subroutine load_paths(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePATHS)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------
    logical                             :: rst
    ! --------------------------------------------------------------------------

! simple init and allocation --------------------
    cv_item%ctype         = 'PATHS' 
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)
    cv_item%requirepath   = .true.

! read path name  --------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'path',cv_item%srcpathname)) then
        call pmf_utils_exit(PMF_OUT,1,'[PATHS] The PATH name (path) is not provided!')
    end if

    write(PMF_OUT,10) trim(cv_item%srcpathname)

! read alpha -------------------------------------
    cv_item%alpha = 0.0d0
    rst =  prmfile_get_real8_by_key(prm_fin,'alpha',cv_item%alpha)

    if( cv_item%alpha .gt. 0.0d0 ) then
        write(PMF_OUT,20) cv_item%alpha
    else
        write(PMF_OUT,30) 
    end if   

    ! zero value will be handled in attach_cv2path_paths()

! read offset
    cv_item%ioffset       = 1
    if( prmfile_get_integer_by_key(prm_fin,'ioffset',cv_item%ioffset) ) then
        write(PMF_OUT,40) cv_item%ioffset
    end if

10 format('   ** Path               : ',A)
20 format('   ** Alpha              : ',F6.3)
30 format('   ** Alpha              : *auto*')
40 format('   ** IOffset            : ',I2)

end subroutine load_paths

!===============================================================================
! Subroutine:  join_cv2path_paths
!===============================================================================

subroutine join_cv2path_paths(cv_item)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePATHS)  :: cv_item
    ! --------------------------------------------
    integer             :: i,j,alloc_failed
    real(PMFDP)         :: ndl,dl2,s,min,max,valn1,valn2
    ! --------------------------------------------------------------------------

! get the path
    cv_item%srcpathindx = pmf_paths_find_path(cv_item%srcpathname)
    cv_item%srcpath => PathList(cv_item%srcpathindx)%path

! init algebraic indexes
    cv_item%isalgebraic = .true.

    allocate(cv_item%algebraicidxs(cv_item%srcpath%ncvs), &
             cv_item%dsc(cv_item%srcpath%ncvs), stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for algebraicidxs and/or dsc in join_cv2path_paths!')
    endif

    do i=1,cv_item%srcpath%ncvs
        cv_item%algebraicidxs(i) = cv_item%srcpath%cvindxs(i)
        if( cv_item%srcpath%cvs(i)%cv%gradforanycrd .neqv. .true. ) then
            call pmf_utils_exit(PMF_OUT,1,'Collective variable ''' // &
                            trim(cv_item%srcpath%cvs(i)%cv%name) // &
                            ''' does not have gradforanycrd == true in join_cv2path_paths!')
        end if
    end do

! calculate alpha if requested
    if( cv_item%alpha .le. 0.0d0 ) then
        ndl = 0
        do i=1,cv_item%srcpath%nbeads-1
            dl2 = 0.0d0
            do j=1,cv_item%srcpath%ncvs
                min = cv_item%srcpath%minvalues(j)
                max = cv_item%srcpath%maxvalues(j)
                valn1 = cv_item%srcpath%points(i+0,j)
                valn2 = cv_item%srcpath%points(i+1,j)
                s = (valn2 - valn1) / (max - min)
                dl2 = dl2 + s**2
            end do
            cv_item%alpha = cv_item%alpha + sqrt(dl2)
            ndl = ndl + 1.0d0
        end do
        if( ndl .gt. 0.0d0 ) then
            cv_item%alpha = cv_item%alpha / ndl
        end if
    end if

    write(PMF_OUT,20) cv_item%alpha

20 format('   ** Alpha              : ',F6.3)

end subroutine join_cv2path_paths

!===============================================================================
! Subroutine:  calculate_paths
!===============================================================================

subroutine calculate_paths(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypePATHS)      :: cv_item
    real(PMFDP)             :: x(:,:)
    type(CVContextType)     :: ctx
    ! -----------------------------------------------
    integer                 :: i,j
    real(PMFDP)             :: r2,vala,valr,s
    real(PMFDP)             :: exponent,max_exponent,weight,sumw,sumkw,meank
    real(PMFDP)             :: range,alpha2,ki
    ! --------------------------------------------------------------------------

    if( cv_item%srcpath%nbeads .le. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PATHS] At least two path beads are required!')
    end if

    alpha2 = cv_item%alpha**2
    if( alpha2 .le. tiny(alpha2) ) then
        call pmf_utils_exit(PMF_OUT,1,'[PATHS] Alpha must be positive and sufficiently large!')
    end if

! Find the largest exponent.  Subtracting it below implements the
! log-sum-exp/softmax stabilization.  The largest shifted weight is then one,
! so the denominator cannot underflow even far away from the path.
    max_exponent = -huge(1.0_PMFDP)

    do i=1,cv_item%srcpath%nbeads
        r2 = 0.0d0
        do j=1,cv_item%srcpath%ncvs
            range = cv_item%srcpath%maxvalues(j) - cv_item%srcpath%minvalues(j)
            if( abs(range) .le. tiny(range) ) then
                call pmf_utils_exit(PMF_OUT,1,'[PATHS] Invalid zero range of a source collective variable!')
            end if
            valr = cv_item%srcpath%points(i,j)
            vala = ctx%CVsValues(cv_item%srcpath%cvs(j)%cv%idx)
            s = (vala - valr) / range
            r2 = r2 + s*s
        end do
        exponent = -r2/alpha2
        max_exponent = max(max_exponent,exponent)
    end do

! Calculate normalized weighted mean of the bead index.  The common factor
! exp(max_exponent) cancels exactly from numerator and denominator.
    sumw  = 0.0_PMFDP
    sumkw = 0.0_PMFDP

    do i=1,cv_item%srcpath%nbeads
        r2 = 0.0_PMFDP
        do j=1,cv_item%srcpath%ncvs
            range = cv_item%srcpath%maxvalues(j) - cv_item%srcpath%minvalues(j)
            valr = cv_item%srcpath%points(i,j)
            vala = ctx%CVsValues(cv_item%srcpath%cvs(j)%cv%idx)
            s = (vala - valr) / range
            r2 = r2 + s*s
        end do

        exponent = -r2/alpha2
        weight = exp(exponent-max_exponent)
        ki = real(i-cv_item%ioffset,PMFDP)
        sumw  = sumw  + weight
        sumkw = sumkw + ki*weight
    end do

    meank = sumkw/sumw
    ctx%CVsValues(cv_item%idx) = meank/real(cv_item%srcpath%nbeads-1,PMFDP)

    if( fdebug ) then
        write(PMF_DEBUG,*) 'PATHS max exponent= ', max_exponent
        write(PMF_DEBUG,*) 'PATHS shifted sumw= ', sumw
        write(PMF_DEBUG,*) 'PATHS mean index= ', meank
        write(PMF_DEBUG,*) 'PATHS v= ', ctx%CVsValues(cv_item%idx)
    end if

! The derivative is evaluated as a covariance under the normalized weights:
!
!   d<k>/dq = sum_i p_i (k_i-<k>) d(log w_i)/dq
!
! This avoids the unstable 1/sumw and 1/sumw**2 terms of the quotient rule.
    cv_item%dsc(:) = 0.0_PMFDP

    do i=1,cv_item%srcpath%nbeads
        r2 = 0.0_PMFDP
        do j=1,cv_item%srcpath%ncvs
            range = cv_item%srcpath%maxvalues(j) - cv_item%srcpath%minvalues(j)
            valr = cv_item%srcpath%points(i,j)
            vala = ctx%CVsValues(cv_item%srcpath%cvs(j)%cv%idx)
            s = (vala - valr) / range
            r2 = r2 + s*s
        end do

        exponent = -r2/alpha2
        weight = exp(exponent-max_exponent)/sumw
        ki = real(i-cv_item%ioffset,PMFDP)

        do j=1,cv_item%srcpath%ncvs
            range = cv_item%srcpath%maxvalues(j) - cv_item%srcpath%minvalues(j)
            valr = cv_item%srcpath%points(i,j)
            vala = ctx%CVsValues(cv_item%srcpath%cvs(j)%cv%idx)
            s = (vala - valr) / range

            cv_item%dsc(j) = cv_item%dsc(j) &
                - 2.0_PMFDP*weight*(ki-meank)*s &
                / (alpha2*range*real(cv_item%srcpath%nbeads-1,PMFDP))
        end do
    end do

    do j=1,cv_item%srcpath%ncvs
        if( fdebug ) then
            write(PMF_DEBUG,*) 'PATHS grd= ', cv_item%dsc(j)
        end if
        ctx%CVsDrvs(:,:,cv_item%idx) = ctx%CVsDrvs(:,:,cv_item%idx) &
            + cv_item%dsc(j)*ctx%CVsDrvs(:,:,cv_item%srcpath%cvs(j)%cv%idx)
    end do

    ! disable unused variable warning
    ignored_arg__ = size(x) .ne. 0

 return

end subroutine calculate_paths

!===============================================================================

end module cv_paths

