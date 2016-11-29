!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2016 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_rswitch

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeRSWITCH

    integer     :: arg_cv
    real(PMFDP) :: offset
    real(PMFDP) :: reference
    integer     :: npow
    integer     :: mpow

    contains
        procedure :: load_cv        => load_rswitch
        procedure :: calculate_cv   => calculate_rswitch
end type CVTypeRSWITCH

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_rswitch
!===============================================================================

subroutine load_rswitch(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_unit

    implicit none
    class(CVTypeRSWITCH)                :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: alloc_failed, i
    character(PMF_MAX_CV_NAME)          :: cv_name
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'RSWITCH'
    cv_item%gradforanycrd = .true.
    cv_item%isalgebraic   = .true.
    call pmf_unit_init(cv_item%unit)
    call cv_common_read_name(cv_item,prm_fin)

    ! load sub specific data  -----------------------
    allocate(cv_item%algebraicidxs(1), stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for algebraicidxs!')
    endif

    if( .not. prmfile_get_string_by_key(prm_fin,'arg',cv_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'Argument CV (arg) is not specified!')
    else
        write(PMF_OUT,5) cv_name
        cv_item%arg_cv = cv_common_find_cv(cv_name)
        cv_item%algebraicidxs(1) = cv_item%arg_cv
    end if

    ! this CV does not have any groups
    cv_item%ngrps = 0

    ! set atoms as a copy from left and right CVs, this information is required fro proper detection of
    ! SHAKE constraints in collisions
    cv_item%natoms = CVList(cv_item%arg_cv)%cv%natoms

    allocate(cv_item%rindexes(cv_item%natoms), &
             cv_item%lindexes(cv_item%natoms), &
             stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for atoms!')
    endif

    do i=1,CVList(cv_item%arg_cv)%cv%natoms
        cv_item%rindexes(i) = CVList(cv_item%arg_cv)%cv%rindexes(i)
    end do

    if( CVList(cv_item%arg_cv)%cv%gradforanycrd .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> Argument CV does not have gradforanycrd == true!')
    end if

    ! load rest of setup
    if( .not. prmfile_get_real8_by_key(prm_fin,'offset',cv_item%offset) ) then
        call pmf_utils_exit(PMF_OUT,1,'offset is not specified!')
    else
        write(PMF_OUT,10) cv_item%offset,trim(pmf_unit_label(CVList(cv_item%arg_cv)%cv%unit))
    end if
    call pmf_unit_conv_to_ivalue(CVList(cv_item%arg_cv)%cv%unit,cv_item%offset)

    if( .not. prmfile_get_real8_by_key(prm_fin,'reference',cv_item%reference) ) then
        call pmf_utils_exit(PMF_OUT,1,'reference is not specified!')
    else
        write(PMF_OUT,20) cv_item%reference,trim(pmf_unit_label(CVList(cv_item%arg_cv)%cv%unit))
    end if
    call pmf_unit_conv_to_ivalue(CVList(cv_item%arg_cv)%cv%unit,cv_item%reference)

    if( cv_item%reference .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'reference cannot be zero!')
    end if

    if( .not. prmfile_get_integer_by_key(prm_fin,'npower',cv_item%npow) ) then
        call pmf_utils_exit(PMF_OUT,1,'npower is not specified!')
    else
        write(PMF_OUT,30) cv_item%npow
    end if

    if( cv_item%npow .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'npower cannot be zero!')
    end if

    if( .not. prmfile_get_integer_by_key(prm_fin,'mpower',cv_item%mpow) ) then
        call pmf_utils_exit(PMF_OUT,1,'mpower is not specified!')
    else
        write(PMF_OUT,40) cv_item%mpow
    end if

    if( cv_item%mpow .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'mpower cannot be zero!')
    end if

    return

  5 format('   ** arg CV             : ',A)
 10 format('   ** distance offset    : ',E14.5,' [',A,']')
 20 format('   ** reference value    : ',E14.5,' [',A,']')
 30 format('   ** n-power            : ',I8)
 40 format('   ** m-power            : ',I8)

end subroutine load_rswitch

!===============================================================================
! Subroutine:  calculate_rswitch
!===============================================================================

subroutine calculate_rswitch(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeRSWITCH)    :: cv_item
    real(PMFDP)             :: x(:,:)
    type(CVContextType)     :: ctx
    ! -----------------------------------------------
    real(PMFDP)             :: dv,iref,up,dn,dm
    ! --------------------------------------------------------------------------

    iref = 1.0d0 / cv_item%reference

    dm = ctx%CVsValues(cv_item%arg_cv) - cv_item%offset
    if( dm .le. 0 ) then
        ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + 1.0d0
        ctx%CVsDrvs(:,:,cv_item%idx) = 0.0d0
    else
        if( dm .ne. cv_item%reference ) then
            up = (dm*iref)**cv_item%npow
            dn = (dm*iref)**cv_item%mpow
            ctx%CVsValues(cv_item%idx) = (1.0d0 - up) / (1.0d0 - dn)
            dv =    - cv_item%npow*((dm*iref)**(cv_item%npow-1))*iref*(1.0d0 - dn)
            dv = dv + cv_item%mpow*((dm*iref)**(cv_item%mpow-1))*iref*(1.0d0 - up)
            dv = dv / (1.0d0 - dn)**2
            ctx%CVsDrvs(:,:,cv_item%idx) = ctx%CVsDrvs(:,:,cv_item%idx) + dv*ctx%CVsDrvs(:,:,cv_item%arg_cv)
        else
            ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + real(cv_item%npow) / real(cv_item%mpow)
            ! TODO derivatives ?
            call pmf_utils_exit(PMF_OUT,1,'not implemented B!')
        end if
    end if

    ! disable unused variable warning
    ignored_arg__ = size(x) .ne. 0

end subroutine calculate_rswitch

!===============================================================================

end module cv_rswitch

