!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@enzim.hu
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

! worm position

module cv_wormp

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeWORMP

    ! plane definition
    integer             :: x_direction  ! plane x-direction
    integer             :: y_direction  ! plane y-direction
    real(PMFDP)         :: seldist      ! selector radius
    real(PMFDP)         :: steepness    ! selector slope

    ! worm setup
    integer             :: nsegs        ! number of segments
    real(PMFDP)         :: totlen       ! total length of worm
    real(PMFDP),pointer :: alphas(:)    ! parametric positions of segments

    ! intermediate results
    real(PMFDP),pointer :: coms(:,:)   ! 3,nsegs+1
    real(PMFDP),pointer :: totmass(:)   ! nsegs+1
    real(PMFDP),pointer :: odist(:)     ! nsegs
    real(PMFDP),pointer :: wdist(:)     ! nsegs
    real(PMFDP),pointer :: wd(:)        ! nsegs

    contains
        procedure :: load_cv        => load_wormp
        procedure :: calculate_cv   => calculate_wormp
end type CVTypeWORMP

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_wormp
!===============================================================================

subroutine load_wormp(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypeWORMP)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: m, i
    logical                             :: found
    character(len=PRMFILE_MAX_LINE)     :: mask
    integer,parameter                   :: group_index = 96   ! ascii code of 'a' - 1
    integer                             :: alloc_failed
    type(UnitType)                      :: steepnessunit
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'WORM'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load number of segments
    if( .not. prmfile_get_integer_by_key(prm_fin,'nsegs',cv_item%nsegs) ) then
        call pmf_utils_exit(PMF_OUT,1,'nsegs - number of segments is not specified!')
    end if

    ! init groups -----------------------------------
    cv_item%ngrps = cv_item%nsegs + 1   ! +1 for the plane
    call cv_common_init_groups_I(cv_item,prm_fin)
    mask = 'plane'
    call cv_common_init_groups_II(cv_item,prm_fin,1,mask)
    do i=1,cv_item%nsegs
        mask = 'segment_'//achar(group_index+i)
        call cv_common_init_groups_II(cv_item,prm_fin,1+i,mask)
    end do
    call cv_common_init_groups_III(cv_item,prm_fin)

! read plane ------------------------------------
    write(PMF_OUT,50)
    mask = 'plane'
    call cv_common_read_group_by_name(cv_item,prm_fin,1,mask)

    ! load x_direction ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'x_direction',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'x_direction atom is not specified!')
    end if
    write(PMF_OUT,70) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%x_direction)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=1,cv_item%grps(1)
        if( cv_item%rindexes(m) .eq. cv_item%x_direction ) then
            cv_item%x_direction = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'x_direction atom has to be member of group_a!')
    end if

    ! load y_direction ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'y_direction',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'y_direction atom is not specified!')
    end if
    write(PMF_OUT,80) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%y_direction)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=1,cv_item%grps(1)
        if( cv_item%rindexes(m) .eq. cv_item%y_direction ) then
            cv_item%y_direction = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'y_direction atom has to be part of group_a!')
    end if

    ! load rest of setup
    if( .not. prmfile_get_real8_by_key(prm_fin,'seldist',cv_item%seldist) ) then
        call pmf_utils_exit(PMF_OUT,1,'seldist is not specified!')
    else
        write(PMF_OUT,90) cv_item%seldist,trim(pmf_unit_label(LengthUnit))
    end if
    call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%seldist)

    steepnessunit = pmf_unit_power_unit(LengthUnit,-1)
    if( .not. prmfile_get_real8_by_key(prm_fin,'steepness',cv_item%steepness) ) then
        call pmf_utils_exit(PMF_OUT,1,'steepness is not specified!')
    else
        write(PMF_OUT,100) cv_item%steepness,trim(pmf_unit_label(steepnessunit))
    end if
    call pmf_unit_conv_to_ivalue(steepnessunit,cv_item%steepness)

! allocate data
    allocate(cv_item%alphas(cv_item%nsegs), cv_item%coms(3,cv_item%nsegs+1), cv_item%totmass(cv_item%nsegs+1), &
             cv_item%odist(cv_item%nsegs), cv_item%wdist(cv_item%nsegs), &
             cv_item%wd(cv_item%nsegs), stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                       'Unable to allocate memory for CV ('//trim(cv_item%name)//')!')
    end if

! read worm ----------------------------------
    write(PMF_OUT,200)
    write(PMF_OUT,210) cv_item%nsegs

    if( .not. prmfile_get_real8_by_key(prm_fin,'totlen',cv_item%totlen) ) then
        call pmf_utils_exit(PMF_OUT,1,'totlen is not specified!')
    else
        write(PMF_OUT,220) cv_item%totlen,trim(pmf_unit_label(LengthUnit))
    end if
    call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%totlen)

    do i=1,cv_item%nsegs
        ! read segment mask
        mask = 'segment_'//achar(group_index+i)
        call cv_common_read_group_by_name(cv_item,prm_fin,i+1,mask)

        ! read alpha
        if( .not. prmfile_get_real8_by_key(prm_fin,'alpha_'//achar(group_index+i),cv_item%alphas(i)) ) then
            call pmf_utils_exit(PMF_OUT,1, &
                            ''//'alpha_'//achar(group_index+i)//' parameter is not specified!')
        end if
        write(PMF_OUT,230) achar(group_index+i),cv_item%alphas(i)
    end do

    return

 50 format('   == Plane ======================================')
 70 format('   ** x-direction atom   = ',A)
 80 format('   ** y-direction atom   = ',A)
 90 format('   ** selector distance  : ',E14.5,' [',A,']')
100 format('   ** selector steepness : ',E14.5,' [',A,']')

200 format('   == Worm =======================================')
210 format('   ** Number of segments : ',I6)
220 format('   ** Total length       : ',F10.3,' [',A,']')
230 format('   ** alpha_',A1,'            : ',F10.3)

end subroutine load_wormp

!===============================================================================
! Subroutine:  calculate_wormp
!===============================================================================

subroutine calculate_wormp(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeWORMP)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,g,s,ai,m,info,orient,mi,mj
    real(PMFDP)    :: dx(3),dzx(3),dzy(3),dzz(3)
    real(PMFDP)    :: a(3,3),a11,a22,a33,a12,a13,a23,eigenvalues(3)
    real(PMFDP)    :: top,down,amass,msign,ac,dzz_s
    real(PMFDP)    :: work(26*3)
    real(PMFDP)    :: v(3,3),api(3,3),cij(3),xij(3,3,3),bint(3,3)
    ! -----------------------------------------------------------------------------

    ! for each group, calculate coms and totmass
    do g = 1, cv_item%ngrps
        call cv_get_group_com(cv_item,g,x,cv_item%coms(:,i),cv_item%totmass(i))
    end do

    ! calculate parameters of plane -----------------
    a11 = 0.0d0
    a22 = 0.0d0
    a33 = 0.0d0
    a12 = 0.0d0
    a13 = 0.0d0
    a23 = 0.0d0

    ! construct matrix:
    do m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        a11 = a11 + amass*(x(1,ai) - cv_item%coms(1,1))**2
        a22 = a22 + amass*(x(2,ai) - cv_item%coms(2,1))**2
        a33 = a33 + amass*(x(3,ai) - cv_item%coms(3,1))**2
        a12 = a12 + amass*(x(1,ai) - cv_item%coms(1,1))*(x(2,ai) - cv_item%coms(2,1))
        a13 = a13 + amass*(x(1,ai) - cv_item%coms(1,1))*(x(3,ai) - cv_item%coms(3,1))
        a23 = a23 + amass*(x(2,ai) - cv_item%coms(2,1))*(x(3,ai) - cv_item%coms(3,1))
    end do

    a(1,1) = a11
    a(2,2) = a22
    a(3,3) = a33

    a(1,2) = a12
    a(1,3) = a13
    a(2,3) = a23

    a(2,1) = a12
    a(3,1) = a13
    a(3,2) = a23

    ! calculate eignevalues and eigenvectors of matrix
    eigenvalues(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 3, a, 3, eigenvalues, work, 26*3, info)

    if( info .ne. 0 ) then
    call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_wormp!')
    end if

    ! determine z-axis orientation
    dzx(:) = x(:,cv_item%lindexes(cv_item%x_direction)) - cv_item%coms(:,1)
    dzy(:) = x(:,cv_item%lindexes(cv_item%y_direction)) - cv_item%coms(:,1)

    ! cross product
    dzz(1) = dzx(2)*dzy(3) - dzx(3)*dzy(2)
    dzz(2) = dzx(3)*dzy(1) - dzx(1)*dzy(3)
    dzz(3) = dzx(1)*dzy(2) - dzx(2)*dzy(1)

    ! norm dzz (this is not neccessary, but it is keep for debug)
    dzz_s = sqrt(dzz(1)**2 + dzz(2)**2 + dzz(3)**2)
    dzz(:) = dzz(:)/dzz_s

    ! calculate distance to the plane, abs(orient) define which solution will be used
    orient = 1

    ! angle between z-axis and plane vector
    ac = dzz(1)*a(1,orient) + dzz(2)*a(2,orient) + dzz(3)*a(3,orient)
    ! sign of orient is used to change the sign of final value
    msign = sign(1.0d0,ac)

    top = 0.0d0
    down = 0.0d0

    do s=1,cv_item%nsegs

        ! direction vector
        dx(:) = cv_item%coms(:,s+1) - cv_item%coms(:,1)

        ! oriented distance
        cv_item%odist(s) = msign*(a(1,orient)*dx(1) + a(2,orient)*dx(2) + a(3,orient)*dx(3))

        ! normal distance for selector
        cv_item%wdist(s) = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)

        ! weight
        cv_item%wd(s)  = 1.0d0 / (1.0d0 + exp(cv_item%steepness*(cv_item%wdist(s) - cv_item%seldist)))

        top = top + cv_item%wd(s)*(cv_item%alphas(s)*cv_item%totlen + cv_item%odist(s))
        down = down + cv_item%wd(s)
    end do

    ctx%CVsValues(cv_item%idx) = top / (cv_item%totlen * down)

!    ! first derivatives ------------------------------

!    ! first part, e.g. a*dx'
!    do  m = 1, cv_item%grps(1)
!        ai = cv_item%lindexes(m)
!        amass = mass(ai)
!        ctx%CVsDrvs(:,ai,cv_item%idx) =  ctx%CVsDrvs(:,ai,cv_item%idx) - msign*a(:,orient)*amass/totmass1
!    end do

!    do  m = cv_item%grps(1) + 1, cv_item%grps(2)
!        ai = cv_item%lindexes(m)
!        amass = mass(ai)
!        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + msign*a(:,orient)*amass/totmass2
!    end do

!    ! second part, e.g. a'dx

!    ! eigenvector derivatives -----------------------

!    ! construct pseudoinverse matrix, api
!    v(:,:) = a(:,:)
!    api(:,:) = 0.0d0
!    do i=1,3
!        if( i .ne. orient ) api(i,i) = 1.0d0/(eigenvalues(i) - eigenvalues(orient))
!    end do
!    call dgemm('N','N',3,3,3,1.0d0,v,3,api,3,0.0d0,bint,3)
!    call dgemm('N','T',3,3,3,1.0d0,bint,3,v,3,0.0d0,api,3)

!    ! and solve system of equations
!    xij(:,:,:) = 0.0d0
!    do mi=1,3
!        do mj=1,3
!            ! construct cij
!            cij(:) = 0.0d0
!            cij(mi) = cij(mi) + a(mj,orient)

!            ! find eigenvector derivatives
!            ! xi contains derivatives of eigenvector by A_ij element
!            call dgemv('N',3,3,-1.0d0,api,3,cij,1,0.0d0,xij(:,mi,mj),1)

!            ! multiply by dx
!            xij(:,mi,mj) = xij(:,mi,mj)*dx(:)*msign
!        end do
!    end do

!    ! and finaly gradients --------------------------
!    do m = 1, cv_item%grps(1)
!        ai = cv_item%lindexes(m)
!        amass = mass(ai)

!        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
!                            + amass*(2.0d0*(x(1,ai) - d1(1))*(xij(1,1,1) + xij(2,1,1) + xij(3,1,1)) &
!                            +       (x(2,ai) - d1(2))*(xij(1,1,2) + xij(2,1,2) + xij(3,1,2)) &
!                            +       (x(2,ai) - d1(2))*(xij(1,2,1) + xij(2,2,1) + xij(3,2,1)) &
!                            +       (x(3,ai) - d1(3))*(xij(1,1,3) + xij(2,1,3) + xij(3,1,3)) &
!                            +       (x(3,ai) - d1(3))*(xij(1,3,1) + xij(2,3,1) + xij(3,3,1)))

!        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
!                            +       amass*((x(1,ai) - d1(1))*(xij(1,1,2) + xij(2,1,2) + xij(3,1,2)) &
!                            +       (x(1,ai) - d1(1))*(xij(1,2,1) + xij(2,2,1) + xij(3,2,1)) &
!                            + 2.0d0*(x(2,ai) - d1(2))*(xij(1,2,2) + xij(2,2,2) + xij(3,2,2)) &
!                            +       (x(3,ai) - d1(3))*(xij(1,2,3) + xij(2,2,3) + xij(3,2,3)) &
!                            +       (x(3,ai) - d1(3))*(xij(1,3,2) + xij(2,3,2) + xij(3,3,2)))

!        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
!                            +       amass*((x(1,ai) - d1(1))*(xij(1,1,3) + xij(2,1,3) + xij(3,1,3)) &
!                            +       (x(1,ai) - d1(1))*(xij(1,3,1) + xij(2,3,1) + xij(3,3,1)) &
!                            +       (x(2,ai) - d1(2))*(xij(1,2,3) + xij(2,2,3) + xij(3,2,3)) &
!                            +       (x(2,ai) - d1(2))*(xij(1,3,2) + xij(2,3,2) + xij(3,3,2)) &
!                            + 2.0d0*(x(3,ai) - d1(3))*(xij(1,3,3) + xij(2,3,3) + xij(3,3,3)))
!    end do

    return

end subroutine calculate_wormp

!===============================================================================

end module cv_wormp

