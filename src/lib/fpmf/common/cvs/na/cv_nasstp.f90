!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2016 Ivo Durnik,
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

module cv_nasstp

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type
use cv_math

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNASSTP

    type(XYZFILE_TYPE)  :: xyz_str_a1
    type(SImpStrData)   :: simpdat_a1
    type(XYZFILE_TYPE)  :: xyz_str_b1
    type(SImpStrData)   :: simpdat_b1
    type(XYZFILE_TYPE)  :: xyz_str_a2
    type(SImpStrData)   :: simpdat_a2
    type(XYZFILE_TYPE)  :: xyz_str_b2
    type(SImpStrData)   :: simpdat_b2
    integer             :: lst_par

    contains
        procedure :: load_cv        => load_nasstp
        procedure :: calculate_cv   => calculate_nasstp
end type CVTypeNASSTP

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_nasstp
!===============================================================================

subroutine load_nasstp(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeNASSTP)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'NASSTP'
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 8
    call cv_common_init_groups(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

!===============================================================================

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    if( cv_get_group_natoms(cv_item,1) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_a must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_a',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_a) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_a1)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a1,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a1)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a1)

    if( cv_item%xyz_str_a1%natoms .ne. cv_item%grps(1) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group a and reference structure A differs!')
    end if

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a1%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a1%symbols(i))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,60)
    call cv_common_read_group(cv_item,prm_fin,2)

    if( cv_get_group_natoms(cv_item,2) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_b must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_b',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_b) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_b1)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_b1,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_b1)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_b1)

    if( cv_item%xyz_str_b1%natoms .ne. (cv_item%grps(2) - cv_item%grps(1)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group b and reference structure B differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(1)+1, cv_item%grps(2)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b1%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
                write(tmpstr,110) i-cv_item%grps(1), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b1%symbols(i-cv_item%grps(1)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group c,d ----------------------------------
    write(PMF_OUT,75)
    call cv_common_read_group(cv_item,prm_fin,3)
    if( cv_get_group_natoms(cv_item,3) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_c can contain only one atom!')
    end if

    call cv_common_read_group(cv_item,prm_fin,4)
    if( cv_get_group_natoms(cv_item,4) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_d can contain only one atom!')
    end if

!===============================================================================

    ! read group a ----------------------------------
    write(PMF_OUT,55)
    call cv_common_read_group(cv_item,prm_fin,5)

    if( cv_get_group_natoms(cv_item,5) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_e must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_e',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_e) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_a2)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a2,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a2)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a2)

    if( cv_item%xyz_str_a2%natoms .ne. (cv_item%grps(5) - cv_item%grps(4) ) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group e and reference structure E differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(4)+1, cv_item%grps(5)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a2%symbols(i-cv_item%grps(4)))) .gt. 1.0 ) then
                write(tmpstr,120) i-cv_item%grps(4), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a2%symbols(i-cv_item%grps(4)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,65)
    call cv_common_read_group(cv_item,prm_fin,6)

    if( cv_get_group_natoms(cv_item,6) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_f must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_f',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_f) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_b2)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_b2,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_b2)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_b2)

    if( cv_item%xyz_str_b2%natoms .ne. (cv_item%grps(6) - cv_item%grps(5)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group f and reference structure F differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(5)+1, cv_item%grps(6)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b2%symbols(i-cv_item%grps(5)))) .gt. 1.0 ) then
                write(tmpstr,130) i-cv_item%grps(5), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b2%symbols(i-cv_item%grps(5)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group g,h ----------------------------------
    write(PMF_OUT,76)
    call cv_common_read_group(cv_item,prm_fin,7)
    if( cv_get_group_natoms(cv_item,7) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_g can contain only one atom!')
    end if

    call cv_common_read_group(cv_item,prm_fin,8)
    if( cv_get_group_natoms(cv_item,8) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_h can contain only one atom!')
    end if

!===============================================================================

    ! parameter to be calculated ------------------------
    write(PMF_OUT,80)
    if( .not. prmfile_get_string_by_key(prm_fin,'parameter',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'Type of local base pair step parameter (parameter) is not specified!')
    end if
    write(PMF_OUT,90) trim(tmpstr)

    select case(trim(tmpstr))
        case('shift')
            cv_item%unit    = LengthUnit
            cv_item%lst_par= 1
        case('slide')
            cv_item%unit    = LengthUnit
            cv_item%lst_par= 2
        case('rise')
            cv_item%unit    = LengthUnit
            cv_item%lst_par= 3
        case('tilt')
            cv_item%unit    = AngleUnit
            cv_item%lst_par= 4
        case('roll')
            cv_item%unit    = AngleUnit
            cv_item%lst_par= 5
        case('twist')
            cv_item%unit    = AngleUnit
            cv_item%lst_par= 6
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option!')
    end select

    return

 50 format('   == Base #A1 ===================================')
 60 format('   == Base #B1 ===================================')
 55 format('   == Base #A2 ===================================')
 65 format('   == Base #B2 ===================================')
 70 format('   ** reference structure: ',A)
 75 format('   == y1-axis  ===================================')
 76 format('   == y2-axis  ===================================')
 80 format('   -----------------------------------------------')
 90 format('   ** local step parameter : ',A)

100 format('Atom mismatch between group A1 and reference A1 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B1 and reference B1 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

120 format('Atom mismatch between group A2 and reference A2 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
130 format('Atom mismatch between group B2 and reference B2 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

end subroutine load_nasstp

!===============================================================================
! Subroutine:  calculate_nasstp
!===============================================================================

subroutine calculate_nasstp(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    ! --------------------------------------------------------------------------

    call calculate_nasstp_getbp1(cv_item,x,ua,oa)
    call calculate_nasstp_getbp2(cv_item,x,ub,ob)

    call calculate_nasstp_value_num(cv_item,ctx,ua,oa,ub,ob,a_ua,a_oa,a_ub,a_ob)

    write(*,*) ctx%CVsDrvs(:,:,cv_item%idx)

    call calculate_nasstp_getbp1_der(cv_item,x,ctx,a_ua,a_oa)

    write(*,*) ctx%CVsDrvs(:,:,cv_item%idx)

    call calculate_nasstp_getbp2_der(cv_item,x,ctx,a_ub,a_ob)

    write(*,*) ctx%CVsDrvs(:,:,cv_item%idx)

end subroutine calculate_nasstp

!===============================================================================
! Subroutine:  calculate_nasstp_value_num
!===============================================================================

subroutine calculate_nasstp_value_num(cv_item,ctx,ua,oa,ub,ob,a_ua,a_oa,a_ub,a_ob)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    type(CVContextType) :: ctx
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    ! -----------------------------------------------
    real(PMFDP)         :: l_ua(3,3),l_ub(3,3)
    real(PMFDP)         :: l_oa(3),l_ob(3)
    real(PMFDP)         :: l_a_ua(3,3),l_a_ub(3,3)
    real(PMFDP)         :: l_a_oa(3),l_a_ob(3)
    real(PMFDP)         :: d,v1,v2
    integer             :: i,j
    ! --------------------------------------------------------------------------

    d = 1d-4

    l_ua = ua
    l_oa = oa
    l_ub = ub
    l_ob = ob

    ! ua
    do i=1,3
        do j=1,3
            l_ua = ua
            l_ua(i,j) = l_ua(i,j) - d
            call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
            v1 = ctx%CVsValues(cv_item%idx)
            l_ua = ua
            l_ua(i,j) = l_ua(i,j) + d
            call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
            v2 = ctx%CVsValues(cv_item%idx)
            a_ua(i,j) = (v2-v1)/(2.0d0*d)
        end do
    end do

    l_ua = ua
    l_oa = oa
    l_ub = ub
    l_ob = ob

    ! oa
    do i=1,3
        l_oa = oa
        l_oa(i) = l_oa(i) - d
        call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
        v1 = ctx%CVsValues(cv_item%idx)
        l_oa = oa
        l_oa(i) = l_oa(i) + d
        call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
        v2 = ctx%CVsValues(cv_item%idx)
        a_oa(i) = (v2-v1)/(2.0d0*d)
    end do

    l_ua = ua
    l_oa = oa
    l_ub = ub
    l_ob = ob

    ! ub
    do i=1,3
        do j=1,3
            l_ub = ub
            l_ub(i,j) = l_ub(i,j) - d
            call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
            v1 = ctx%CVsValues(cv_item%idx)
            l_ub = ub
            l_ub(i,j) = l_ub(i,j) + d
            call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
            v2 = ctx%CVsValues(cv_item%idx)
            a_ub(i,j) = (v2-v1)/(2.0d0*d)
        end do
    end do

    l_ua = ua
    l_oa = oa
    l_ub = ub
    l_ob = ob

    ! ob
    do i=1,3
        l_ob = ob
        l_ob(i) = l_ob(i) - d
        call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
        v1 = ctx%CVsValues(cv_item%idx)
        l_ob = ob
        l_ob(i) = l_ob(i) + d
        call calculate_nasstp_value(cv_item,ctx,l_ua,l_oa,l_ub,l_ob,l_a_ua,l_a_oa,l_a_ub,l_a_ob)
        v2 = ctx%CVsValues(cv_item%idx)
        a_ob(i) = (v2-v1)/(2.0d0*d)
    end do

end subroutine calculate_nasstp_value_num

!===============================================================================
! Subroutine:  calculate_nasstp_value
!===============================================================================

subroutine calculate_nasstp_value(cv_item,ctx,ua,oa,ub,ob,a_ua,a_oa,a_ub,a_ob)

    use pmf_dat
    use pmf_pbc
    use pmf_utils
    use cv_math

    implicit none
    class(CVTypeNASSTP) :: cv_item
    type(CVContextType) :: ctx
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    ! --------------------------------------------
    real(PMFDP)         :: g,h(3),nh(3),hlen,d(3),arg,sc,phi
    real(PMFDP)         :: rua(3,3),rub(3,3),mst(3,3)
    ! -----------------------------------------------
    real(PMFDP)         :: a_mst(3,3),a_g,a_h(3),a_nh(3),a_phi,a_rua(3,3),a_rub(3,3)
    real(PMFDP)         :: a_ga,a_gb
    ! --------------------------------------------------------------------------

! determine gamma and hinge axis
    call get_nvangle(ua(:,3),ub(:,3),g)
    call get_cross_product(ua(:,3),ub(:,3),h)

! handle situation with aligned z-axis
    call get_vlen(h,hlen)
    if( (abs(g) .le. PMF_MEPS) .or. (abs(PMF_PI-g) .le. PMF_MEPS) .or. (hlen .le. PMF_MEPS) ) then
        h(:) = ua(:,1) + ub(:,1) + ua(:,2) + ub(:,2)
    end if

    call norm_vec(h,nh)

! rotate ua and ub
    call rotate_ux(nh,+0.5d0*g,ua,rua)
    call rotate_ux(nh,-0.5d0*g,ub,rub)

    ! get mst
    call get_mst(rua,rub,mst)

! derivatives ====================================

! final derivatives
    a_ua(:,:)   = 0.0d0
    a_ub(:,:)   = 0.0d0
    a_oa(:)     = 0.0d0
    a_ob(:)     = 0.0d0

! intermediate derivatives
    a_mst(:,:)  = 0.0d0
    a_nh(:)     = 0.0d0
    a_h(:)      = 0.0d0
    a_g         = 0.0d0
    a_phi       = 0.0d0
    a_rua(:,:)  = 0.0d0
    a_rub(:,:)  = 0.d00

! calculate derivatives
    select case(cv_item%lst_par)
        case(1)
            ! 'shift'
            ! vector between origins
            d(:) = ob(:) - oa(:)
            ctx%CVsValues(cv_item%idx) = d(1)*mst(1,1) + d(2)*mst(2,1) + d(3)*mst(3,1)
            ! with respect to oa and ob
            a_oa(:) = a_oa(:) - mst(:,1)
            a_ob(:) = a_ob(:) + mst(:,1)
            ! with respect to mst
            a_mst(:,1) = a_mst(:,1) + d(:)
        case(2)
            ! 'slide'
            ! vector between origins
            d(:) = ob(:) - oa(:)
            ctx%CVsValues(cv_item%idx) = d(1)*mst(1,2) + d(2)*mst(2,2) + d(3)*mst(3,2)
            ! with respect to oa and ob
            a_oa(:) = a_oa(:) - mst(:,2)
            a_ob(:) = a_ob(:) + mst(:,2)
            ! with respect to mst
            a_mst(:,2) = a_mst(:,2) + d(:)
        case(3)
            ! 'rise'
            ! vector between origins
            d(:) = ob(:) - oa(:)
            ctx%CVsValues(cv_item%idx) = d(1)*mst(1,3) + d(2)*mst(2,3) + d(3)*mst(3,3)
            ! with respect to oa and ob
            a_oa(:) = a_oa(:) - mst(:,3)
            a_ob(:) = a_ob(:) + mst(:,3)
            ! with respect to mst
            a_mst(:,3) = a_mst(:,3) + d(:)
        case(4)
            ! 'tilt'
            call get_nvangle(nh,mst(:,2),phi)
            call get_vtors_sign(nh,mst(:,2),mst(:,3),sc)
            ctx%CVsValues(cv_item%idx) = g * sin(phi*sc)
            ! with respect to g
            a_g = sin(phi*sc)
            ! with respect to phi
            a_phi = g*cos(phi*sc)*sc
            ! with respect to h a mst
            call get_nvangle_der(nh,mst(:,2),a_phi,a_nh,a_mst(:,2))
        case(5)
            ! 'roll'
            call get_nvangle(nh,mst(:,2),phi)
            call get_vtors_sign(nh,mst(:,2),mst(:,3),sc)
            ctx%CVsValues(cv_item%idx) = g * cos(phi*sc)
            ! with respect to g
            a_g = cos(phi*sc)
            ! with respect to phi
            a_phi = -g*sin(phi*sc)*sc
            ! with respect to h a mst
            call get_nvangle_der(nh,mst(:,2),a_phi,a_nh,a_mst(:,2))
        case(6)
            ! 'twist'
            call get_nvangle(rua(:,2),rub(:,2),arg)
            call get_vtors_sign(rua(:,2),rub(:,2),mst(:,3),sc)
            ctx%CVsValues(cv_item%idx) = arg * sc
            ! with respect to h a mst
            a_phi = sc
            call get_nvangle_der(rua(:,2),rub(:,2),a_phi,a_rua(:,2),a_rub(:,2))
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option in calculate_nasstp!')
    end select

    !! get mst
    !    call get_mst(rua,rub,mst)
    ! with respect to rua a rub
    call get_mst_der(rua,rub,a_mst,a_rua,a_rub)

    !! rotate ua and ub
    !    call rotate_ux(h,+0.5d0*g,ua,rua)
    !    call rotate_ux(h,-0.5d0*g,ub,rub)

    a_ga = 0.0d0
    a_gb = 0.0d0
    call rotate_ux_der(nh,+0.5d0*g,ua,a_rua,a_nh,a_ga,a_ua)
    call rotate_ux_der(nh,-0.5d0*g,ub,a_rub,a_nh,a_gb,a_ub)
    a_g = a_g + a_ga * 0.5d0 - a_gb * 0.5d0

!! handle situation with aligned z-axis
!    call norm_vec(h)
!    call get_vlen(h,hlen)
!    if( (abs(g) .le. PMF_MEPS) .or. (abs(PMF_PI-g) .le. PMF_MEPS) .or. (hlen .le. PMF_MEPS) ) then
!        h(:) = ua(:,1) + ub(:,1) + ua(:,2) + ub(:,2);
!        call norm_vec(h)
!    end if
!! determine gamma and hinge axis
!    call get_cross_product(ua(:,3),ub(:,3),h)
!

    call norm_vec_der(h,a_nh,a_h)

    if( (abs(g) .le. PMF_MEPS) .or. (abs(PMF_PI-g) .le. PMF_MEPS) .or. (hlen .le. PMF_MEPS) ) then
        ! FIXME
        stop
    else
        call get_cross_product_der(ua(:,3),ub(:,3),a_h,a_ua(:,3),a_ub(:,3))
    end if

!    call get_nvangle(ua(:,3),ub(:,3),g)
    call get_nvangle_der(ua(:,3),ub(:,3),a_g,a_ua(:,3),a_ub(:,3))

end subroutine calculate_nasstp_value

!===============================================================================
! Subroutine:  calculate_nasstp_getbp1
!===============================================================================

subroutine calculate_nasstp_getbp1(cv_item,x,mst,morg)
    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    real(PMFDP)         :: mst(3,3)
    real(PMFDP)         :: morg(3)
    ! -----------------------------------------------
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: y1(3),y2(3)
    integer             :: ai
    ! --------------------------------------------------------------------------

! superimpose bases ==============================
    call superimpose_str(cv_item,0,              cv_item%grps(1),x,cv_item%xyz_str_a1,cv_item%simpdat_a1,ua,oa)
    call superimpose_str(cv_item,cv_item%grps(1),cv_item%grps(2),x,cv_item%xyz_str_b1,cv_item%simpdat_b1,ub,ob)

! assemble mst and morg ==========================
    ai = cv_item%lindexes(cv_item%grps(3))
    y1(:) = x(:,ai)
    ai = cv_item%lindexes(cv_item%grps(4))
    y2(:) = x(:,ai)
    call get_mst_morg(ua,oa,ub,ob,y1,y2,mst,morg)

end subroutine calculate_nasstp_getbp1

!===============================================================================
! Subroutine:  calculate_nasstp_getbp1_der
!===============================================================================

subroutine calculate_nasstp_getbp1_der(cv_item,x,ctx,a_mst,a_morg)
    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    real(PMFDP)         :: a_mst(3,3)
    real(PMFDP)         :: a_morg(3)
    ! -----------------------------------------------
    integer             :: ai
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: y1(3),y2(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    real(PMFDP)         :: a_y1(3),a_y2(3)
    ! --------------------------------------------------------------------------

! input data
    ua(:,:) = cv_item%simpdat_a1%u(:,:)
    ub(:,:) = cv_item%simpdat_b1%u(:,:)

    ai = cv_item%lindexes(cv_item%grps(3))
    y1(:) = x(:,ai)
    ai = cv_item%lindexes(cv_item%grps(4))
    y2(:) = x(:,ai)

    a_ua(:,:) = 0.0d0
    a_oa(:) = 0.0d0
    a_ub(:,:) = 0.0d0
    a_ob(:) = 0.0d0
    a_y1(:) = 0.0d0
    a_y2(:) = 0.0d0

    call get_mst_morg_der(ua,ub,y1,y2,a_mst,a_morg,a_ua,a_oa,a_ub,a_ob,a_y1,a_y2)

! derivatives for superimposed bases
    call superimpose_str_der(cv_item,0,              cv_item%grps(1),ctx,cv_item%xyz_str_a1,cv_item%simpdat_a1,a_ua,a_oa)
    call superimpose_str_der(cv_item,cv_item%grps(1),cv_item%grps(2),ctx,cv_item%xyz_str_b1,cv_item%simpdat_b1,a_ub,a_ob)

! finally gradients for group_c, group_d
    ai = cv_item%lindexes(cv_item%grps(3))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y1(:)

    ai = cv_item%lindexes(cv_item%grps(4))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y2(:)

end subroutine calculate_nasstp_getbp1_der

!===============================================================================
! Subroutine:  calculate_nasstp_getbp2
!===============================================================================

subroutine calculate_nasstp_getbp2(cv_item,x,mst,morg)
    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    real(PMFDP)         :: mst(3,3)
    real(PMFDP)         :: morg(3)
    ! -----------------------------------------------
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: y1(3),y2(3)
    integer             :: ai
    ! --------------------------------------------------------------------------

! superimpose bases ==============================
    call superimpose_str(cv_item,cv_item%grps(4),cv_item%grps(5),x,cv_item%xyz_str_a2,cv_item%simpdat_a2,ua,oa)
    call superimpose_str(cv_item,cv_item%grps(5),cv_item%grps(6),x,cv_item%xyz_str_b2,cv_item%simpdat_b2,ub,ob)

! assemble mst and morg ==========================
    ai = cv_item%lindexes(cv_item%grps(7))
    y1(:) = x(:,ai)
    ai = cv_item%lindexes(cv_item%grps(8))
    y2(:) = x(:,ai)
    call get_mst_morg(ua,oa,ub,ob,y1,y2,mst,morg)

end subroutine calculate_nasstp_getbp2

!===============================================================================
! Subroutine:  calculate_nasstp_getbp2_der
!===============================================================================

subroutine calculate_nasstp_getbp2_der(cv_item,x,ctx,a_mst,a_morg)
    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    real(PMFDP)         :: a_mst(3,3)
    real(PMFDP)         :: a_morg(3)
    ! -----------------------------------------------
    integer             :: ai
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: y1(3),y2(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    real(PMFDP)         :: a_y1(3),a_y2(3)
    ! --------------------------------------------------------------------------

! input data
    ua(:,:) = cv_item%simpdat_a2%u(:,:)
    ub(:,:) = cv_item%simpdat_b2%u(:,:)

    ai = cv_item%lindexes(cv_item%grps(7))
    y1(:) = x(:,ai)
    ai = cv_item%lindexes(cv_item%grps(8))
    y2(:) = x(:,ai)

    a_ua(:,:) = 0.0d0
    a_oa(:) = 0.0d0
    a_ub(:,:) = 0.0d0
    a_ob(:) = 0.0d0
    a_y1(:) = 0.0d0
    a_y2(:) = 0.0d0

    call get_mst_morg_der(ua,ub,y1,y2,a_mst,a_morg,a_ua,a_oa,a_ub,a_ob,a_y1,a_y2)

! derivatives for superimposed bases =============
    call superimpose_str_der(cv_item,cv_item%grps(4),cv_item%grps(5),ctx,cv_item%xyz_str_a2,cv_item%simpdat_a2,a_ua,a_oa)
    call superimpose_str_der(cv_item,cv_item%grps(5),cv_item%grps(6),ctx,cv_item%xyz_str_b2,cv_item%simpdat_b2,a_ub,a_ob)

! finally gradients for group_c, group_d ========
    ai = cv_item%lindexes(cv_item%grps(7))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y1(:)

    ai = cv_item%lindexes(cv_item%grps(8))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y2(:)

end subroutine calculate_nasstp_getbp2_der

!===============================================================================

end module cv_nasstp

