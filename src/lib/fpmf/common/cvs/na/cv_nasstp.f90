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

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNASSTP

    type(XYZFILE_TYPE)  :: xyz_str_a1
    type(XYZFILE_TYPE)  :: xyz_str_b1
    type(XYZFILE_TYPE)  :: xyz_str_a2
    type(XYZFILE_TYPE)  :: xyz_str_b2
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

    call calculate_nasstp_value(cv_item,ctx,ua,oa,ub,ob,a_ua,a_oa,a_ub,a_ob)

    call calculate_nasstp_getbp1_der(cv_item,x,ctx,a_ua,a_oa)
    call calculate_nasstp_getbp2_der(cv_item,x,ctx,a_ub,a_ob)

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

    implicit none
    class(CVTypeNASSTP) :: cv_item
    type(CVContextType) :: ctx
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    ! -----------------------------------------------
    real(PMFDP)         :: zsc,xsc,o_zaxis2,o_zaxis,o_xaxis2,o_xaxis
    real(PMFDP)         :: zaxisr(3),zaxis(3)
    real(PMFDP)         :: x0axis(3),xaxisr(3),xaxis(3),yaxis(3)
    real(PMFDP)         :: d(3),tmp1(3),tmp2(3),a_tmp1(3),a_tmp2(3)
    real(PMFDP)         :: tmp12,tmp22,o_tmp12,o_tmp22,o_tmp1vtmp2v,arg,asc
    real(PMFDP)         :: f1,argo_tmp12,argo_tmp22,t1,t2
    real(PMFDP)         :: a_xaxis(3),a_yaxis(3),a_zaxis(3),a_xaxisr(3)
    real(PMFDP)         :: a_x0axis(3)
    ! --------------------------------------------------------------------------

! z-axis ===================================================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = ua(:,3) + zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

! x-axis ===================================================================
    xsc = sign(1.0d0,ua(1,1)*ub(1,1)+ua(2,1)*ub(2,1)+ua(3,1)*ub(3,1))
    ! get x-axis as average of two axes
    x0axis(:) = ua(:,1) + xsc*ub(:,1)
    ! remove projections to z-axis
    xaxisr(:) = x0axis(:) - (x0axis(1)*zaxis(1)+x0axis(2)*zaxis(2)+x0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_xaxis2 = 1.0d0 / (xaxisr(1)**2 + xaxisr(2)**2 + xaxisr(3)**2)
    o_xaxis  = sqrt(o_xaxis2)
    xaxis(:) = xaxisr(:) * o_xaxis

! y-axis ===================================================================
    ! is cross product of z and x vectors
    yaxis(1) = zaxis(2)*xaxis(3) - zaxis(3)*xaxis(2)
    yaxis(2) = zaxis(3)*xaxis(1) - zaxis(1)*xaxis(3)
    yaxis(3) = zaxis(1)*xaxis(2) - zaxis(2)*xaxis(1)

! vector between origins
    d(:) = ob(:) - oa(:)

    select case(cv_item%lst_par)
        case(1)
            ! 'shift'
            ctx%CVsValues(cv_item%idx) = d(1)*xaxis(1) + d(2)*xaxis(2) + d(3)*xaxis(3)
        case(2)
            ! 'slide'
            ctx%CVsValues(cv_item%idx) = d(1)*yaxis(1) + d(2)*yaxis(2) + d(3)*yaxis(3)
        case(3)
            ! 'rise'
            ctx%CVsValues(cv_item%idx) = d(1)*zaxis(1) + d(2)*zaxis(2) + d(3)*zaxis(3)
        case(4)
            ! 'tilt'
            tmp1(1) = ua(2,3)*xaxis(3) - ua(3,3)*xaxis(2)
            tmp1(2) = ua(3,3)*xaxis(1) - ua(1,3)*xaxis(3)
            tmp1(3) = ua(1,3)*xaxis(2) - ua(2,3)*xaxis(1)

            tmp2(1) = zsc*ub(2,3)*xaxis(3) - zsc*ub(3,3)*xaxis(2)
            tmp2(2) = zsc*ub(3,3)*xaxis(1) - zsc*ub(1,3)*xaxis(3)
            tmp2(3) = zsc*ub(1,3)*xaxis(2) - zsc*ub(2,3)*xaxis(1)

            tmp12 = tmp1(1)**2 + tmp1(2)**2 + tmp1(3)**2
            tmp22 = tmp2(1)**2 + tmp2(2)**2 + tmp2(3)**2

            o_tmp12 = 1.0d0 / tmp12
            o_tmp22 = 1.0d0 / tmp22

            o_tmp1vtmp2v = sqrt( o_tmp12 * o_tmp22 )

            arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3)) * o_tmp1vtmp2v

            asc = (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*xaxis(1) + &
                  (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*xaxis(2) + &
                  (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*xaxis(3)

            if ( arg .gt.  1.0 ) then
                arg =  1.0
            else if ( arg .lt. -1.0 ) then
                arg = -1.0
            end if

            ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)

        case(5)
            ! 'roll'
            tmp1(1) = ua(2,3)*yaxis(3) - ua(3,3)*yaxis(2)
            tmp1(2) = ua(3,3)*yaxis(1) - ua(1,3)*yaxis(3)
            tmp1(3) = ua(1,3)*yaxis(2) - ua(2,3)*yaxis(1)

            tmp2(1) = zsc*ub(2,3)*yaxis(3) - zsc*ub(3,3)*yaxis(2)
            tmp2(2) = zsc*ub(3,3)*yaxis(1) - zsc*ub(1,3)*yaxis(3)
            tmp2(3) = zsc*ub(1,3)*yaxis(2) - zsc*ub(2,3)*yaxis(1)

            tmp12 = tmp1(1)**2 + tmp1(2)**2 + tmp1(3)**2
            tmp22 = tmp2(1)**2 + tmp2(2)**2 + tmp2(3)**2

            o_tmp12 = 1.0d0 / tmp12
            o_tmp22 = 1.0d0 / tmp22

            o_tmp1vtmp2v = sqrt( o_tmp12 * o_tmp22 )

            arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3)) * o_tmp1vtmp2v

            asc = (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*yaxis(1) + &
                  (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*yaxis(2) + &
                  (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*yaxis(3)

            if ( arg .gt.  1.0 ) then
                arg =  1.0
            else if ( arg .lt. -1.0 ) then
                arg = -1.0
            end if

            ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)
        case(6)
            ! 'twist'
            tmp1(1) = ua(2,2)*zaxis(3) - ua(3,2)*zaxis(2)
            tmp1(2) = ua(3,2)*zaxis(1) - ua(1,2)*zaxis(3)
            tmp1(3) = ua(1,2)*zaxis(2) - ua(2,2)*zaxis(1)

            tmp2(1) = zsc*ub(2,2)*zaxis(3) - zsc*ub(3,2)*zaxis(2)
            tmp2(2) = zsc*ub(3,2)*zaxis(1) - zsc*ub(1,2)*zaxis(3)
            tmp2(3) = zsc*ub(1,2)*zaxis(2) - zsc*ub(2,2)*zaxis(1)

            tmp12 = tmp1(1)**2 + tmp1(2)**2 + tmp1(3)**2
            tmp22 = tmp2(1)**2 + tmp2(2)**2 + tmp2(3)**2

            o_tmp12 = 1.0d0 / tmp12
            o_tmp22 = 1.0d0 / tmp22

            o_tmp1vtmp2v = sqrt( o_tmp12 * o_tmp22 )

            arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3)) * o_tmp1vtmp2v

            asc = (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*zaxis(1) + &
                  (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*zaxis(2) + &
                  (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*zaxis(3)

            if ( arg .gt.  1.0 ) then
                arg =  1.0
            else if ( arg .lt. -1.0 ) then
                arg = -1.0
            end if

            ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option in calculate_nasstp!')
    end select

! derivatives ====================================

    select case(cv_item%lst_par)
        case(1,2,3)
            ! nothing to be here
        case(4,5,6)
            argo_tmp12 = arg*o_tmp12;
            argo_tmp22 = arg*o_tmp22;

            f1 = sin(ctx%CVsValues(cv_item%idx))
            if( abs(f1) .lt. 1.e-12 ) then
                ! avoid division by zero
                f1 = -1.e12
            else
                f1 = -1.0d0 / f1
            end if

        ! with respect to tmp1 and tmp2
            a_tmp1(:) = f1*(tmp2(:)*o_tmp1vtmp2v - tmp1(:)*argo_tmp12)
            a_tmp2(:) = f1*(tmp1(:)*o_tmp1vtmp2v - tmp2(:)*argo_tmp22)

        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option in calculate_nasbpp!')
    end select

    a_xaxis(:)  = 0.0d0
    a_yaxis(:)  = 0.0d0
    a_zaxis(:)  = 0.0d0
! final derivatives
    a_ua(:,:)   = 0.0d0
    a_ub(:,:)   = 0.0d0
    a_oa(:)     = 0.0d0
    a_ob(:)     = 0.0d0

    select case(cv_item%lst_par)
!---------------------------------------------------------------------------------------------------
        case(1)
            ! 'shift'
        ! with respect to xaxis and d
            a_xaxis(:) = d(:)
            a_oa(:) = - xaxis(:)
            a_ob(:) =   xaxis(:)

! with respect to xaxis
            t1 = xaxisr(1)*a_xaxis(1) + xaxisr(2)*a_xaxis(2) + xaxisr(3)*a_xaxis(3)
            a_xaxisr(:) =   a_xaxis(:)*o_xaxis - o_xaxis*o_xaxis2*xaxisr(:)*t1
            t1 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_x0axis(:) = a_xaxisr(:) - zaxis(:)*t1
            a_ua(:,1) = a_x0axis(:)
            a_ub(:,1) = xsc*a_x0axis(:)

! with respect to zaxis
            t1 = x0axis(1)*zaxis(1)+x0axis(2)*zaxis(2)+x0axis(3)*zaxis(3)
            t2 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_zaxis(:) = - t1*a_xaxisr(:) - x0axis(:)*t2

            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

!---------------------------------------------------------------------------------------------------
        case(2)
            ! 'slide'
!           ctx%CVsValues(cv_item%idx) = d(1)*yaxis(1) + d(2)*yaxis(2) + d(3)*yaxis(3)
        ! with respect to yaxis and d
            a_yaxis(:) = d(:)
            a_oa(:) = - yaxis(:)
            a_ob(:) =   yaxis(:)

!    ! is cross product of z and x vectors
!    yaxis(1) = zaxis(2)*xaxis(3) - zaxis(3)*xaxis(2)
!    yaxis(2) = zaxis(3)*xaxis(1) - zaxis(1)*xaxis(3)
!    yaxis(3) = zaxis(1)*xaxis(2) - zaxis(2)*xaxis(1)

        ! with respect to xaxis
            a_xaxis(1) = - zaxis(2)*a_yaxis(3) + zaxis(3)*a_yaxis(2)
            a_xaxis(2) = - zaxis(3)*a_yaxis(1) + zaxis(1)*a_yaxis(3)
            a_xaxis(3) = - zaxis(1)*a_yaxis(2) + zaxis(2)*a_yaxis(1)

        ! with respect to xaxis
            t1 = xaxisr(1)*a_xaxis(1) + xaxisr(2)*a_xaxis(2) + xaxisr(3)*a_xaxis(3)
            a_xaxisr(:) =   a_xaxis(:)*o_xaxis - o_xaxis*o_xaxis2*xaxisr(:)*t1
            t1 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_x0axis(:) = a_xaxisr(:) - zaxis(:)*t1
            a_ua(:,1) = a_x0axis(:)
            a_ub(:,1) = xsc*a_x0axis(:)

        ! with respect to zaxis
            t1 = x0axis(1)*zaxis(1)+x0axis(2)*zaxis(2)+x0axis(3)*zaxis(3)
            t2 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_zaxis(:) = - t1*a_xaxisr(:) - x0axis(:)*t2

            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

        ! with respect to zaxis
            a_zaxis(1) = - a_yaxis(2)*xaxis(3) + a_yaxis(3)*xaxis(2)
            a_zaxis(2) = - a_yaxis(3)*xaxis(1) + a_yaxis(1)*xaxis(3)
            a_zaxis(3) = - a_yaxis(1)*xaxis(2) + a_yaxis(2)*xaxis(1)

        ! with respect to zaxis
            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_ua(:,3) + a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = a_ub(:,3) + zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

!---------------------------------------------------------------------------------------------------
        case(3)
            ! 'rise'
        ! with respect to zaxis and d
            a_zaxis(:) = d(:)
            a_oa(:) = - zaxis(:)
            a_ob(:) =   zaxis(:)

        ! with respect to zaxis
            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

!---------------------------------------------------------------------------------------------------
        case(4)
            ! 'tilt'
        ! with respect to ua, ub, and xaxis
            a_ua(1,3) = xaxis(2)*a_tmp1(3) - a_tmp1(2)*xaxis(3)
            a_ua(2,3) = xaxis(3)*a_tmp1(1) - a_tmp1(3)*xaxis(1)
            a_ua(3,3) = xaxis(1)*a_tmp1(2) - a_tmp1(1)*xaxis(2)

            a_ub(1,3) = zsc*xaxis(2)*a_tmp2(3) - zsc*a_tmp2(2)*xaxis(3)
            a_ub(2,3) = zsc*xaxis(3)*a_tmp2(1) - zsc*a_tmp2(3)*xaxis(1)
            a_ub(3,3) = zsc*xaxis(1)*a_tmp2(2) - zsc*a_tmp2(1)*xaxis(2)

            a_xaxis(1) = a_tmp1(2)*ua(3,3) - ua(2,3)*a_tmp1(3) + zsc*a_tmp2(2)*ub(3,3) - zsc*ub(2,3)*a_tmp2(3)
            a_xaxis(2) = a_tmp1(3)*ua(1,3) - ua(3,3)*a_tmp1(1) + zsc*a_tmp2(3)*ub(1,3) - zsc*ub(3,3)*a_tmp2(1)
            a_xaxis(3) = a_tmp1(1)*ua(2,3) - ua(1,3)*a_tmp1(2) + zsc*a_tmp2(1)*ub(2,3) - zsc*ub(1,3)*a_tmp2(2)

! with respect to xaxis
            t1 = xaxisr(1)*a_xaxis(1) + xaxisr(2)*a_xaxis(2) + xaxisr(3)*a_xaxis(3)
            a_xaxisr(:) =   a_xaxis(:)*o_xaxis - o_xaxis*o_xaxis2*xaxisr(:)*t1
            t1 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_x0axis(:) = a_xaxisr(:) - zaxis(:)*t1
            a_ua(:,1) = a_x0axis(:)
            a_ub(:,1) = xsc*a_x0axis(:)

! with respect to zaxis
            t1 = x0axis(1)*zaxis(1)+x0axis(2)*zaxis(2)+x0axis(3)*zaxis(3)
            t2 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_zaxis(:) = - t1*a_xaxisr(:) - x0axis(:)*t2

            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_ua(:,3) + a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = a_ub(:,3) + zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

!---------------------------------------------------------------------------------------------------
        case(5)
            ! 'roll'
        ! with respect to ua, ub, and yaxis
            a_ua(1,3) = yaxis(2)*a_tmp1(3) - a_tmp1(2)*yaxis(3)
            a_ua(2,3) = yaxis(3)*a_tmp1(1) - a_tmp1(3)*yaxis(1)
            a_ua(3,3) = yaxis(1)*a_tmp1(2) - a_tmp1(1)*yaxis(2)

            a_ub(1,3) = zsc*yaxis(2)*a_tmp2(3) - zsc*a_tmp2(2)*yaxis(3)
            a_ub(2,3) = zsc*yaxis(3)*a_tmp2(1) - zsc*a_tmp2(3)*yaxis(1)
            a_ub(3,3) = zsc*yaxis(1)*a_tmp2(2) - zsc*a_tmp2(1)*yaxis(2)

            a_yaxis(1) = a_tmp1(2)*ua(3,3) - ua(2,3)*a_tmp1(3) + zsc*a_tmp2(2)*ub(3,3) - zsc*ub(2,3)*a_tmp2(3)
            a_yaxis(2) = a_tmp1(3)*ua(1,3) - ua(3,3)*a_tmp1(1) + zsc*a_tmp2(3)*ub(1,3) - zsc*ub(3,3)*a_tmp2(1)
            a_yaxis(3) = a_tmp1(1)*ua(2,3) - ua(1,3)*a_tmp1(2) + zsc*a_tmp2(1)*ub(2,3) - zsc*ub(1,3)*a_tmp2(2)

!    ! is cross product of z and x vectors
!    yaxis(1) = zaxis(2)*xaxis(3) - zaxis(3)*xaxis(2)
!    yaxis(2) = zaxis(3)*xaxis(1) - zaxis(1)*xaxis(3)
!    yaxis(3) = zaxis(1)*xaxis(2) - zaxis(2)*xaxis(1)

        ! with respect to xaxis
            a_xaxis(1) = - zaxis(2)*a_yaxis(3) + zaxis(3)*a_yaxis(2)
            a_xaxis(2) = - zaxis(3)*a_yaxis(1) + zaxis(1)*a_yaxis(3)
            a_xaxis(3) = - zaxis(1)*a_yaxis(2) + zaxis(2)*a_yaxis(1)

        ! with respect to xaxis
            t1 = xaxisr(1)*a_xaxis(1) + xaxisr(2)*a_xaxis(2) + xaxisr(3)*a_xaxis(3)
            a_xaxisr(:) =   a_xaxis(:)*o_xaxis - o_xaxis*o_xaxis2*xaxisr(:)*t1
            t1 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_x0axis(:) = a_xaxisr(:) - zaxis(:)*t1
            a_ua(:,1) = a_x0axis(:)
            a_ub(:,1) = xsc*a_x0axis(:)

        ! with respect to zaxis
            t1 = x0axis(1)*zaxis(1)+x0axis(2)*zaxis(2)+x0axis(3)*zaxis(3)
            t2 = zaxis(1)*a_xaxisr(1) + zaxis(2)*a_xaxisr(2) + zaxis(3)*a_xaxisr(3)
            a_zaxis(:) = - t1*a_xaxisr(:) - x0axis(:)*t2

            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_ua(:,3) + a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = a_ub(:,3) + zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

        ! with respect to zaxis
            a_zaxis(1) = - a_yaxis(2)*xaxis(3) + a_yaxis(3)*xaxis(2)
            a_zaxis(2) = - a_yaxis(3)*xaxis(1) + a_yaxis(1)*xaxis(3)
            a_zaxis(3) = - a_yaxis(1)*xaxis(2) + a_yaxis(2)*xaxis(1)

        ! with respect to zaxis
            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_ua(:,3) + a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = a_ub(:,3) + zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

!---------------------------------------------------------------------------------------------------
        case(6)
            ! twist
        ! with respect to ua, ub, and zaxis
            a_ua(1,2) = zaxis(2)*a_tmp1(3) - a_tmp1(2)*zaxis(3)
            a_ua(2,2) = zaxis(3)*a_tmp1(1) - a_tmp1(3)*zaxis(1)
            a_ua(3,2) = zaxis(1)*a_tmp1(2) - a_tmp1(1)*zaxis(2)

            a_ub(1,2) = zsc*zaxis(2)*a_tmp2(3) - zsc*a_tmp2(2)*zaxis(3)
            a_ub(2,2) = zsc*zaxis(3)*a_tmp2(1) - zsc*a_tmp2(3)*zaxis(1)
            a_ub(3,2) = zsc*zaxis(1)*a_tmp2(2) - zsc*a_tmp2(1)*zaxis(2)

            a_zaxis(1) = a_tmp1(2)*ua(3,2) - ua(2,2)*a_tmp1(3) + zsc*a_tmp2(2)*ub(3,2) - zsc*ub(2,2)*a_tmp2(3)
            a_zaxis(2) = a_tmp1(3)*ua(1,2) - ua(3,2)*a_tmp1(1) + zsc*a_tmp2(3)*ub(1,2) - zsc*ub(3,2)*a_tmp2(1)
            a_zaxis(3) = a_tmp1(1)*ua(2,2) - ua(1,2)*a_tmp1(2) + zsc*a_tmp2(1)*ub(2,2) - zsc*ub(1,2)*a_tmp2(2)

        ! with respect to zaxis
            t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
            a_ua(:,3) = a_zaxis(:)*o_zaxis - o_zaxis*o_zaxis2*zaxisr(:)*t1
            a_ub(:,3) = zsc*a_zaxis(:)*o_zaxis - zsc*o_zaxis*o_zaxis2*zaxisr(:)*t1

        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option in calculate_nasstp!')
    end select

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
    integer             :: i,ai,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: o_zaxis2,o_zaxis
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zaxisr(3),zsc,y0axis(3)
    real(PMFDP)         :: tmp1(3),ingra,ingrb,oa(3),ob(3)
    real(PMFDP)         :: yaxisr(3),o_yaxis,o_yaxis2
    ! --------------------------------------------------------------------------

    ! inverse number of atoms
    ingra = 1.0d0 / cv_item%grps(1)
    ingrb = 1.0d0 / (cv_item%grps(2)-cv_item%grps(1))

    ! calculate geometrical centres (source and target) -------------------
    xsa(:) = 0.0d0
    xra(:) = 0.0d0

    do  i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        ! source
        xsa(:) = xsa(:) + x(:,ai)

        ! reference
        xra(:) = xra(:) + cv_item%xyz_str_a1%cvs(:,i)
    end do

    xsa(:) = xsa(:) * ingra
    xra(:) = xra(:) * ingra

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))
    end do

    r11 = r11 * ingra
    r12 = r12 * ingra
    r13 = r13 * ingra

    r21 = r21 * ingra
    r22 = r22 * ingra
    r23 = r23 * ingra

    r31 = r31 * ingra
    r32 = r32 * ingra
    r33 = r33 * ingra

    ! construct matrix for quaterion fitting
    fa(1,1) =  r11 + r22 + r33
    fa(1,2) =  r23 - r32
    fa(1,3) =  r31 - r13
    fa(1,4) =  r12 - r21

    fa(2,1) =  r23 - r32
    fa(2,2) =  r11 - r22 - r33
    fa(2,3) =  r12 + r21
    fa(2,4) =  r13 + r31

    fa(3,1) =  r31 - r13
    fa(3,2) =  r12 + r21
    fa(3,3) = -r11 + r22 - r33
    fa(3,4) =  r23 + r32

    fa(4,1) =  r12 - r21
    fa(4,2) =  r13 + r31
    fa(4,3) =  r23 + r32
    fa(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesa(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fa, 4, eigenvaluesa, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point a!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    xsb(:) = 0.0d0
    xrb(:) = 0.0d0

    do  i = cv_item%grps(1)+1,cv_item%grps(2)
        ai = cv_item%lindexes(i)

        ! source
        xsb(:) = xsb(:) + x(:,ai)

        ! reference
        xrb(:) = xrb(:) + cv_item%xyz_str_b1%cvs(:,i-cv_item%grps(1))
    end do

    xsb(:) = xsb(:) * ingrb
    xrb(:) = xrb(:) * ingrb

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = cv_item%grps(1)+1,cv_item%grps(2)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))
    end do

    r11 = r11 * ingrb
    r12 = r12 * ingrb
    r13 = r13 * ingrb

    r21 = r21 * ingrb
    r22 = r22 * ingrb
    r23 = r23 * ingrb

    r31 = r31 * ingrb
    r32 = r32 * ingrb
    r33 = r33 * ingrb

    ! construct matrix for quaterion fitting
    fb(1,1) =  r11 + r22 + r33
    fb(1,2) =  r23 - r32
    fb(1,3) =  r31 - r13
    fb(1,4) =  r12 - r21

    fb(2,1) =  r23 - r32
    fb(2,2) =  r11 - r22 - r33
    fb(2,3) =  r12 + r21
    fb(2,4) =  r13 + r31

    fb(3,1) =  r31 - r13
    fb(3,2) =  r12 + r21
    fb(3,3) = -r11 + r22 - r33
    fb(3,4) =  r23 + r32

    fb(4,1) =  r12 - r21
    fb(4,2) =  r13 + r31
    fb(4,3) =  r23 + r32
    fb(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesb(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fb, 4, eigenvaluesb, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
    ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
    ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
    ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
    ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
    ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
    ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    ! rotation matrix b  ------------------------------
    ub(1,1) = fb(1,best)**2 + fb(2,best)**2 - fb(3,best)**2 - fb(4,best)**2
    ub(2,1) = 2.0d0*( fb(2,best)*fb(3,best) - fb(1,best)*fb(4,best) )
    ub(3,1) = 2.0d0*( fb(2,best)*fb(4,best) + fb(1,best)*fb(3,best) )

    ub(1,2) = 2.0d0*( fb(2,best)*fb(3,best) + fb(1,best)*fb(4,best) )
    ub(2,2) = fb(1,best)**2 - fb(2,best)**2 + fb(3,best)**2 - fb(4,best)**2
    ub(3,2) = 2.0d0*( fb(3,best)*fb(4,best) - fb(1,best)*fb(2,best) )

    ub(1,3) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
    ub(2,3) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
    ub(3,3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2

    ! z-axis ===================================================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_yaxis2 = 1.0d0 / (yaxisr(1)**2 + yaxisr(2)**2 + yaxisr(3)**2)
    o_yaxis  = sqrt(o_yaxis2)
    yaxis(:) = yaxisr(:) * o_yaxis

    ! x-axis ===================================================================
    ! is cross product of y and z axes
    xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
    xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
    xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)

    ! get origins of bases
    oa(:) = 0.0d0
    ! move reference point to origin
    oa(:) = oa(:) - xra(:)
    ! rotate
    tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
    tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
    tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
    ! move origin to new reference point (experimental structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

    ! position of bp origin
    morg(:) = 0.5d0*(oa(:) + ob(:))

    ! axis
    mst(:,1) = xaxis(:)
    mst(:,2) = yaxis(:)
    mst(:,3) = zaxis(:)

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
    integer             :: i,ai,info,best,mi,mj
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: o_zaxis2,o_zaxis
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zaxisr(3),zsc,y0axis(3)
    real(PMFDP)         :: tmp1(3),ingra,ingrb,oa(3),ob(3)
    real(PMFDP)         :: yaxisr(3),o_yaxis,o_yaxis2

    real(PMFDP)         :: a_fa(4),a_fb(4),a_zaxis(3),a_rij(4,4),a_xaxis(3),a_y0axis(3)
    real(PMFDP)         :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4),a_xsa(3),a_xsb(3)
    real(PMFDP)         :: a_yaxisr(3),a_ua(3,3),a_ub(3,3),a_yaxis(3)
    real(PMFDP)         :: t1,t2
    ! --------------------------------------------------------------------------

! inverse number of atoms
    ingra = 1.0d0 / cv_item%grps(1)
    ingrb = 1.0d0 / (cv_item%grps(2)-cv_item%grps(1))

    ! calculate geometrical centres (source and target) -------------------
    xsa(:) = 0.0d0
    xra(:) = 0.0d0

    do  i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        ! source
        xsa(:) = xsa(:) + x(:,ai)

        ! reference
        xra(:) = xra(:) + cv_item%xyz_str_a1%cvs(:,i)
    end do

    xsa(:) = xsa(:) * ingra
    xra(:) = xra(:) * ingra

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))
    end do

    r11 = r11 * ingra
    r12 = r12 * ingra
    r13 = r13 * ingra

    r21 = r21 * ingra
    r22 = r22 * ingra
    r23 = r23 * ingra

    r31 = r31 * ingra
    r32 = r32 * ingra
    r33 = r33 * ingra

    ! construct matrix for quaterion fitting
    fa(1,1) =  r11 + r22 + r33
    fa(1,2) =  r23 - r32
    fa(1,3) =  r31 - r13
    fa(1,4) =  r12 - r21

    fa(2,1) =  r23 - r32
    fa(2,2) =  r11 - r22 - r33
    fa(2,3) =  r12 + r21
    fa(2,4) =  r13 + r31

    fa(3,1) =  r31 - r13
    fa(3,2) =  r12 + r21
    fa(3,3) = -r11 + r22 - r33
    fa(3,4) =  r23 + r32

    fa(4,1) =  r12 - r21
    fa(4,2) =  r13 + r31
    fa(4,3) =  r23 + r32
    fa(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesa(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fa, 4, eigenvaluesa, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point a!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    xsb(:) = 0.0d0
    xrb(:) = 0.0d0

    do  i = cv_item%grps(1)+1,cv_item%grps(2)
        ai = cv_item%lindexes(i)

        ! source
        xsb(:) = xsb(:) + x(:,ai)

        ! reference
        xrb(:) = xrb(:) + cv_item%xyz_str_b1%cvs(:,i-cv_item%grps(1))
    end do

    xsb(:) = xsb(:) * ingrb
    xrb(:) = xrb(:) * ingrb

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = cv_item%grps(1)+1,cv_item%grps(2)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))
    end do

    r11 = r11 * ingrb
    r12 = r12 * ingrb
    r13 = r13 * ingrb

    r21 = r21 * ingrb
    r22 = r22 * ingrb
    r23 = r23 * ingrb

    r31 = r31 * ingrb
    r32 = r32 * ingrb
    r33 = r33 * ingrb

    ! construct matrix for quaterion fitting
    fb(1,1) =  r11 + r22 + r33
    fb(1,2) =  r23 - r32
    fb(1,3) =  r31 - r13
    fb(1,4) =  r12 - r21

    fb(2,1) =  r23 - r32
    fb(2,2) =  r11 - r22 - r33
    fb(2,3) =  r12 + r21
    fb(2,4) =  r13 + r31

    fb(3,1) =  r31 - r13
    fb(3,2) =  r12 + r21
    fb(3,3) = -r11 + r22 - r33
    fb(3,4) =  r23 + r32

    fb(4,1) =  r12 - r21
    fb(4,2) =  r13 + r31
    fb(4,3) =  r23 + r32
    fb(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesb(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fb, 4, eigenvaluesb, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
    ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
    ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
    ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
    ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
    ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
    ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    ! rotation matrix b  ------------------------------
    ub(1,1) = fb(1,best)**2 + fb(2,best)**2 - fb(3,best)**2 - fb(4,best)**2
    ub(2,1) = 2.0d0*( fb(2,best)*fb(3,best) - fb(1,best)*fb(4,best) )
    ub(3,1) = 2.0d0*( fb(2,best)*fb(4,best) + fb(1,best)*fb(3,best) )

    ub(1,2) = 2.0d0*( fb(2,best)*fb(3,best) + fb(1,best)*fb(4,best) )
    ub(2,2) = fb(1,best)**2 - fb(2,best)**2 + fb(3,best)**2 - fb(4,best)**2
    ub(3,2) = 2.0d0*( fb(3,best)*fb(4,best) - fb(1,best)*fb(2,best) )

    ub(1,3) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
    ub(2,3) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
    ub(3,3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2

    ! z-axis ===================================================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_yaxis2 = 1.0d0 / (yaxisr(1)**2 + yaxisr(2)**2 + yaxisr(3)**2)
    o_yaxis  = sqrt(o_yaxis2)
    yaxis(:) = yaxisr(:) * o_yaxis

    ! x-axis ===================================================================
    ! is cross product of y and z axes
    xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
    xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
    xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)

    ! get origins of bases
    oa(:) = 0.0d0
    ! move reference point to origin
    oa(:) = oa(:) - xra(:)
    ! rotate
    tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
    tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
    tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
    ! move origin to new reference point (experimental structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

! ==============================================================================

! final derivatives
    a_ua(:,:)   = 0.0d0
    a_ub(:,:)   = 0.0d0
    a_y0axis(:) = 0.0d0
    a_xsa(:)    = 0.0d0
    a_xsb(:)    = 0.0d0

! ************
! ==== a_xaxis
! ************

    a_xaxis(:) = a_mst(:,1)

! a_xaxis with respect to ua and ub
! with respect to yaxis
!     xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
!     xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
!     xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)
    a_yaxis(1) = - a_xaxis(2)*zaxis(3) + zaxis(2)*a_xaxis(3)
    a_yaxis(2) = - a_xaxis(3)*zaxis(1) + zaxis(3)*a_xaxis(1)
    a_yaxis(3) = - a_xaxis(1)*zaxis(2) + zaxis(1)*a_xaxis(2)

! with respect to yaxisr
    t1 = yaxisr(1)*a_yaxis(1) + yaxisr(2)*a_yaxis(2) + yaxisr(3)*a_yaxis(3)
    a_yaxisr(1) =   a_yaxis(1)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(1)*t1
    a_yaxisr(2) =   a_yaxis(2)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(2)*t1
    a_yaxisr(3) =   a_yaxis(3)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(3)*t1

! with respect to y0axis
    t1 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_y0axis(1) = a_yaxisr(1) - zaxis(1)*t1
    a_y0axis(2) = a_yaxisr(2) - zaxis(2)*t1
    a_y0axis(3) = a_yaxisr(3) - zaxis(3)*t1

! with respect to zaxis
    t1 = y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3)
    t2 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_zaxis(1) = - t1*a_yaxisr(1) - y0axis(1)*t2
    a_zaxis(2) = - t1*a_yaxisr(2) - y0axis(2)*t2
    a_zaxis(3) = - t1*a_yaxisr(3) - y0axis(3)*t2

! with respect to zaxis
!     xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
!     xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
!     xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)
    a_zaxis(1) = a_zaxis(1) - yaxis(2)*a_xaxis(3) + yaxis(3)*a_xaxis(2)
    a_zaxis(2) = a_zaxis(2) - yaxis(3)*a_xaxis(1) + yaxis(1)*a_xaxis(3)
    a_zaxis(3) = a_zaxis(3) - yaxis(1)*a_xaxis(2) + yaxis(2)*a_xaxis(1)

    t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
    a_ua(1,3) = a_ua(1,3) + 0.5d0*a_zaxis(1)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ua(2,3) = a_ua(2,3) + 0.5d0*a_zaxis(2)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ua(3,3) = a_ua(3,3) + 0.5d0*a_zaxis(3)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(3)*t1

    a_ub(1,3) = a_ub(1,3) + 0.5d0*zsc*a_zaxis(1)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ub(2,3) = a_ub(2,3) + 0.5d0*zsc*a_zaxis(2)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ub(3,3) = a_ub(3,3) + 0.5d0*zsc*a_zaxis(3)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(3)*t1


! ************
! ==== a_yaxis
! ************

    a_yaxis(:) = a_mst(:,2)

! with respect to yaxisr
    t1 = yaxisr(1)*a_yaxis(1) + yaxisr(2)*a_yaxis(2) + yaxisr(3)*a_yaxis(3)
    a_yaxisr(1) =   a_yaxis(1)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(1)*t1
    a_yaxisr(2) =   a_yaxis(2)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(2)*t1
    a_yaxisr(3) =   a_yaxis(3)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(3)*t1
! with respect to y0axis
    t1 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_y0axis(1) = a_y0axis(1) + a_yaxisr(1) - zaxis(1)*t1
    a_y0axis(2) = a_y0axis(2) + a_yaxisr(2) - zaxis(2)*t1
    a_y0axis(3) = a_y0axis(3) + a_yaxisr(3) - zaxis(3)*t1

! with respect to zaxis
    t1 = y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3)
    t2 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_zaxis(1) = - t1*a_yaxisr(1) - y0axis(1)*t2
    a_zaxis(2) = - t1*a_yaxisr(2) - y0axis(2)*t2
    a_zaxis(3) = - t1*a_yaxisr(3) - y0axis(3)*t2

    t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
    a_ua(1,3) = a_ua(1,3) + 0.5d0*a_zaxis(1)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ua(2,3) = a_ua(2,3) + 0.5d0*a_zaxis(2)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ua(3,3) = a_ua(3,3) + 0.5d0*a_zaxis(3)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(3)*t1

    a_ub(1,3) = a_ub(1,3) + 0.5d0*zsc*a_zaxis(1)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ub(2,3) = a_ub(2,3) + 0.5d0*zsc*a_zaxis(2)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ub(3,3) = a_ub(3,3) + 0.5d0*zsc*a_zaxis(3)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(3)*t1

! ************
! ==== a_zaxis
! ************

    a_zaxis(:) = a_mst(:,3)

! with respect to zaxis
    t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
    a_ua(1,3) = a_ua(1,3) + 0.5d0*a_zaxis(1)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ua(2,3) = a_ua(2,3) + 0.5d0*a_zaxis(2)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ua(3,3) = a_ua(3,3) + 0.5d0*a_zaxis(3)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(3)*t1

    a_ub(1,3) = a_ub(1,3) + 0.5d0*zsc*a_zaxis(1)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ub(2,3) = a_ub(2,3) + 0.5d0*zsc*a_zaxis(2)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ub(3,3) = a_ub(3,3) + 0.5d0*zsc*a_zaxis(3)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(3)*t1

! ************
! ==== origin
! ************

!     ! get origins of bases
!     oa(:) = 0.0d0
!     ! move reference point to origin
!     oa(:) = oa(:) - xra(:)
!     ! rotate
!     tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
!     tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
!     tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
!     ! move origin to new reference point (experiemntal structure)
!     oa(:) = tmp1(:) + xsa(:)
!
!     ob(:) = 0.0d0
!     ob(:) = ob(:) - xrb(:)
!     tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
!     tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
!     tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
!     ob(:) = tmp1(:) + xsb(:)

! a_morg
!    ! position of bp origin
!    morg(:) = 0.5d0*(oa(:) + ob(:))

! a_morg with respect to ua and ub
    a_ua(1,1) = a_ua(1,1) - 0.5d0*xra(1)*a_morg(1)
    a_ua(2,1) = a_ua(2,1) - 0.5d0*xra(1)*a_morg(2)
    a_ua(3,1) = a_ua(3,1) - 0.5d0*xra(1)*a_morg(3)
    a_ua(1,2) = a_ua(1,2) - 0.5d0*xra(2)*a_morg(1)
    a_ua(2,2) = a_ua(2,2) - 0.5d0*xra(2)*a_morg(2)
    a_ua(3,2) = a_ua(3,2) - 0.5d0*xra(2)*a_morg(3)
    a_ua(1,3) = a_ua(1,3) - 0.5d0*xra(3)*a_morg(1)
    a_ua(2,3) = a_ua(2,3) - 0.5d0*xra(3)*a_morg(2)
    a_ua(3,3) = a_ua(3,3) - 0.5d0*xra(3)*a_morg(3)

    a_ub(1,1) = a_ub(1,1) - 0.5d0*xrb(1)*a_morg(1)
    a_ub(2,1) = a_ub(2,1) - 0.5d0*xrb(1)*a_morg(2)
    a_ub(3,1) = a_ub(3,1) - 0.5d0*xrb(1)*a_morg(3)
    a_ub(1,2) = a_ub(1,2) - 0.5d0*xrb(2)*a_morg(1)
    a_ub(2,2) = a_ub(2,2) - 0.5d0*xrb(2)*a_morg(2)
    a_ub(3,2) = a_ub(3,2) - 0.5d0*xrb(2)*a_morg(3)
    a_ub(1,3) = a_ub(1,3) - 0.5d0*xrb(3)*a_morg(1)
    a_ub(2,3) = a_ub(2,3) - 0.5d0*xrb(3)*a_morg(2)
    a_ub(3,3) = a_ub(3,3) - 0.5d0*xrb(3)*a_morg(3)

    ! a_d with respect to xsa, xsb
    a_xsa(:) = + 0.5d0*ingra*a_morg(:)
    a_xsb(:) = + 0.5d0*ingrb*a_morg(:)

! rotation matrix a ------------------------------
!     ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
!     ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
!     ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    a_fa(1) = 2.0d0*( fa(1,best)*a_ua(1,1) - fa(4,best)*a_ua(2,1) + fa(3,best)*a_ua(3,1))
    a_fa(2) = 2.0d0*( fa(2,best)*a_ua(1,1) + fa(3,best)*a_ua(2,1) + fa(4,best)*a_ua(3,1))
    a_fa(3) = 2.0d0*(-fa(3,best)*a_ua(1,1) + fa(2,best)*a_ua(2,1) + fa(1,best)*a_ua(3,1))
    a_fa(4) = 2.0d0*(-fa(4,best)*a_ua(1,1) - fa(1,best)*a_ua(2,1) + fa(2,best)*a_ua(3,1))

    a_fb(1) = 2.0d0*( fb(1,best)*a_ub(1,1) - fb(4,best)*a_ub(2,1) + fb(3,best)*a_ub(3,1))
    a_fb(2) = 2.0d0*( fb(2,best)*a_ub(1,1) + fb(3,best)*a_ub(2,1) + fb(4,best)*a_ub(3,1))
    a_fb(3) = 2.0d0*(-fb(3,best)*a_ub(1,1) + fb(2,best)*a_ub(2,1) + fb(1,best)*a_ub(3,1))
    a_fb(4) = 2.0d0*(-fb(4,best)*a_ub(1,1) - fb(1,best)*a_ub(2,1) + fb(2,best)*a_ub(3,1))

!     ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
!     ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
!     ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    a_fa(1) = a_fa(1) + 2.0d0*(fa(4,best)*a_ua(1,2) + fa(1,best)*a_ua(2,2) - fa(2,best)*a_ua(3,2))
    a_fa(2) = a_fa(2) + 2.0d0*(fa(3,best)*a_ua(1,2) - fa(2,best)*a_ua(2,2) - fa(1,best)*a_ua(3,2))
    a_fa(3) = a_fa(3) + 2.0d0*(fa(2,best)*a_ua(1,2) + fa(3,best)*a_ua(2,2) + fa(4,best)*a_ua(3,2))
    a_fa(4) = a_fa(4) + 2.0d0*(fa(1,best)*a_ua(1,2) - fa(4,best)*a_ua(2,2) + fa(3,best)*a_ua(3,2))

    a_fb(1) = a_fb(1) + 2.0d0*(fb(4,best)*a_ub(1,2) + fb(1,best)*a_ub(2,2) - fb(2,best)*a_ub(3,2))
    a_fb(2) = a_fb(2) + 2.0d0*(fb(3,best)*a_ub(1,2) - fb(2,best)*a_ub(2,2) - fb(1,best)*a_ub(3,2))
    a_fb(3) = a_fb(3) + 2.0d0*(fb(2,best)*a_ub(1,2) + fb(3,best)*a_ub(2,2) + fb(4,best)*a_ub(3,2))
    a_fb(4) = a_fb(4) + 2.0d0*(fb(1,best)*a_ub(1,2) - fb(4,best)*a_ub(2,2) + fb(3,best)*a_ub(3,2))

!     ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
!     ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
!     ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    a_fa(1) = a_fa(1) + 2.0d0*(-fa(3,best)*a_ua(1,3) + fa(2,best)*a_ua(2,3) + fa(1,best)*a_ua(3,3))
    a_fa(2) = a_fa(2) + 2.0d0*( fa(4,best)*a_ua(1,3) + fa(1,best)*a_ua(2,3) - fa(2,best)*a_ua(3,3))
    a_fa(3) = a_fa(3) + 2.0d0*(-fa(1,best)*a_ua(1,3) + fa(4,best)*a_ua(2,3) - fa(3,best)*a_ua(3,3))
    a_fa(4) = a_fa(4) + 2.0d0*( fa(2,best)*a_ua(1,3) + fa(3,best)*a_ua(2,3) + fa(4,best)*a_ua(3,3))

    a_fb(1) = a_fb(1) + 2.0d0*(-fb(3,best)*a_ub(1,3) + fb(2,best)*a_ub(2,3) + fb(1,best)*a_ub(3,3))
    a_fb(2) = a_fb(2) + 2.0d0*( fb(4,best)*a_ub(1,3) + fb(1,best)*a_ub(2,3) - fb(2,best)*a_ub(3,3))
    a_fb(3) = a_fb(3) + 2.0d0*(-fb(1,best)*a_ub(1,3) + fb(4,best)*a_ub(2,3) - fb(3,best)*a_ub(3,3))
    a_fb(4) = a_fb(4) + 2.0d0*( fb(2,best)*a_ub(1,3) + fb(3,best)*a_ub(2,3) + fb(4,best)*a_ub(3,3))

! derivatives of fa with respect to matrix elements
    v(:,:) = fa(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvaluesa(i) - eigenvaluesa(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + fa(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

! merge xij with a_fa, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fa(1)*xij(1,mi,mj)+a_fa(2)*xij(2,mi,mj)+a_fa(3)*xij(3,mi,mj)+a_fa(4)*xij(4,mi,mj))*ingra
        end do
    end do

! finaly gradients for group_a
    do i = 1, cv_item%grps(1)

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3)) &
                + a_xsa(1)

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3)) &
                + a_xsa(2)

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3)) &
                + a_xsa(3)

    end do

! derivatives of fb with respect to matrix elements
    v(:,:) = fb(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvaluesb(i) - eigenvaluesb(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + fb(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

! merge xij with a_fb, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fb(1)*xij(1,mi,mj)+a_fb(2)*xij(2,mi,mj)+a_fb(3)*xij(3,mi,mj)+a_fb(4)*xij(4,mi,mj))*ingrb
        end do
    end do

! finaly gradients for group_b
    do i = cv_item%grps(1) + 1, cv_item%grps(2)

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3)) &
                + a_xsb(1)

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3)) &
                + a_xsb(2)

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3)) &
                + a_xsb(3)
    end do

! finally gradients for group_c, group_d
    ai = cv_item%lindexes(cv_item%grps(3))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y0axis(:)

    ai = cv_item%lindexes(cv_item%grps(4))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - a_y0axis(:)

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
    integer             :: i,ai,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: o_zaxis2,o_zaxis
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zaxisr(3),zsc,y0axis(3)
    real(PMFDP)         :: tmp1(3),ingra,ingrb,oa(3),ob(3)
    real(PMFDP)         :: yaxisr(3),o_yaxis,o_yaxis2
    ! --------------------------------------------------------------------------

    ! inverse number of atoms
    ingra = 1.0d0 / (cv_item%grps(5)-cv_item%grps(4))
    ingrb = 1.0d0 / (cv_item%grps(6)-cv_item%grps(5))

    ! calculate geometrical centres (source and target) -------------------
    xsa(:) = 0.0d0
    xra(:) = 0.0d0

    do  i = cv_item%grps(4)+1,cv_item%grps(5)
        ai = cv_item%lindexes(i)
        ! source
        xsa(:) = xsa(:) + x(:,ai)

        ! reference
        xra(:) = xra(:) + cv_item%xyz_str_a2%cvs(:,i-cv_item%grps(4))
    end do

    xsa(:) = xsa(:) * ingra
    xra(:) = xra(:) * ingra

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = cv_item%grps(4)+1,cv_item%grps(5)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3))
    end do

    r11 = r11 * ingra
    r12 = r12 * ingra
    r13 = r13 * ingra

    r21 = r21 * ingra
    r22 = r22 * ingra
    r23 = r23 * ingra

    r31 = r31 * ingra
    r32 = r32 * ingra
    r33 = r33 * ingra

    ! construct matrix for quaterion fitting
    fa(1,1) =  r11 + r22 + r33
    fa(1,2) =  r23 - r32
    fa(1,3) =  r31 - r13
    fa(1,4) =  r12 - r21

    fa(2,1) =  r23 - r32
    fa(2,2) =  r11 - r22 - r33
    fa(2,3) =  r12 + r21
    fa(2,4) =  r13 + r31

    fa(3,1) =  r31 - r13
    fa(3,2) =  r12 + r21
    fa(3,3) = -r11 + r22 - r33
    fa(3,4) =  r23 + r32

    fa(4,1) =  r12 - r21
    fa(4,2) =  r13 + r31
    fa(4,3) =  r23 + r32
    fa(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesa(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fa, 4, eigenvaluesa, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point a!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    xsb(:) = 0.0d0
    xrb(:) = 0.0d0

    do  i = cv_item%grps(5)+1,cv_item%grps(6)
        ai = cv_item%lindexes(i)

        ! source
        xsb(:) = xsb(:) + x(:,ai)

        ! reference
        xrb(:) = xrb(:) + cv_item%xyz_str_b2%cvs(:,i-cv_item%grps(5))
    end do

    xsb(:) = xsb(:) * ingrb
    xrb(:) = xrb(:) * ingrb

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = cv_item%grps(5)+1,cv_item%grps(6)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3))
    end do

    r11 = r11 * ingrb
    r12 = r12 * ingrb
    r13 = r13 * ingrb

    r21 = r21 * ingrb
    r22 = r22 * ingrb
    r23 = r23 * ingrb

    r31 = r31 * ingrb
    r32 = r32 * ingrb
    r33 = r33 * ingrb

    ! construct matrix for quaterion fitting
    fb(1,1) =  r11 + r22 + r33
    fb(1,2) =  r23 - r32
    fb(1,3) =  r31 - r13
    fb(1,4) =  r12 - r21

    fb(2,1) =  r23 - r32
    fb(2,2) =  r11 - r22 - r33
    fb(2,3) =  r12 + r21
    fb(2,4) =  r13 + r31

    fb(3,1) =  r31 - r13
    fb(3,2) =  r12 + r21
    fb(3,3) = -r11 + r22 - r33
    fb(3,4) =  r23 + r32

    fb(4,1) =  r12 - r21
    fb(4,2) =  r13 + r31
    fb(4,3) =  r23 + r32
    fb(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesb(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fb, 4, eigenvaluesb, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
    ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
    ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
    ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
    ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
    ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
    ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    ! rotation matrix b  ------------------------------
    ub(1,1) = fb(1,best)**2 + fb(2,best)**2 - fb(3,best)**2 - fb(4,best)**2
    ub(2,1) = 2.0d0*( fb(2,best)*fb(3,best) - fb(1,best)*fb(4,best) )
    ub(3,1) = 2.0d0*( fb(2,best)*fb(4,best) + fb(1,best)*fb(3,best) )

    ub(1,2) = 2.0d0*( fb(2,best)*fb(3,best) + fb(1,best)*fb(4,best) )
    ub(2,2) = fb(1,best)**2 - fb(2,best)**2 + fb(3,best)**2 - fb(4,best)**2
    ub(3,2) = 2.0d0*( fb(3,best)*fb(4,best) - fb(1,best)*fb(2,best) )

    ub(1,3) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
    ub(2,3) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
    ub(3,3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2

    ! z-axis ===================================================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(7))) - x(:,cv_item%lindexes(cv_item%grps(8)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_yaxis2 = 1.0d0 / (yaxisr(1)**2 + yaxisr(2)**2 + yaxisr(3)**2)
    o_yaxis  = sqrt(o_yaxis2)
    yaxis(:) = yaxisr(:) * o_yaxis

    ! x-axis ===================================================================
    ! is cross product of y and z axes
    xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
    xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
    xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)

    ! get origins of bases
    oa(:) = 0.0d0
    ! move reference point to origin
    oa(:) = oa(:) - xra(:)
    ! rotate
    tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
    tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
    tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
    ! move origin to new reference point (experimental structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

    ! position of bp origin
    morg(:) = 0.5d0*(oa(:) + ob(:))

    ! axis
    mst(:,1) = xaxis(:)
    mst(:,2) = yaxis(:)
    mst(:,3) = zaxis(:)

end subroutine calculate_nasstp_getbp2

!===============================================================================

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
    integer             :: i,ai,info,best,mi,mj
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: o_zaxis2,o_zaxis
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zaxisr(3),zsc,y0axis(3)
    real(PMFDP)         :: tmp1(3),ingra,ingrb,oa(3),ob(3)
    real(PMFDP)         :: yaxisr(3),o_yaxis,o_yaxis2

    real(PMFDP)         :: a_fa(4),a_fb(4),a_zaxis(3),a_rij(4,4),a_xaxis(3),a_y0axis(3)
    real(PMFDP)         :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4),a_xsa(3),a_xsb(3)
    real(PMFDP)         :: a_yaxisr(3),a_ua(3,3),a_ub(3,3),a_yaxis(3)
    real(PMFDP)         :: t1,t2
    ! --------------------------------------------------------------------------

! inverse number of atoms
    ingra = 1.0d0 / (cv_item%grps(5)-cv_item%grps(4))
    ingrb = 1.0d0 / (cv_item%grps(6)-cv_item%grps(5))

    ! calculate geometrical centres (source and target) -------------------
    xsa(:) = 0.0d0
    xra(:) = 0.0d0

    do  i = cv_item%grps(4)+1,cv_item%grps(5)
        ai = cv_item%lindexes(i)
        ! source
        xsa(:) = xsa(:) + x(:,ai)

        ! reference
        xra(:) = xra(:) + cv_item%xyz_str_a2%cvs(:,i-cv_item%grps(4))
    end do

    xsa(:) = xsa(:) * ingra
    xra(:) = xra(:) * ingra

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = cv_item%grps(4)+1,cv_item%grps(5)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3))
    end do

    r11 = r11 * ingra
    r12 = r12 * ingra
    r13 = r13 * ingra

    r21 = r21 * ingra
    r22 = r22 * ingra
    r23 = r23 * ingra

    r31 = r31 * ingra
    r32 = r32 * ingra
    r33 = r33 * ingra

    ! construct matrix for quaterion fitting
    fa(1,1) =  r11 + r22 + r33
    fa(1,2) =  r23 - r32
    fa(1,3) =  r31 - r13
    fa(1,4) =  r12 - r21

    fa(2,1) =  r23 - r32
    fa(2,2) =  r11 - r22 - r33
    fa(2,3) =  r12 + r21
    fa(2,4) =  r13 + r31

    fa(3,1) =  r31 - r13
    fa(3,2) =  r12 + r21
    fa(3,3) = -r11 + r22 - r33
    fa(3,4) =  r23 + r32

    fa(4,1) =  r12 - r21
    fa(4,2) =  r13 + r31
    fa(4,3) =  r23 + r32
    fa(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesa(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fa, 4, eigenvaluesa, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point a!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    xsb(:) = 0.0d0
    xrb(:) = 0.0d0

    do  i = cv_item%grps(5)+1,cv_item%grps(6)
        ai = cv_item%lindexes(i)

        ! source
        xsb(:) = xsb(:) + x(:,ai)

        ! reference
        xrb(:) = xrb(:) + cv_item%xyz_str_b2%cvs(:,i-cv_item%grps(5))
    end do

    xsb(:) = xsb(:) * ingrb
    xrb(:) = xrb(:) * ingrb

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = cv_item%grps(5)+1,cv_item%grps(6)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3))
    end do

    r11 = r11 * ingrb
    r12 = r12 * ingrb
    r13 = r13 * ingrb

    r21 = r21 * ingrb
    r22 = r22 * ingrb
    r23 = r23 * ingrb

    r31 = r31 * ingrb
    r32 = r32 * ingrb
    r33 = r33 * ingrb

    ! construct matrix for quaterion fitting
    fb(1,1) =  r11 + r22 + r33
    fb(1,2) =  r23 - r32
    fb(1,3) =  r31 - r13
    fb(1,4) =  r12 - r21

    fb(2,1) =  r23 - r32
    fb(2,2) =  r11 - r22 - r33
    fb(2,3) =  r12 + r21
    fb(2,4) =  r13 + r31

    fb(3,1) =  r31 - r13
    fb(3,2) =  r12 + r21
    fb(3,3) = -r11 + r22 - r33
    fb(3,4) =  r23 + r32

    fb(4,1) =  r12 - r21
    fb(4,2) =  r13 + r31
    fb(4,3) =  r23 + r32
    fb(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesb(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fb, 4, eigenvaluesb, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
    ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
    ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
    ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
    ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
    ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
    ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    ! rotation matrix b  ------------------------------
    ub(1,1) = fb(1,best)**2 + fb(2,best)**2 - fb(3,best)**2 - fb(4,best)**2
    ub(2,1) = 2.0d0*( fb(2,best)*fb(3,best) - fb(1,best)*fb(4,best) )
    ub(3,1) = 2.0d0*( fb(2,best)*fb(4,best) + fb(1,best)*fb(3,best) )

    ub(1,2) = 2.0d0*( fb(2,best)*fb(3,best) + fb(1,best)*fb(4,best) )
    ub(2,2) = fb(1,best)**2 - fb(2,best)**2 + fb(3,best)**2 - fb(4,best)**2
    ub(3,2) = 2.0d0*( fb(3,best)*fb(4,best) - fb(1,best)*fb(2,best) )

    ub(1,3) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
    ub(2,3) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
    ub(3,3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2

    ! z-axis ===================================================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(7))) - x(:,cv_item%lindexes(cv_item%grps(8)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_yaxis2 = 1.0d0 / (yaxisr(1)**2 + yaxisr(2)**2 + yaxisr(3)**2)
    o_yaxis  = sqrt(o_yaxis2)
    yaxis(:) = yaxisr(:) * o_yaxis

    ! x-axis ===================================================================
    ! is cross product of y and z axes
    xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
    xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
    xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)

    ! get origins of bases
    oa(:) = 0.0d0
    ! move reference point to origin
    oa(:) = oa(:) - xra(:)
    ! rotate
    tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
    tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
    tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
    ! move origin to new reference point (experimental structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

! ==============================================================================

! final derivatives
    a_ua(:,:)   = 0.0d0
    a_ub(:,:)   = 0.0d0
    a_y0axis(:) = 0.0d0
    a_xsa(:)    = 0.0d0
    a_xsb(:)    = 0.0d0

! ************
! ==== a_xaxis
! ************

    a_xaxis(:) = a_mst(:,1)

! a_xaxis with respect to ua and ub
! with respect to yaxis
!     xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
!     xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
!     xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)
    a_yaxis(1) = - a_xaxis(2)*zaxis(3) + zaxis(2)*a_xaxis(3)
    a_yaxis(2) = - a_xaxis(3)*zaxis(1) + zaxis(3)*a_xaxis(1)
    a_yaxis(3) = - a_xaxis(1)*zaxis(2) + zaxis(1)*a_xaxis(2)

! with respect to yaxisr
    t1 = yaxisr(1)*a_yaxis(1) + yaxisr(2)*a_yaxis(2) + yaxisr(3)*a_yaxis(3)
    a_yaxisr(1) =   a_yaxis(1)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(1)*t1
    a_yaxisr(2) =   a_yaxis(2)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(2)*t1
    a_yaxisr(3) =   a_yaxis(3)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(3)*t1

! with respect to y0axis
    t1 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_y0axis(1) = a_yaxisr(1) - zaxis(1)*t1
    a_y0axis(2) = a_yaxisr(2) - zaxis(2)*t1
    a_y0axis(3) = a_yaxisr(3) - zaxis(3)*t1

! with respect to zaxis
    t1 = y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3)
    t2 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_zaxis(1) = - t1*a_yaxisr(1) - y0axis(1)*t2
    a_zaxis(2) = - t1*a_yaxisr(2) - y0axis(2)*t2
    a_zaxis(3) = - t1*a_yaxisr(3) - y0axis(3)*t2

! with respect to zaxis
!     xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
!     xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
!     xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)
    a_zaxis(1) = a_zaxis(1) - yaxis(2)*a_xaxis(3) + yaxis(3)*a_xaxis(2)
    a_zaxis(2) = a_zaxis(2) - yaxis(3)*a_xaxis(1) + yaxis(1)*a_xaxis(3)
    a_zaxis(3) = a_zaxis(3) - yaxis(1)*a_xaxis(2) + yaxis(2)*a_xaxis(1)

    t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
    a_ua(1,3) = a_ua(1,3) + 0.5d0*a_zaxis(1)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ua(2,3) = a_ua(2,3) + 0.5d0*a_zaxis(2)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ua(3,3) = a_ua(3,3) + 0.5d0*a_zaxis(3)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(3)*t1

    a_ub(1,3) = a_ub(1,3) + 0.5d0*zsc*a_zaxis(1)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ub(2,3) = a_ub(2,3) + 0.5d0*zsc*a_zaxis(2)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ub(3,3) = a_ub(3,3) + 0.5d0*zsc*a_zaxis(3)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(3)*t1


! ************
! ==== a_yaxis
! ************

    a_yaxis(:) = a_mst(:,2)

! with respect to yaxisr
    t1 = yaxisr(1)*a_yaxis(1) + yaxisr(2)*a_yaxis(2) + yaxisr(3)*a_yaxis(3)
    a_yaxisr(1) =   a_yaxis(1)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(1)*t1
    a_yaxisr(2) =   a_yaxis(2)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(2)*t1
    a_yaxisr(3) =   a_yaxis(3)*o_yaxis - o_yaxis*o_yaxis2*yaxisr(3)*t1
! with respect to y0axis
    t1 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_y0axis(1) = a_y0axis(1) + a_yaxisr(1) - zaxis(1)*t1
    a_y0axis(2) = a_y0axis(2) + a_yaxisr(2) - zaxis(2)*t1
    a_y0axis(3) = a_y0axis(3) + a_yaxisr(3) - zaxis(3)*t1

! with respect to zaxis
    t1 = y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3)
    t2 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_zaxis(1) = - t1*a_yaxisr(1) - y0axis(1)*t2
    a_zaxis(2) = - t1*a_yaxisr(2) - y0axis(2)*t2
    a_zaxis(3) = - t1*a_yaxisr(3) - y0axis(3)*t2

    t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
    a_ua(1,3) = a_ua(1,3) + 0.5d0*a_zaxis(1)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ua(2,3) = a_ua(2,3) + 0.5d0*a_zaxis(2)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ua(3,3) = a_ua(3,3) + 0.5d0*a_zaxis(3)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(3)*t1

    a_ub(1,3) = a_ub(1,3) + 0.5d0*zsc*a_zaxis(1)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ub(2,3) = a_ub(2,3) + 0.5d0*zsc*a_zaxis(2)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ub(3,3) = a_ub(3,3) + 0.5d0*zsc*a_zaxis(3)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(3)*t1

! ************
! ==== a_zaxis
! ************

    a_zaxis(:) = a_mst(:,3)

! with respect to zaxis
    t1 = zaxisr(1)*a_zaxis(1) + zaxisr(2)*a_zaxis(2) + zaxisr(3)*a_zaxis(3)
    a_ua(1,3) = a_ua(1,3) + 0.5d0*a_zaxis(1)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ua(2,3) = a_ua(2,3) + 0.5d0*a_zaxis(2)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ua(3,3) = a_ua(3,3) + 0.5d0*a_zaxis(3)*o_zaxis - 0.5d0*o_zaxis*o_zaxis2*zaxisr(3)*t1

    a_ub(1,3) = a_ub(1,3) + 0.5d0*zsc*a_zaxis(1)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(1)*t1
    a_ub(2,3) = a_ub(2,3) + 0.5d0*zsc*a_zaxis(2)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(2)*t1
    a_ub(3,3) = a_ub(3,3) + 0.5d0*zsc*a_zaxis(3)*o_zaxis - 0.5d0*zsc*o_zaxis*o_zaxis2*zaxisr(3)*t1

! ************
! ==== origin
! ************

!     ! get origins of bases
!     oa(:) = 0.0d0
!     ! move reference point to origin
!     oa(:) = oa(:) - xra(:)
!     ! rotate
!     tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
!     tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
!     tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
!     ! move origin to new reference point (experiemntal structure)
!     oa(:) = tmp1(:) + xsa(:)
!
!     ob(:) = 0.0d0
!     ob(:) = ob(:) - xrb(:)
!     tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
!     tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
!     tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
!     ob(:) = tmp1(:) + xsb(:)

! a_morg
!    ! position of bp origin
!    morg(:) = 0.5d0*(oa(:) + ob(:))

! a_morg with respect to ua and ub
    a_ua(1,1) = a_ua(1,1) - 0.5d0*xra(1)*a_morg(1)
    a_ua(2,1) = a_ua(2,1) - 0.5d0*xra(1)*a_morg(2)
    a_ua(3,1) = a_ua(3,1) - 0.5d0*xra(1)*a_morg(3)
    a_ua(1,2) = a_ua(1,2) - 0.5d0*xra(2)*a_morg(1)
    a_ua(2,2) = a_ua(2,2) - 0.5d0*xra(2)*a_morg(2)
    a_ua(3,2) = a_ua(3,2) - 0.5d0*xra(2)*a_morg(3)
    a_ua(1,3) = a_ua(1,3) - 0.5d0*xra(3)*a_morg(1)
    a_ua(2,3) = a_ua(2,3) - 0.5d0*xra(3)*a_morg(2)
    a_ua(3,3) = a_ua(3,3) - 0.5d0*xra(3)*a_morg(3)

    a_ub(1,1) = a_ub(1,1) - 0.5d0*xrb(1)*a_morg(1)
    a_ub(2,1) = a_ub(2,1) - 0.5d0*xrb(1)*a_morg(2)
    a_ub(3,1) = a_ub(3,1) - 0.5d0*xrb(1)*a_morg(3)
    a_ub(1,2) = a_ub(1,2) - 0.5d0*xrb(2)*a_morg(1)
    a_ub(2,2) = a_ub(2,2) - 0.5d0*xrb(2)*a_morg(2)
    a_ub(3,2) = a_ub(3,2) - 0.5d0*xrb(2)*a_morg(3)
    a_ub(1,3) = a_ub(1,3) - 0.5d0*xrb(3)*a_morg(1)
    a_ub(2,3) = a_ub(2,3) - 0.5d0*xrb(3)*a_morg(2)
    a_ub(3,3) = a_ub(3,3) - 0.5d0*xrb(3)*a_morg(3)

    ! a_d with respect to xsa, xsb
    a_xsa(:) = + 0.5d0*ingra*a_morg(:)
    a_xsb(:) = + 0.5d0*ingrb*a_morg(:)

! rotation matrix a ------------------------------
!     ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
!     ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
!     ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    a_fa(1) = 2.0d0*( fa(1,best)*a_ua(1,1) - fa(4,best)*a_ua(2,1) + fa(3,best)*a_ua(3,1))
    a_fa(2) = 2.0d0*( fa(2,best)*a_ua(1,1) + fa(3,best)*a_ua(2,1) + fa(4,best)*a_ua(3,1))
    a_fa(3) = 2.0d0*(-fa(3,best)*a_ua(1,1) + fa(2,best)*a_ua(2,1) + fa(1,best)*a_ua(3,1))
    a_fa(4) = 2.0d0*(-fa(4,best)*a_ua(1,1) - fa(1,best)*a_ua(2,1) + fa(2,best)*a_ua(3,1))

    a_fb(1) = 2.0d0*( fb(1,best)*a_ub(1,1) - fb(4,best)*a_ub(2,1) + fb(3,best)*a_ub(3,1))
    a_fb(2) = 2.0d0*( fb(2,best)*a_ub(1,1) + fb(3,best)*a_ub(2,1) + fb(4,best)*a_ub(3,1))
    a_fb(3) = 2.0d0*(-fb(3,best)*a_ub(1,1) + fb(2,best)*a_ub(2,1) + fb(1,best)*a_ub(3,1))
    a_fb(4) = 2.0d0*(-fb(4,best)*a_ub(1,1) - fb(1,best)*a_ub(2,1) + fb(2,best)*a_ub(3,1))

!     ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
!     ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
!     ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    a_fa(1) = a_fa(1) + 2.0d0*(fa(4,best)*a_ua(1,2) + fa(1,best)*a_ua(2,2) - fa(2,best)*a_ua(3,2))
    a_fa(2) = a_fa(2) + 2.0d0*(fa(3,best)*a_ua(1,2) - fa(2,best)*a_ua(2,2) - fa(1,best)*a_ua(3,2))
    a_fa(3) = a_fa(3) + 2.0d0*(fa(2,best)*a_ua(1,2) + fa(3,best)*a_ua(2,2) + fa(4,best)*a_ua(3,2))
    a_fa(4) = a_fa(4) + 2.0d0*(fa(1,best)*a_ua(1,2) - fa(4,best)*a_ua(2,2) + fa(3,best)*a_ua(3,2))

    a_fb(1) = a_fb(1) + 2.0d0*(fb(4,best)*a_ub(1,2) + fb(1,best)*a_ub(2,2) - fb(2,best)*a_ub(3,2))
    a_fb(2) = a_fb(2) + 2.0d0*(fb(3,best)*a_ub(1,2) - fb(2,best)*a_ub(2,2) - fb(1,best)*a_ub(3,2))
    a_fb(3) = a_fb(3) + 2.0d0*(fb(2,best)*a_ub(1,2) + fb(3,best)*a_ub(2,2) + fb(4,best)*a_ub(3,2))
    a_fb(4) = a_fb(4) + 2.0d0*(fb(1,best)*a_ub(1,2) - fb(4,best)*a_ub(2,2) + fb(3,best)*a_ub(3,2))

!     ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
!     ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
!     ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    a_fa(1) = a_fa(1) + 2.0d0*(-fa(3,best)*a_ua(1,3) + fa(2,best)*a_ua(2,3) + fa(1,best)*a_ua(3,3))
    a_fa(2) = a_fa(2) + 2.0d0*( fa(4,best)*a_ua(1,3) + fa(1,best)*a_ua(2,3) - fa(2,best)*a_ua(3,3))
    a_fa(3) = a_fa(3) + 2.0d0*(-fa(1,best)*a_ua(1,3) + fa(4,best)*a_ua(2,3) - fa(3,best)*a_ua(3,3))
    a_fa(4) = a_fa(4) + 2.0d0*( fa(2,best)*a_ua(1,3) + fa(3,best)*a_ua(2,3) + fa(4,best)*a_ua(3,3))

    a_fb(1) = a_fb(1) + 2.0d0*(-fb(3,best)*a_ub(1,3) + fb(2,best)*a_ub(2,3) + fb(1,best)*a_ub(3,3))
    a_fb(2) = a_fb(2) + 2.0d0*( fb(4,best)*a_ub(1,3) + fb(1,best)*a_ub(2,3) - fb(2,best)*a_ub(3,3))
    a_fb(3) = a_fb(3) + 2.0d0*(-fb(1,best)*a_ub(1,3) + fb(4,best)*a_ub(2,3) - fb(3,best)*a_ub(3,3))
    a_fb(4) = a_fb(4) + 2.0d0*( fb(2,best)*a_ub(1,3) + fb(3,best)*a_ub(2,3) + fb(4,best)*a_ub(3,3))

! derivatives of fa with respect to matrix elements
    v(:,:) = fa(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvaluesa(i) - eigenvaluesa(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + fa(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

! merge xij with a_fa, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fa(1)*xij(1,mi,mj)+a_fa(2)*xij(2,mi,mj)+a_fa(3)*xij(3,mi,mj)+a_fa(4)*xij(4,mi,mj))*ingra
        end do
    end do

! finaly gradients for group_a
    do i = cv_item%grps(4) + 1, cv_item%grps(5)

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3)) &
                + a_xsa(1)

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3)) &
                + a_xsa(2)

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a2%cvs(1,i-cv_item%grps(4)) - xra(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a2%cvs(2,i-cv_item%grps(4)) - xra(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_a2%cvs(3,i-cv_item%grps(4)) - xra(3)) &
                + a_xsa(3)

    end do

! derivatives of fb with respect to matrix elements
    v(:,:) = fb(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvaluesb(i) - eigenvaluesb(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + fb(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

! merge xij with a_fb, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fb(1)*xij(1,mi,mj)+a_fb(2)*xij(2,mi,mj)+a_fb(3)*xij(3,mi,mj)+a_fb(4)*xij(4,mi,mj))*ingrb
        end do
    end do

! finaly gradients for group_b
    do i = cv_item%grps(5) + 1, cv_item%grps(6)

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3)) &
                + a_xsb(1)

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3)) &
                + a_xsb(2)

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_b2%cvs(1,i-cv_item%grps(5)) - xrb(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_b2%cvs(2,i-cv_item%grps(5)) - xrb(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_b2%cvs(3,i-cv_item%grps(5)) - xrb(3)) &
                + a_xsb(3)
    end do

! finally gradients for group_c, group_d
    ai = cv_item%lindexes(cv_item%grps(7))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y0axis(:)

    ai = cv_item%lindexes(cv_item%grps(8))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - a_y0axis(:)

end subroutine calculate_nasstp_getbp2_der

end module cv_nasstp

