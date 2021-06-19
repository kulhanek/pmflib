!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2018-2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_naspsxy

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type
use cv_math

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNASPSXY

    type(XYZFILE_TYPE)  :: xyz_str_a
    type(SImpStrData)   :: simpdat_a
    type(XYZFILE_TYPE)  :: xyz_str_b
    type(SImpStrData)   :: simpdat_b
    real(PMFDP)         :: opfac

    contains
        procedure :: load_cv        => load_naspsxy
        procedure :: calculate_cv   => calculate_naspsxy
end type CVTypeNASPSXY

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_naspsxy
!===============================================================================

subroutine load_naspsxy(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeNASPSXY)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'NASPSXY'
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 4
    call cv_common_init_groups(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

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

    call init_xyz(cv_item%xyz_str_a)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a)

    if( cv_item%xyz_str_a%natoms .ne. cv_item%grps(1) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group A and reference structure A differs!')
    end if

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))
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

    call init_xyz(cv_item%xyz_str_b)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_b,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_b)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_b)

    if( cv_item%xyz_str_b%natoms .ne. (cv_item%grps(2) - cv_item%grps(1)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group B and reference structure B differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(1)+1, cv_item%grps(2)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
                write(tmpstr,110) i-cv_item%grps(1), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b%symbols(i-cv_item%grps(1)))
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

    ! FIXME - tunable?
    write(PMF_OUT,80)
    cv_item%opfac = 2.0d0
    write(PMF_OUT,85) cv_item%opfac

    return

 50 format('   == Base #1 ====================================')
 60 format('   == Base #2 ====================================')
 70 format('   ** reference structure: ',A)
 75 format('   == y-axis  ====================================')
 80 format('   -----------------------------------------------')
 85 format('   ** Scaling of Opening:  ',F10.2)
100 format('Atom mismatch between group A and reference A atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B and reference B atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

end subroutine load_naspsxy

!===============================================================================
! Subroutine:  calculate_naspsxy
!===============================================================================

subroutine calculate_naspsxy(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASPSXY) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3)
    real(PMFDP)         :: yaxisr(3),zaxisr(3),y0axis(3)
    real(PMFDP)         :: a_xaxis(3),a_yaxis(3),a_zaxis(3)
    real(PMFDP)         :: a_yaxisr(3),a_zaxisr(3),a_y0axis(3)
    real(PMFDP)         :: d(3),t1,zsc
    real(PMFDP)         :: sx,sy,op,cop,sop,a_op
    integer             :: ai
    ! --------------------------------------------------------------------------

    ! ALL a_xxx must be initialized to zero due to additive function of _der methods from cv_math
    a_ua(:,:)   = 0.0d0
    a_ub(:,:)   = 0.0d0
    a_oa(:)     = 0.0d0
    a_ob(:)     = 0.0d0
    a_xaxis(:)  = 0.0d0
    a_yaxis(:)  = 0.0d0
    a_zaxis(:)  = 0.0d0
    a_yaxisr(:) = 0.0d0
    a_zaxisr(:) = 0.0d0
    a_y0axis(:) = 0.0d0
    a_op        = 0.0d0

! superimpose bases ==============================
    call superimpose_str(cv_item,0,              cv_item%grps(1),x,cv_item%xyz_str_a,cv_item%simpdat_a,ua,oa)
    call superimpose_str(cv_item,cv_item%grps(1),cv_item%grps(2),x,cv_item%xyz_str_b,cv_item%simpdat_b,ub,ob)

! z-axis =========================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    call norm_vec(zaxisr,zaxis)

! y-axis =========================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    call norm_vec(yaxisr,yaxis)

! x-axis =========================================
    call get_cross_product(yaxis,zaxis,xaxis)

! get value and derivatives ======================
    d(:) = oa(:) - ob(:)
    ! shear
    sx = d(1)*xaxis(1) + d(2)*xaxis(2) + d(3)*xaxis(3)
    ! stretch
    sy = d(1)*yaxis(1) + d(2)*yaxis(2) + d(3)*yaxis(3)
    ! opening
    call get_vtors(ua(:,2),ub(:,2),zaxis,op)
    op = cv_item%opfac * op

    ! final value
    cop = cos(op)
    sop = sin(op)
    ctx%CVsValues(cv_item%idx) = sx * cop + sy * sop

    ! derivative with respect to sx
    a_oa(:) = a_oa(:) + xaxis(:)*cop
    a_ob(:) = a_ob(:) - xaxis(:)*cop
    a_xaxis(:) = a_xaxis(:) + d(:)*cop

    ! derivative with respect to sy
    a_oa(:) = a_oa(:) + yaxis(:)*sop
    a_ob(:) = a_ob(:) - yaxis(:)*sop
    a_yaxis(:) = a_yaxis(:) + d(:)*sop

    ! derivative with respect to op
    a_op = (- sx * sop + sy * cop) * cv_item%opfac

    call get_vtors_der(ua(:,2),ub(:,2),zaxis,a_op,a_ua(:,2),a_ub(:,2),a_zaxis)

! derivatives with respect to axes ===============
!! z-axis
!    ! mutual orientation of two z-axis
!    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
!    ! get z-axis as average of two axes
!    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
!    call norm_vec(zaxisr,zaxis)
!
!! y-axis
!    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
!    ! remove projections to z-axis
!    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
!    ! normalize
!    call norm_vec(yaxisr,yaxis)
!
!! x-axis
!    call get_cross_product(yaxis,zaxis,xaxis)

! x-axis
    call get_cross_product_der(yaxis,zaxis,a_xaxis,a_yaxis,a_zaxis)

! y-axis
    call norm_vec_der(yaxisr,a_yaxis,a_yaxisr)

    ! with respect to y0axis
    t1 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_y0axis(1) = a_yaxisr(1) - zaxis(1)*t1
    a_y0axis(2) = a_yaxisr(2) - zaxis(2)*t1
    a_y0axis(3) = a_yaxisr(3) - zaxis(3)*t1

    ! with respect to zaxis
    a_zaxis(1) = a_zaxis(1) - a_yaxisr(1) * ( 2.0d0*y0axis(1)*zaxis(1) + y0axis(2)*zaxis(2) + y0axis(3)*zaxis(3) ) &
                            - a_yaxisr(2) * y0axis(1)*zaxis(2) &
                            - a_yaxisr(3) * y0axis(1)*zaxis(3)

    a_zaxis(2) = a_zaxis(2) - a_yaxisr(1) * y0axis(2)*zaxis(1) &
                            - a_yaxisr(2) * ( y0axis(1)*zaxis(1) + 2.0d0*y0axis(2)*zaxis(2) + y0axis(3)*zaxis(3) ) &
                            - a_yaxisr(3) * y0axis(2)*zaxis(3)

    a_zaxis(3) = a_zaxis(3) - a_yaxisr(1) * y0axis(3)*zaxis(1) &
                            - a_yaxisr(2) * y0axis(3)*zaxis(2) &
                            - a_yaxisr(3) * ( y0axis(1)*zaxis(1) + y0axis(2)*zaxis(2) + 2.0d0*y0axis(3)*zaxis(3) )

! z-axis
    call norm_vec_der(zaxisr,a_zaxis,a_zaxisr)

    a_ua(:,3) = a_ua(:,3) + 0.5d0*a_zaxisr(:)
    a_ub(:,3) = a_ub(:,3) + 0.5d0*zsc*a_zaxisr(:)

! derivatives for superimposed bases =============
    call superimpose_str_der(cv_item,0,              cv_item%grps(1),ctx,cv_item%xyz_str_a,cv_item%simpdat_a,a_ua,a_oa)

    call superimpose_str_der(cv_item,cv_item%grps(1),cv_item%grps(2),ctx,cv_item%xyz_str_b,cv_item%simpdat_b,a_ub,a_ob)

! finaly gradients for group_c, group_d ==========
    ai = cv_item%lindexes(cv_item%grps(3))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_y0axis(:)
    ai = cv_item%lindexes(cv_item%grps(4))
    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - a_y0axis(:)

end subroutine calculate_naspsxy

!===============================================================================

end module cv_naspsxy

