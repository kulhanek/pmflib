!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz 
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_dih2

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeDIH2
    contains
        procedure :: load_cv        => load_dih2
        procedure :: calculate_cv   => calculate_dih2
        ! PBC related methods
        procedure   :: is_periodic_cv       => is_periodic_cv_dih2
        procedure   :: get_min_cv_value     => get_min_cv_value_dih2
        procedure   :: get_max_cv_value     => get_max_cv_value_dih2
end type CVTypeDIH2

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_dih2
!===============================================================================

subroutine load_dih2(cv_item,prm_fin)

    use prmfile

    implicit none
    class(CVTypeDIH2)                    :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'DIH2'
    cv_item%unit          = AngleUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 6
    call cv_common_read_groups(cv_item,prm_fin)

end subroutine load_dih2

!===============================================================================
! Function:  is_periodic_cv_dih2
!===============================================================================

logical function is_periodic_cv_dih2(cv_item)

    implicit none
    class(CVTypeDIH2) :: cv_item
    ! --------------------------------------------------------------------------

    is_periodic_cv_dih2 = .true.

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function is_periodic_cv_dih2

!===============================================================================
! Function:  get_min_cv_value_dih2
!===============================================================================

real(PMFDP) function get_min_cv_value_dih2(cv_item)

    implicit none
    class(CVTypeDIH2) :: cv_item
    ! --------------------------------------------------------------------------

    get_min_cv_value_dih2 = -PMF_PI

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_min_cv_value_dih2

!===============================================================================
! Function:  get_max_cv_value_dih2
!===============================================================================

real(PMFDP) function get_max_cv_value_dih2(cv_item)

    implicit none
    class(CVTypeDIH2) :: cv_item
    ! --------------------------------------------------------------------------

    get_max_cv_value_dih2 = PMF_PI

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_max_cv_value_dih2

!===============================================================================
! Subroutine:  calculate_dih2
!===============================================================================

subroutine calculate_dih2(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeDIH2)    :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: m,ai
    real(PMFDP)    :: totmass1,totmass2,totmass3,totmass4,totmass5,totmass6,amass,tmp
    real(PMFDP)    :: x1,y1,z1,x2,y2,z2
    real(PMFDP)    :: x3,y3,z3,x4,y4,z4
    real(PMFDP)    :: x5,y5,z5,x6,y6,z6
    real(PMFDP)    :: rijx,rijy,rijz
    real(PMFDP)    :: rkjx,rkjy,rkjz,rkj2,rkj,one_rkj2,one_rkj4
    real(PMFDP)    :: rklx,rkly,rklz
    real(PMFDP)    :: dx,dy,dz,d2,one_d2
    real(PMFDP)    :: gx,gy,gz,g2,one_g2
    real(PMFDP)    :: rkjijx,rkjijy,rkjijz
    real(PMFDP)    :: rkjklx,rkjkly,rkjklz
    real(PMFDP)    :: a_xix,a_xiy,a_xiz
    real(PMFDP)    :: a_xjx,a_xjy,a_xjz
    real(PMFDP)    :: a_xkx,a_xky,a_xkz
    real(PMFDP)    :: a_xlx,a_xly,a_xlz
    real(PMFDP)    :: rkj_d2
    real(PMFDP)    :: mrkj_g2
    real(PMFDP)    :: rijorkj, rkjorkl
    real(PMFDP)    :: WjA,WjB,WkA,WkB
    real(PMFDP)    :: s1, targ
    ! -----------------------------------------------------------------------------

    ! first point --------
    totmass1 = 0
    x1 = 0
    y1 = 0
    z1 = 0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        x1 = x1 + x(1,ai)*amass
        y1 = y1 + x(2,ai)*amass
        z1 = z1 + x(3,ai)*amass
        totmass1 = totmass1 + amass
    end do
    x1 = x1 / totmass1
    y1 = y1 / totmass1
    z1 = z1 / totmass1

    ! second point --------
    totmass2 = 0
    x2 = 0
    y2 = 0
    z2 = 0
    do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        x2 = x2 + x(1,ai)*amass
        y2 = y2 + x(2,ai)*amass
        z2 = z2 + x(3,ai)*amass
        totmass2 = totmass2 + amass
    end do
    x2 = x2 / totmass2
    y2 = y2 / totmass2
    z2 = z2 / totmass2

    ! third point --------
    totmass3 = 0
    x3 = 0
    y3 = 0
    z3 = 0
    do  m = cv_item%grps(2) + 1 , cv_item%grps(3)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        x3 = x3 + x(1,ai)*amass
        y3 = y3 + x(2,ai)*amass
        z3 = z3 + x(3,ai)*amass
        totmass3 = totmass3 + amass
    end do
    x3 = x3 / totmass3
    y3 = y3 / totmass3
    z3 = z3 / totmass3

    ! fourth point --------
    totmass4 = 0
    x4 = 0
    y4 = 0
    z4 = 0
    do  m = cv_item%grps(3) + 1 , cv_item%grps(4)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        x4 = x4 + x(1,ai)*amass
        y4 = y4 + x(2,ai)*amass
        z4 = z4 + x(3,ai)*amass
        totmass4 = totmass4 + amass
    end do
    x4 = x4 / totmass4
    y4 = y4 / totmass4
    z4 = z4 / totmass4

    ! fith point --------
    totmass5 = 0
    x5 = 0
    y5 = 0
    z5 = 0
    do  m = cv_item%grps(4) + 1 , cv_item%grps(5)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        x5 = x5 + x(1,ai)*amass
        y5 = y5 + x(2,ai)*amass
        z5 = z5 + x(3,ai)*amass
        totmass5 = totmass5 + amass
    end do
    x5 = x5 / totmass5
    y5 = y5 / totmass5
    z5 = z5 / totmass5

    ! sixth point --------
    totmass6 = 0
    x6 = 0
    y6 = 0
    z6 = 0
    do  m = cv_item%grps(5) + 1 , cv_item%grps(6)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        x4 = x6 + x(1,ai)*amass
        y4 = y6 + x(2,ai)*amass
        z4 = z6 + x(3,ai)*amass
        totmass6 = totmass6 + amass
    end do
    x6 = x6 / totmass6
    y6 = y6 / totmass6
    z6 = z6 / totmass6

    ! calc dih2edral angle

    rijx = x1 - x2
    rijy = y1 - y2
    rijz = z1 - z2

    rkjx = x3 - x4
    rkjy = y3 - y4
    rkjz = z3 - z4

    rklx = x5 - x6
    rkly = y5 - y6
    rklz = z5 - z6

    if( fenable_pbc ) then
        call pmf_pbc_image_vector3(rijx,rijy,rijz)
        call pmf_pbc_image_vector3(rkjx,rkjy,rkjz)
        call pmf_pbc_image_vector3(rklx,rkly,rklz)
    end if

    rkjijx = rkjx - rijx
    rkjijy = rkjy - rijy
    rkjijz = rkjz - rijz

    rkjklx = rkjx - rklx
    rkjkly = rkjy - rkly
    rkjklz = rkjz - rklz

    ! coordinate definition
    !
    ! d = rij x rkj
    ! g = rkj x rkl
    !
    ! s1 = (rkjy*rklz - rkjz*rkly)*rijx + (rkjz*rklx - rkjx*rklz)*rijy + (rkjx*rkly - rkjy*rklx)*rijz
    !
    ! ksi = sign(s1)arccos(d.g/(|d|.|g|))
    !

    ! d = rij x rkj
    dx = rijy*rkjz - rijz*rkjy
    dy = rijz*rkjx - rijx*rkjz
    dz = rijx*rkjy - rijy*rkjx

    d2 = dx**2 + dy**2 + dz**2

    one_d2 = 1.0 / d2

    ! g = rkj x rkl
    gx = rkjy*rklz - rkjz*rkly
    gy = rkjz*rklx - rkjx*rklz
    gz = rkjx*rkly - rkjy*rklx

    g2 = gx**2 + gy**2 + gz**2

    one_g2 = 1.0 / g2

    ! ctx%CVsValues(cv_item%idx) of coordinate

    s1 = (rkjy*rklz - rkjz*rkly)*rijx + (rkjz*rklx - rkjx*rklz)*rijy + (rkjx*rkly - rkjy*rklx)*rijz

    targ = (dx*gx + dy*gy + dz*gz)/sqrt(d2*g2)
    if( targ .gt. 1.0d0 ) then
        targ = 1.0d0
    end if
    if( targ .lt. -1.0d0 ) then
        targ = -1.0d0
    end if

    ctx%CVsValues(cv_item%idx) = sign(1.0d0,s1)*acos( targ )

    ! and it's first derivatives --------------------------

    rkj2 = rkjx**2 + rkjy**2 + rkjz**2

    rkj = sqrt(rkj2)

    one_rkj2 = 1.0 / rkj2
    one_rkj4 = one_rkj2 * one_rkj2

    rkj_d2 = rkj / d2
    mrkj_g2 = -rkj / g2

    a_xix = rkj_d2 * dx
    a_xiy = rkj_d2 * dy
    a_xiz = rkj_d2 * dz

    a_xlx = mrkj_g2 * gx
    a_xly = mrkj_g2 * gy
    a_xlz = mrkj_g2 * gz

    rijorkj = rijx * rkjx + rijy * rkjy + rijz * rkjz
    rkjorkl = rkjx * rklx + rkjy * rkly + rkjz * rklz

    WjA = rijorkj / rkj2 - 1
    WjB = rkjorkl / rkj2

    WkA = rkjorkl / rkj2 - 1
    WkB = rijorkj / rkj2

    a_xjx = WjA * a_xix - WjB * a_xlx
    a_xjy = WjA * a_xiy - WjB * a_xly
    a_xjz = WjA * a_xiz - WjB * a_xlz

    a_xkx = WkA * a_xlx - WkB * a_xix
    a_xky = WkA * a_xly - WkB * a_xiy
    a_xkz = WkA * a_xlz - WkB * a_xiz

    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        tmp = mass(ai) / totmass1
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + a_xix * tmp
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + a_xiy * tmp
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + a_xiz * tmp
    end do

    do  m = cv_item%grps(1) + 1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        tmp = mass(ai) / totmass2
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + a_xjx * tmp
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + a_xjy * tmp
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + a_xjz * tmp
    end do

    do  m = cv_item%grps(2) + 1, cv_item%grps(3)
        ai = cv_item%lindexes(m)
        tmp = mass(ai) / totmass3
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + a_xkx * tmp
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + a_xky * tmp
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + a_xkz * tmp
    end do

    do  m = cv_item%grps(3) + 1, cv_item%grps(4)
        ai = cv_item%lindexes(m)
        tmp = mass(ai) / totmass4
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + a_xlx * tmp
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + a_xly * tmp
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + a_xlz * tmp
    end do

    do  m = cv_item%grps(4) + 1, cv_item%grps(5)
        ai = cv_item%lindexes(m)
        tmp = mass(ai) / totmass5
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + a_xkx * tmp
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + a_xky * tmp
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + a_xkz * tmp
    end do

    do  m = cv_item%grps(5) + 1, cv_item%grps(6)
        ai = cv_item%lindexes(m)
        tmp = mass(ai) / totmass6
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + a_xlx * tmp
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + a_xly * tmp
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + a_xlz * tmp
    end do

    return

end subroutine calculate_dih2

!===============================================================================

end module cv_dih2

