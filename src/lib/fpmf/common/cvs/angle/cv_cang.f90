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

! cos(angle) - three points (i-j-k)

module cv_cang

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeCANG
    contains
        procedure :: load_cv        => load_cang
        procedure :: calculate_cv   => calculate_cang
end type CVTypeCANG

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_cang
!===============================================================================

subroutine load_cang(cv_item,prm_fin)

    use prmfile

    implicit none
    class(CVTypeCANG)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'CANG'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 3
    call cv_common_read_groups(cv_item,prm_fin)

end subroutine load_cang

!===============================================================================
! Subroutine:  calculate_cang
!===============================================================================

subroutine calculate_cang(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeCANG)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: ai,m
    real(PMFDP)    :: arg
    real(PMFDP)    :: x1,y1,z1,x2,y2,z2,x3,y3,z3
    real(PMFDP)    :: rijx,rijy,rijz,rij,rij2
    real(PMFDP)    :: rkjx,rkjy,rkjz,rkj,rkj2
    real(PMFDP)    :: one_rij2,one_rij
    real(PMFDP)    :: one_rkj2,one_rkj
    real(PMFDP)    :: one_rij2rkj2,one_rijrkj
    real(PMFDP)    :: one_rij3rkj,one_rijrkj3
    real(PMFDP)    :: argone_rij2,argone_rkj2
    real(PMFDP)    :: a_xix,a_xiy,a_xiz
    real(PMFDP)    :: a_xjx,a_xjy,a_xjz
    real(PMFDP)    :: a_xkx,a_xky,a_xkz
    real(PMFDP)    :: tmp, amass, totmass1, totmass2, totmass3
    ! -----------------------------------------------------------------------------

    ! first point ----------
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
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_cang!')
    end if
    x1 = x1 / totmass1;
    y1 = y1 / totmass1;
    z1 = z1 / totmass1;

    ! second point ----------
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
    if( totmass2 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_cang!')
    end if
    x2 = x2 / totmass2;
    y2 = y2 / totmass2;
    z2 = z2 / totmass2;

    ! third point ----------
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
    if( totmass3 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass3 is zero in calculate_cang!')
    end if
    x3 = x3 / totmass3;
    y3 = y3 / totmass3;
    z3 = z3 / totmass3;

    ! calc. cos of angle ----------

    rijx = x1 - x2
    rijy = y1 - y2
    rijz = z1 - z2

    rkjx = x3 - x2
    rkjy = y3 - y2
    rkjz = z3 - z2

    if( fenable_pbc ) then
        call pmf_pbc_image_vector3(rijx,rijy,rijz)
        call pmf_pbc_image_vector3(rkjx,rkjy,rkjz)
    end if

    rij2 = rijx*rijx + rijy*rijy + rijz*rijz
    rkj2 = rkjx*rkjx + rkjy*rkjy + rkjz*rkjz

    rij = sqrt(rij2);
    rkj = sqrt(rkj2);

    one_rij2 = 1.0d0 / rij2
    one_rij  = 1.0d0 / rij
    one_rkj2 = 1.0d0 / rkj2
    one_rkj  = 1.0d0 / rkj

    one_rij2rkj2= one_rij2*one_rkj2
    one_rijrkj  = one_rij*one_rkj
    one_rij3rkj = one_rij2rkj2*one_rij*rkj
    one_rijrkj3 = one_rij2rkj2*one_rkj*rij

    arg = (rijx*rkjx + rijy*rkjy + rijz*rkjz)*one_rijrkj

    ! coordinate ctx%CVsValues(cv_item%idx) ------------------------------------------------------------

    if ( arg .gt.  1.0 ) then
        arg =  1.0
    else if ( arg .lt. -1.0 ) then
        arg = -1.0
    end if

    ctx%CVsValues(cv_item%idx) = arg

    ! ------------------------------------------------------------------------------

    argone_rij2 = arg*one_rij2;
    argone_rkj2 = arg*one_rkj2;

    ! first derivatives -----------------------------------------------------------

    a_xix = (rkjx*one_rijrkj - argone_rij2*rijx)
    a_xiy = (rkjy*one_rijrkj - argone_rij2*rijy)
    a_xiz = (rkjz*one_rijrkj - argone_rij2*rijz)

    a_xkx = (rijx*one_rijrkj - argone_rkj2*rkjx)
    a_xky = (rijy*one_rijrkj - argone_rkj2*rkjy)
    a_xkz = (rijz*one_rijrkj - argone_rkj2*rkjz)

    a_xjx = -(a_xix + a_xkx)
    a_xjy = -(a_xiy + a_xky)
    a_xjz = -(a_xiz + a_xkz)

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

return

end subroutine calculate_cang

!===============================================================================

end module cv_cang

