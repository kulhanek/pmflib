!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module cv_dd

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeDD
    contains
        procedure :: load_cv        => load_dd
        procedure :: calculate_cv   => calculate_dd
end type CVTypeDD

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_dd
!===============================================================================

subroutine load_dd(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common

    implicit none
    class(CVTypeDD)                    :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'DD'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 4
    call cv_common_read_groups(cv_item,prm_fin)

end subroutine load_dd

!===============================================================================
! Subroutine:  calculate_dd
!===============================================================================

subroutine calculate_dd(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeDD)     :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: ai,m
    real(PMFDP)    :: d1(3),d2(3),d3(3),d4(3),dx(3),dy(3)
    real(PMFDP)    :: totmass1,totmass2,totmass3,totmass4,amass
    real(PMFDP)    :: sqvp1,sqvp2,sc1,sc2
    ! -----------------------------------------------------------------------------

    ! calculate first distance ----------------------
    totmass1 = 0.0d0
    d1(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d1(:) = d1(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_dd!')
    end if
    d1(:) = d1(:) / totmass1

    totmass2 = 0.0d0
    d2(:) = 0.0d0
    do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d2(:) = d2(:) + x(:,ai)*amass
        totmass2 = totmass2 + amass
    end do
    if( totmass2 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_dd!')
    end if
    d2(:) = d2(:) / totmass2

    dx(:) = d1(:) - d2(:)
    if( fenable_pbc ) then
        call pmf_pbc_image_vector(dx)
    end if
    sqvp1 = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)

    ! calculate second distance  --------------------
    totmass3 = 0.0d0
    d3(:) = 0.0
    do  m = cv_item%grps(2) + 1, cv_item%grps(3)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d3(:) = d3(:) + x(:,ai)*amass
        totmass3 = totmass3 + amass
    end do
    if( totmass3 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass3 is zero in calculate_dd!')
    end if
    d3(:) = d3(:) / totmass3

    totmass4 = 0.0d0
    d4(:) = 0.0d0
    do  m = cv_item%grps(3) + 1 , cv_item%grps(4)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d4(:) = d4(:) + x(:,ai)*amass
        totmass4 = totmass4 + amass
    end do
    if( totmass4 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass4 is zero in calculate_dd!')
    end if
    d4(:) = d4(:) / totmass4

    dy(:) = d3(:) - d4(:)
    if( fenable_pbc ) then
        call pmf_pbc_image_vector(dy)
    end if
    sqvp2 = sqrt(dy(1)**2 + dy(2)**2 + dy(3)**2)

    ctx%CVsValues(cv_item%idx)  = sqvp1 - sqvp2

    ! calculate forces ------------------------------
    if( sqvp1 .gt. 1e-7 ) then
        sc1 = 1.0d0 / sqvp1
    else
        sc1 = 0.0
    end if

    if( sqvp2 .gt. 1e-7 ) then
        sc2 = 1.0d0 / sqvp2
    else
        sc2 = 0.0
    end if

    ! warning - groups can overlap - it is therefore important to add gradients

    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + sc1*dx(:)*amass/totmass1
    end do

    do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - sc1*dx(:)*amass/totmass2
    end do

    do  m = cv_item%grps(2) + 1 , cv_item%grps(3)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - sc2*dy(:)*amass/totmass3
    end do

    do  m = cv_item%grps(3) + 1 , cv_item%grps(4)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + sc2*dy(:)*amass/totmass4
    end do

end subroutine calculate_dd

!===============================================================================

end module cv_dd

