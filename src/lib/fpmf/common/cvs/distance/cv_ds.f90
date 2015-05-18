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

module cv_ds

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeDS
    contains
        procedure :: load_cv        => load_ds
        procedure :: calculate_cv   => calculate_ds
end type CVTypeDS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_ds
!===============================================================================

subroutine load_ds(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use pmf_unit
    use cv_common

    implicit none
    class(CVTypeDS)                     :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'DS'
    cv_item%unit = pmf_unit_power_unit(LengthUnit,2)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_read_groups(cv_item,prm_fin)

end subroutine load_ds

!===============================================================================
! Subroutine:  calculate_ds
!===============================================================================

subroutine calculate_ds(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeDS)     :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: ai,m
    real(PMFDP)         :: d1(3),d2(3),dx(3)
    real(PMFDP)         :: totmass1,totmass2,amass
    ! -----------------------------------------------------------------------------

    ! calculate actual value
    totmass1 = 0.0d0
    d1(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d1(:) = d1(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_ds!')
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
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_ds!')
    end if
    d2(:) = d2(:) / totmass2

    dx(:) = d1(:) - d2(:)
    ctx%CVsValues(cv_item%idx) = dx(1)**2 + dx(2)**2 + dx(3)**2

    ! ------------------------------------------------

    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + 2.0d0*dx(:)*amass/totmass1
    end do

    do  m = cv_item%grps(1) + 1 , cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - 2.0d0*dx(:)*amass/totmass2
    end do

    return

end subroutine calculate_ds

!===============================================================================

end module cv_ds

