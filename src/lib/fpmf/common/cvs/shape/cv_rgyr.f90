!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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

module cv_rgyr

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeRGYR
    contains
        procedure :: load_cv        => load_rgyr
        procedure :: calculate_cv   => calculate_rgyr
end type CVTypeRGYR

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_rgyr
!===============================================================================

subroutine load_rgyr(cv_item,prm_fin)

    use prmfile

    implicit none
    class(CVTypeRGYR)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'RGYR'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_read_groups(cv_item,prm_fin)

end subroutine load_rgyr

!===============================================================================
! Subroutine:  calculate_rgyr
!===============================================================================

subroutine calculate_rgyr(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeRGYR)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,ai
    real(PMFDP)    :: x1,y1,z1,x2,y2,z2,rsum,dv,dvv,amass,totmass,itotmass
    ! -----------------------------------------------------------------------------

    ! calculate centre of mass -------------------
    x1 = 0.0d0
    y1 = 0.0d0
    z1 = 0.0d0
    totmass = 0.0d0

    do  i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        x1 = x1 + x(1,ai)*amass
        y1 = y1 + x(2,ai)*amass
        z1 = z1 + x(3,ai)*amass
        totmass = totmass + amass
    end do
    itotmass = 1.0d0 / totmass
    x1 = x1 * itotmass
    y1 = y1 * itotmass
    z1 = z1 * itotmass

    ! calculate rgyr --------------------------------
    rsum = 0.0d0
    do  i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        x2 = x(1,ai) - x1
        y2 = x(2,ai) - y1
        z2 = x(3,ai) - z1
        rsum = rsum + (x2**2 + y2**2 + z2**2)*amass
    end do

    ctx%CVsValues(cv_item%idx) = sqrt(rsum  * itotmass)

    if( abs(ctx%CVsValues(cv_item%idx)) .lt. 1.0e-7 ) then
        call pmf_utils_exit(PMF_OUT,1,'Rgyr is smaller than 1.0e-7 in calculate_rgyr!')
    end if

    ! calculate derivatives -------------------------
    dv = itotmass / ctx%CVsValues(cv_item%idx)

    do  i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        dvv = dv*amass
        x2 = x(1,ai) - x1
        y2 = x(2,ai) - y1
        z2 = x(3,ai) - z1
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + dvv*x2
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + dvv*y2
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + dvv*z2
    end do

    return

end subroutine calculate_rgyr

!===============================================================================

end module cv_rgyr

