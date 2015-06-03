!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module gap_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  gap_core_main
!===============================================================================

subroutine gap_core_main

 use gap_dat
 use gap_output

 ! ------------------------------------------------------------------------------------------

 call gap_core_force
 call gap_output_write

end subroutine gap_core_main

!===============================================================================
! Subroutine:  gap_core_force
!===============================================================================

subroutine gap_core_force

 use pmf_dat
 use pmf_cvs
 use gp_dat_mod
 use gp_basic_mod
 use gap_dat

 implicit none
 integer                :: k, i, ci, n, j
 real(PMFDP)            :: gf
 ! ------------------------------------------------------------------------------------------

 TotalGAPEnergy = 0.0d0

 do k=1,NumOfGAPGroups
    do i=1,GAPGroupList(k)%nindexes
       ci = GAPGroupList(k)%gapcvindx(i)
       GAPGroupList(k)%values(i) = CVContext%CVsValues(GAPCVList(ci)%cvindx)
    end do

    if( gap_core_is_value_in_range(GAPGroupList(k)%values(:), k) ) then

        GAPGroupList(k)%energy = f_predict(GAPGroupList(k)%gp,GAPGroupList(k)%values,SE_kernel_r_rr)
        GAPGroupList(k)%forces = -f_predict_grad(GAPGroupList(k)%gp,GAPGroupList(k)%values,SE_kernel_r_rr)

        TotalGAPEnergy = TotalGAPEnergy + GAPGroupList(k)%energy
    else
        GAPGroupList(k)%energy = 0.0d0
        GAPGroupList(k)%forces = 0.0d0
    end if
 end do 

 do i=1,NumOfGAPCVs
    ci = GAPCVList(i)%cvindx
    gf = GAPGroupList(GAPCVList(i)%groupid)%forces(GAPCVList(i)%indexid)
    do n=1,GAPCVList(i)%cv%nindatoms
       j = GAPCVList(i)%cv%indlindexes(n)
       Frc(:,j) = Frc(:,j) + gf*CVContext%CVsDrvs(:,j,ci)
    end do
    CVContext%CVsFrc(ci) = CVContext%CVsFrc(ci) + gf
 end do 

end subroutine gap_core_force

!===============================================================================
! Function:  gap_core_is_value_in_range
!===============================================================================

logical function gap_core_is_value_in_range(values, groupid)

    use gap_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: values(:)
    integer                :: groupid
    ! -----------------------------------------------
    integer                :: i, ci
    real(PMFDP)            :: val
    ! -----------------------------------------------

    gap_core_is_value_in_range = .true.

    do i=1,GAPGroupList(groupid)%nindexes
       ci = GAPGroupList(groupid)%gapcvindx(i)
       if( .not. GAPCVList(ci)%cv%is_periodic_cv() ) then
           if( (values(i) < GAPCVList(ci)%min_value) .or. (values(i) > GAPCVList(ci)%max_value) ) then
               gap_core_is_value_in_range = .false.
               return
           end if
       else
           ! move value into default range
           val = values(i) - GAPCVList(ci)%cv%get_period_cv_value() / 2.0d0 - GAPCVList(ci)%cv%get_min_cv_value()
           val = val - nint(values(i)/GAPCVList(ci)%cv%get_period_cv_value())*GAPCVList(ci)%cv%get_period_cv_value()
           val = val + GAPCVList(ci)%cv%get_period_cv_value() / 2.0d0 + GAPCVList(ci)%cv%get_min_cv_value()
           ! min_value < max_value
           if( GAPCVList(ci)%min_value < GAPCVList(ci)%max_value ) then
               if( (val < GAPCVList(ci)%min_value) .or. (val > GAPCVList(ci)%max_value) ) then
                   gap_core_is_value_in_range = .false.
                   return
               end if
           else
           ! min_value > max_value
               if( (val < GAPCVList(ci)%min_value) .and. (val > GAPCVList(ci)%max_value) ) then
                   gap_core_is_value_in_range = .false.
                   return
               end if
           end if
       end if
    end do

    return

end function gap_core_is_value_in_range

!===============================================================================

end module gap_core
