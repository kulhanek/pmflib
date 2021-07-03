!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module rst_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  rst_core_main
!===============================================================================

subroutine rst_core_main

    use rst_dat
    use rst_output
    use rst_restart

    implicit none
    ! --------------------------------------------------------------------------

    call rst_core_force
    call rst_output_write_data
    call rst_restart_update

end subroutine rst_core_main

!===============================================================================
! Subroutine:  rst_core_force
! plain restraints
!===============================================================================

subroutine rst_core_force

    use rst_dat
    use pmf_dat
    use rst_restraints
    use pmf_utils
    use pmf_cvs
    use rst_accu

    implicit none
    integer                :: i, ci
    real(PMFDP)            :: rvalue
    ! --------------------------------------------------------------------------

    ! process increments ----------------------------
    do i=1,NumOfRSTCVs
        call rst_restraints_increment(RSTCVList(i))
    end do

    ! calculate energy and gradients ----------------
    TotalRstEnergy = 0.0
    do i=1,NumOfRSTCVs
        ci = RSTCVList(i)%cvindx
        rvalue = CVContext%CVsValues(ci)
        if( RSTCVList(i)%mode .ne. 'W' ) then
            RSTCVList(i)%deviation = RSTCVList(i)%cv%get_deviation(rvalue,RSTCVList(i)%target_value)
        else
            RSTCVList(i)%deviation = 0.0d0
            if( rvalue .lt. RSTCVList(i)%left_value ) then
                RSTCVList(i)%deviation = RSTCVList(i)%cv%get_deviation(rvalue,RSTCVList(i)%left_value)
            end if
            if( rvalue .gt. RSTCVList(i)%right_value ) then
                RSTCVList(i)%deviation = RSTCVList(i)%cv%get_deviation(rvalue,RSTCVList(i)%right_value)
            end if
        end if
        RSTCVList(i)%energy = 0.5d0*RSTCVList(i)%force_constant*RSTCVList(i)%deviation**2
        TotalRstEnergy = TotalRstEnergy + RSTCVList(i)%energy

        ! correct forces -----------------------------
        Frc(:,:) = Frc(:,:) - RSTCVList(i)%force_constant*RSTCVList(i)%deviation*CVContext%CVsDrvs(:,:,RSTCVList(i)%cvindx)
    end do

    call rst_accu_add_sample(CVContext%CVsValues)

end subroutine rst_core_force

!===============================================================================

end module rst_core

