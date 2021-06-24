!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtd_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_core_main
!===============================================================================

subroutine mtd_core_main

    use mtd_output
    use mtd_client
    use mtd_restart
    use mtd_trajectory
    ! --------------------------------------------------------------------------

    call mtd_core_force
    call mtd_trajectory_write_snapshot
    call mtd_restart_update
    call mtd_client_exchange_data(.false.)

end subroutine mtd_core_main

!===============================================================================
! Subroutine:  mtd_core_force
!===============================================================================

subroutine mtd_core_force

    use pmf_utils
    use mtd_dat

    implicit none
    ! -----------------------------------------------------------------------------

    select case(fmode)
        case(0)
            return
        case(1)
            call mtd_core_grid()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmode!')
    end select

end subroutine mtd_core_force

!===============================================================================
! Subroutine:  mtd_core_grid
! grid metadynamics algorithm
!===============================================================================

subroutine mtd_core_grid

    use pmf_dat
    use mtd_dat
    use mtd_accu

    implicit none
    integer                    :: i,j,ci
    ! --------------------------------------------------------------------------

    do i = 1,NumOfMTDCVs
        ci = MTDCVList(i)%cvindx
        CVValues(i) = CVContext%CVsValues(ci)
    end do

    ! get current data
    call mtd_accu_get_data(CVValues, TotalMTDEnergy, MTDForce)

    ! put new hill if it is a time
    if( fstep .eq. meta_next_fstep ) then
        call mtd_core_add_new_hill(CVValues, TotalMTDEnergy)
    end if

    ! now update forces =======================================
    do i = 1,NumOfMTDCVs
        ci = MTDCVList(i)%cvindx
        do j = 1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + CVContext%CVsDrvs(:,j,ci) * MTDForce(i)
        end do
    end do

end subroutine mtd_core_grid

!===============================================================================
! Subroutine:  mtd_core_add_new_hill
!===============================================================================

subroutine mtd_core_add_new_hill(cvs,mtdpot)

    use pmf_dat
    use mtd_dat
    use mtd_accu
    use pmf_utils
    use mtd_output

    implicit none
    real(PMFDP)     :: cvs(:)
    real(PMFDP)     :: mtdpot
    ! --------------------------------------------
    real(PMFDP)     :: hill_height
    ! --------------------------------------------------------------------------

    ! default
    hill_height = fheight
    meta_next_fstep = fstep + fmetastep
    meta_step = meta_step + 1

    ! well-tempered MTD
    if( fmetatemp .gt. 0 ) then
        ! what to change?
        select case(fmetavary)
            case(0)
                ! no metavary
            case(1)
                ! vary the height
                hill_height = fheight * exp( - mtdpot / (PMF_Rgas*fmetatemp) )
            case(2)
                ! vary the step
                meta_next_fstep = fstep + nint( fmetastep * exp( mtdpot / (PMF_Rgas*fmetatemp) ) )
            case default
                call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmetavary in mtd_core_add_new_hill!')
        end select
    end if

    ! put the hill
    call mtd_accu_add_data(cvs,hill_height,CVWidths)

    ! report if required
    call mtd_output_write(cvs,hill_height,CVWidths)

end subroutine mtd_core_add_new_hill

!===============================================================================

end module mtd_core
