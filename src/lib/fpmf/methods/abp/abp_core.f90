!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more detajls.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module abp_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_core_main
! this is leap-frog and velocity-verlet ABP version
!===============================================================================

subroutine abp_core_main

    use abp_trajectory
    use abp_restart
    use abp_output
    use abp_dat

    ! --------------------------------------------------------------------------

    call abp_core_force
    call abp_output_write
    call abp_trajectory_write_snapshot
    call abp_restart_update

end subroutine abp_core_main

!===============================================================================
! Subroutine:  abp_core_force
! this is leap-frog and velocity-verlet ABP version
!===============================================================================

subroutine abp_core_force

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abp_dat
    use abp_accumulator

    implicit none
    integer                :: i,j,ci
    ! --------------------------------------------------------------------------

    ! save coordinate value to history --------------
    do i=1,NumOfABPCVs
        ci = ABPCVList(i)%cvindx
        cvvalues(i) = CVContext%CVsValues(ci)
    end do

    ! calculate abp force to be applied -------------
    select case(feimode)
        case(1)
            call abp_accumulator_get_la_hramp
        case default
            call pmf_utils_exit(PMF_OUT,1, &
                        '[ABP] Not implemented extrapolation/interpolation mode!')
    end select

    ! project abp force along coordinate ------------
    do j=1,NumOfLAtoms
        do i=1,NumOfABPCVs
            ci = ABPCVList(i)%cvindx
            Frc(:,j) = Frc(:,j) + ftemp * PMF_Rgas * la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! rest of ABP stuff -----------------------------
    select case(fmode)
        case(1)
            call abp_accumulator_update_direct
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Not implemented fmode in abp_core_force!')
    end select

end subroutine abp_core_force

!===============================================================================

end module abp_core
