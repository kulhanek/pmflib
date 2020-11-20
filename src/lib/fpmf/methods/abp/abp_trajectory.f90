!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module abp_trajectory

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_trajectory_open
!===============================================================================

subroutine abp_trajectory_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abp_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    call pmf_utils_open(ABP_TRJ,fabptrj,'R')

    write(ABP_TRJ,10)

    return

    10 format('# ABPTRAJ')

end subroutine abp_trajectory_open

!===============================================================================
! Subroutine:  abp_trajectory_write_snapshot
!===============================================================================

subroutine abp_trajectory_write_snapshot

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abp_dat
    use abp_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    if( mod(fstep,ftrjsample) .ne. 0 ) return

    ! write time
    write(ABP_TRJ,10) fstep

    ! write accumulator
    call abp_accumulator_write(ABP_TRJ)

    return

    10 format('# ABPSNAP ',I7)

end subroutine abp_trajectory_write_snapshot

!===============================================================================
! Subroutine:  abp_trajectory_close
!===============================================================================

subroutine abp_trajectory_close

    use pmf_constants
    use pmf_dat
    use abp_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    close(ABP_TRJ)

    return

end subroutine abp_trajectory_close

!===============================================================================

end module abp_trajectory
