!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf2_trajectory

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf2_trajectory_open
!===============================================================================

subroutine abf2_trajectory_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abf2_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only if ftrjsample > 0

    call pmf_utils_open(ABF_TRJ,fabftrj,'R')

    write(ABF_TRJ,10)

    return

    10 format('# ABFTRAJ')

end subroutine abf2_trajectory_open

!===============================================================================
! Subroutine:  abf2_trajectory_write_snapshot
!===============================================================================

subroutine abf2_trajectory_write_snapshot

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abf2_dat
    use abf2_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    if( mod(fstep,ftrjsample) .ne. 0 ) return

    ! write time
    write(ABF_TRJ,10) fstep

    ! write accumulator
    call abf2_accumulator_write(ABF_TRJ)

    return

10 format('# ABFSNAP',I7)

end subroutine abf2_trajectory_write_snapshot

!===============================================================================
! Subroutine:  abf2_trajectory_close
!===============================================================================

subroutine abf2_trajectory_close

    use pmf_constants
    use pmf_dat
    use abf2_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    close(ABF_TRJ)

    return

end subroutine abf2_trajectory_close

!===============================================================================

end module abf2_trajectory
