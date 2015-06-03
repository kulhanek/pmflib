!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf_trajectory

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_trajectory_open
!===============================================================================

subroutine abf_trajectory_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abf_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only if ftrjsample > 0

    call pmf_utils_open(ABF_TRJ,fabftrj,'R')

    write(ABF_TRJ,10)

    return

    10 format('# ABFTRAJ')

end subroutine abf_trajectory_open

!===============================================================================
! Subroutine:  abf_trajectory_write_snapshot
!===============================================================================

subroutine abf_trajectory_write_snapshot

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abf_dat
    use abf_accumulator

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    if( mod(fstep,ftrjsample) .ne. 0 ) return

    ! write time
    write(ABF_TRJ,10) fstep

    ! write accumulator
    call abf_accumulator_write(ABF_TRJ)

    return

10 format('# ABFSNAP',I7)

end subroutine abf_trajectory_write_snapshot

!===============================================================================
! Subroutine:  abf_trajectory_close
!===============================================================================

subroutine abf_trajectory_close

    use pmf_constants
    use pmf_dat
    use abf_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    close(ABF_TRJ)

    return

end subroutine abf_trajectory_close

!===============================================================================

end module abf_trajectory
