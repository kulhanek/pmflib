!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module cst_trajectory

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_trajectory_open
!===============================================================================

subroutine cst_trajectory_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only if ftrjsample > 0

    call pmf_utils_open(CST_TRJ,fcsttrj,'R')
    write(CST_TRJ,10)

    return

10 format('#ACCUTRAJ')

end subroutine cst_trajectory_open

!===============================================================================
! Subroutine:  cst_trajectory_write_snapshot
!===============================================================================

subroutine cst_trajectory_write_snapshot

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use cst_dat
    use cst_accu

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    if( mod(fstep,ftrjsample) .ne. 0 ) return

    ! write time
    write(CST_TRJ,10) fstep

    ! write accumulator
    call cst_accu_write(CST_TRJ)

    return

10 format('#ACCUSNAP',1X,I10)

end subroutine cst_trajectory_write_snapshot

!===============================================================================
! Subroutine:  cst_trajectory_close
!===============================================================================

subroutine cst_trajectory_close

    use pmf_constants
    use pmf_dat
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    close(CST_TRJ)

    return

end subroutine cst_trajectory_close

!===============================================================================

end module cst_trajectory

