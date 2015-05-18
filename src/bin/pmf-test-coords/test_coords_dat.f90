! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

module test_coords_dat

use pmf_sizes
use pmf_constants
use pmf_dat
use prmfile

implicit none

!===============================================================================

character(len=PMF_MAX_PATH) :: ControlFile
type(PRMFILE_TYPE)          :: ControlPrmfile
character(len=PMF_MAX_PATH) :: CoordFile
type(PRMFILE_TYPE)          :: CoordPrmfile
character(len=PMF_MAX_PATH) :: SnapshotFile
logical                     :: UseExternalSnapshosts = .false.

!===============================================================================

integer,parameter       :: TEST_DERIV       = 1     ! numerical versus analytical derivatives
integer,parameter       :: TEST_IDENTITY    = 2     ! identity between two coordinates
integer,parameter       :: TEST_VALUE       = 3     ! value comparison

! test setups ------------------------------------
integer                 :: test_type        = TEST_DERIV
integer                 :: max_atoms        = 10
real(PMFDP)             :: max_radius       = 20.0d0    ! maximum radius for coordinates
real(PMFDP)             :: max_mass         = 25.0d0    ! max mass of atom
integer                 :: num_of_tests     = 10
real(PMFDP)             :: num_diff         = 0.0001d0   ! difference for numeric differentiation
real(PMFDP)             :: alarm_treshold   = 0.00001d0   ! treshold for signaling problem
logical                 :: stop_when_error  = .false.
logical                 :: verbose_print    = .false.
real(PMFDP)             :: uniform_mass     = 0.0d0     ! zero means random masses

! results ----------------------------------------
integer                 :: total_num_of_tests = 0
integer                 :: num_of_failed_tests = 0

! data -------------------------------------------
class(CVType),pointer   :: rc1              ! reaction coordinate
class(CVType),pointer   :: rc2              ! reaction coordinate


real(PMFDP),allocatable         :: lx(:,:)          ! position in time t
real(PMFDP),allocatable         :: loc_x(:,:)       ! copy of positions for numerical
real(PMFDP),allocatable         :: fdx_ana(:,:)     ! first derivative of coordinate
real(PMFDP),allocatable         :: fdx_num(:,:)     ! first derivative of coordinate
type(CVContextType)             :: tmp_context
character(len=3),allocatable    :: symbols(:)       ! symbols

!===============================================================================

end module test_coords_dat
