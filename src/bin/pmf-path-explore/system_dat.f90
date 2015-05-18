! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_system_dat

use pmf_sizes
use pmf_constants
use prmfile
use smf_xyzfile_type
use smf_sizes
use pmf_paths

implicit none

! ------------------------------------------------------------------------------
! input files
character(len=PMF_MAX_PATH) :: ControlFile
type(PRMFILE_TYPE)          :: ControlPrmfile

character(len=PMF_MAX_PATH) :: InputFESFile
character(len=PMF_MAX_PATH) :: InputPath                = '{PATHS}'
character(len=PMF_MAX_PATH) :: OutputPathFile           = '_stm.path'
character(len=PMF_MAX_PATH) :: OutputPathSummaryFile    = '_stm.results'
character(len=PMF_MAX_PATH) :: TrajectoryFile           = '_stm.traj'

integer,parameter           :: IO_FES               = 235
integer,parameter           :: IO_PATH              = 236
integer,parameter           :: IO_SUM               = 337
integer,parameter           :: IO_TRJ               = 338

! [dynamics] ===================================================================
integer                     :: nsteps               = 1000
real(PMFDP)                 :: stepsize             = 0.0001    ! step size
real(PMFDP)                 :: smoothing            = 0.01      ! path smoothing

! [intervals] ==================================================================
integer                     :: print_freq           = 10        ! how often to update opt info
integer                     :: output_freq          = 10        ! how often to update result files
integer                     :: traj_freq            = 10        ! trajectory write freq
integer                     :: smooth_freq          = 1         ! smoothing frequency
integer                     :: reparam_freq         = 1         ! reparametrization frequency

!-------------------------------------------------------------------------------
! free energy surface

type coord_info_rec
    character(PMF_MAX_TYPE)     :: ctype           ! type of definition (bond,angle,..)
    character(PMF_MAX_CV_NAME)  :: name            ! name
    integer                     :: nbins           ! number of bins
    real(PMFDP)                 :: min_value       ! left boundary of coordinate
    real(PMFDP)                 :: max_value       ! left boundary of coordinate
    real(PMFDP)                 :: bin_width       ! (right-left)/numbins
    real(PMFDP)                 :: width           ! right - left
end type coord_info_rec

integer                             :: ncvs         ! number of coordinates
type(coord_info_rec), allocatable   :: cvs(:)       ! coordinates
integer                             :: tot_nbins    ! total number of bins

real(PMFDP),allocatable     :: fes(:)        ! free energy surface
real(PMFDP),allocatable     :: dfes(:,:)     ! free energy surface derivatives

!-------------------------------------------------------------------------------

! general variables
integer                     :: istep                = 0
integer                     :: nbeads
real(PMFDP), allocatable    :: values(:)            ! values along path
real(PMFDP), allocatable    :: deriv(:,:)           ! energy derivatives along path

class(PathType),pointer     :: path_item
real(PMFDP), allocatable    :: uposition(:,:)       ! updated positions
real(PMFDP), allocatable    :: sposition(:,:)       ! smoothed positions
real(PMFDP), allocatable    :: rposition(:,:)       ! reparametrized positions

real(PMFDP)                 :: old_path_len
real(PMFDP)                 :: path_len_change
real(PMFDP)                 :: MaxMovement
real(PMFDP)                 :: AveMovement
integer                     :: MaxMovementBead

!-------------------------------------------------------------------------------

integer :: INITIALIZATION_TIMER     =  -2
integer :: DERIVATIVES_TIMER        =  -3
integer :: CORE_TIMER               =  -4
integer :: FORCES_TIMER             =  -5
integer :: UPDATE_TIMER             =  -6
integer :: SMOOTH_TIMER             =  -7
integer :: REPARAM_TIMER            =  -8
integer :: FINALIZATION_TIMER       =  -9

! ------------------------------------------------------------------------------

end module pmf_system_dat
