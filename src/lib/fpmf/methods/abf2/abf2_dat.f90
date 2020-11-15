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

module abf2_dat

use pmf_sizes
use pmf_dat

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable ABF
                            ! 1 - enabled ABF (standard algorithm)
                            ! 2 - enabled ABF (numerical algorithm)
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! 1 - restart job with previous data, 0 - otherwise not
integer     :: frstupdate   ! how often is restart file written
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
integer     :: fmask_mode   ! 0 - disable ABF mask, 1 - enable ABF mask
logical     :: fapply_abf   ! on - apply ABF, off - do not apply ABF
integer     :: fblock_size  ! block size for ABF force pre-accumulation
integer     :: fintrpl      ! ABF force interpolation: 0 - no interpolation, 1 - linear interpolation
integer     :: fsmoothlimit = 500
integer     :: fsmoothfac   = 1

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
character(PMF_MAX_PATH) :: fpassword            ! abf-server password
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails

! abf server -----------------
integer                 :: client_id            ! abf walker client ID
logical                 :: use_key              ! is abf-server enabled?
integer                 :: failure_counter      ! current number of MWA failures

! item list --------------------------------------------------------------------
type CVTypeABF
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    integer                 :: set              ! coordinate set index
    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins
    real(PMFDP)             :: maxforce         ! maximum force to be aplied, zero = no treshold
end type CVTypeABF

! ----------------------
integer                     :: NumOfABFCVs      ! number of CVs in a group
type(CVTypeABF),allocatable :: ABFCVList(:)     ! definition of CVs

! accu types -------------------------------------------------------------------
type CVInfoTypeABF
    integer                  :: nbins           ! number of accumulator bins
    real(PMFDP)              :: min_value       ! left boundary of coordinate
    real(PMFDP)              :: max_value       ! right boundary of coordinate
    real(PMFDP)              :: bin_width       ! (max-min)/nbins
    real(PMFDP)              :: width           ! max-min
end type CVInfoTypeABF

! ----------------------

type ABFAccuType
     integer                       :: tot_cvs       ! total number of independent CVs
     type(CVInfoTypeABF), pointer  :: sizes(:)      ! accumulator information
     integer                       :: tot_nbins     ! number of total bins

     ! ABF force
     real(PMFDP),pointer    :: weights(:)               ! mask weights
     integer,pointer        :: nsamples(:)              ! number of hits into bins
     real(PMFDP),pointer    :: abfforce(:,:)            ! accumulated ABF force
     real(PMFDP),pointer    :: abfforce2(:,:)           ! accumulated square of ABF force

     ! ABF force - incremental part for ABF-server
     integer,pointer        :: nisamples(:)             ! number of hits into bins
     real(PMFDP),pointer    :: iabfforce(:,:)           ! accumulated ABF force
     real(PMFDP),pointer    :: iabfforce2(:,:)          ! accumulated square of ABF force
     real(PMFDP),pointer    :: iweights(:)              ! accumulated weights for each bin

     ! ABF force - block pre-sampling
     integer,pointer        :: block_nsamples(:)        ! number of hits into bins
     real(PMFDP),pointer    :: block_abfforce(:,:)      ! accumulated ABF force
end type ABFAccuType

! ----------------------
type(ABFAccuType)           :: accumulator          ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

! global variables for force calculation ---------------------------------------
real(PMFDP)                 :: fdtx                 ! timestep

! global variables for abf - results -------------------------------------------
real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)

! helper arrays -------
real(PMFDP),allocatable     :: a0(:,:)        ! acceleration from previous step (t-dt)
real(PMFDP),allocatable     :: a1(:,:)        ! acceleration in current step (t)

real(PMFDP),allocatable     :: la(:)          ! ABF force in coordinate direction
real(PMFDP),allocatable     :: zd0(:,:,:)     ! ZD0
real(PMFDP),allocatable     :: pxi0(:)        !
real(PMFDP),allocatable     :: pxi1(:)        !
real(PMFDP),allocatable     :: pxip(:)        !
real(PMFDP),allocatable     :: pxim(:)        !
real(PMFDP),allocatable     :: avg_values(:)  ! average values of coordinates at t - 3/2dt

real(PMFDP), allocatable    :: cvaluehist0(:)   ! history of coordinate values
real(PMFDP), allocatable    :: cvaluehist1(:)   ! history of coordinate values
real(PMFDP), allocatable    :: cvaluehist2(:)   ! history of coordinate values
real(PMFDP), allocatable    :: cvaluehist3(:)   ! history of coordinate values

! ------------------------------------------------------------------------------

end module abf2_dat

