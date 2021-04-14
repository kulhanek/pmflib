!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module abf_dat

use pmf_sizes
use pmf_dat

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable ABF
                            ! 1 - enable ABF (standard algorithm)
                            ! 2 - enable ABF (numerical algorithm)
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! 1 - restart job with previous data, 0 - otherwise not
integer     :: frstupdate   ! how often is restart file written
integer     :: feimode      ! extrapolation / interpolation mode
                            ! 1 - linear ramp I
                            ! 2 - linear ramp II
                            ! 3 - block averages
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
integer     :: fmask_mode   ! 0 - disable ABF mask, 1 - enable ABF mask
logical     :: fapply_abf   ! on - apply ABF, off - do not apply ABF
logical     :: fprint_icf   ! T - print instantaneous collective forces (icf), F - do not print
logical     :: fcache_icf   ! T - cache icf into memory and dump them at the end,
                            ! F - write icf immediately at each time step
logical     :: frawicf      ! T - use raw icf data (in internal units), F - transform them to user req. units

logical     :: faccuepot    ! accumulate PotEne
logical     :: faccuekin    ! accumulate KinEne
real(PMFDP) :: ftotoffset

! linear ramp mode I (feimode .eq. 1)
integer     :: fhramp       ! ramp size

! linear ramp mode II (feimode .eq. 2)
integer     :: fhramp_min   ! min value of linear ramp - ABF force is ignored below this value
integer     :: fhramp_max   ! max  value of linear ramp

! block everages (feimode .eq. 3)
integer     :: fblock_size  ! block size for ABF force pre-accumulation

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails

! abf server -----------------
integer                 :: client_id            ! abf walker client ID
integer                 :: failure_counter      ! current number of MWA failures

! item list --------------------------------------------------------------------
type CVTypeABF
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    integer                 :: set              ! coordinate set index
    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins
end type CVTypeABF

! ----------------------
integer                     :: NumOfABFCVs      ! number of CVs in a group
type(CVTypeABF),allocatable :: ABFCVList(:)     ! definition of CVs

! accu types -------------------------------------------------------------------
type CVInfoTypeABF
    integer                  :: nbins           ! number of accumulator bins
    real(PMFDP)              :: min_value       ! left boundary of coordinate
    real(PMFDP)              :: max_value       ! left boundary of coordinate
    real(PMFDP)              :: bin_width       ! (right-left)/numbins
    real(PMFDP)              :: width           ! right - left
end type CVInfoTypeABF

! ----------------------

type ABFAccuType
     integer                       :: tot_cvs       ! total number of independent CVs
     type(CVInfoTypeABF), pointer  :: sizes(:)      ! CV information
     integer                       :: tot_nbins     ! number of total bins

     ! ABF force
     real(PMFDP),pointer    :: weights(:)               ! mask weights
     integer,pointer        :: nsamples(:)              ! number of hits into bins
     real(PMFDP),pointer    :: icfsum(:,:)              ! accumulated ABF force
     real(PMFDP),pointer    :: icfsum2(:,:)             ! accumulated square of ABF force
     real(PMFDP),pointer    :: etotsum(:)               ! accumulated potential energy
     real(PMFDP),pointer    :: etotsum2(:)              ! accumulated square of potential energy
     real(PMFDP),pointer    :: icfetotsum(:,:)          ! accumulated icfsum * etotsum
     real(PMFDP),pointer    :: icfetotsum2(:,:)         ! accumulated square of icfsum * etotsum

     ! ABF force - incremental part for ABF-server
     integer,pointer        :: inc_nsamples(:)          ! number of hits into bins
     real(PMFDP),pointer    :: inc_icfsum(:,:)          ! accumulated ABF force
     real(PMFDP),pointer    :: inc_icfsum2(:,:)         ! accumulated square of ABF force
     real(PMFDP),pointer    :: inc_etotsum(:)           ! accumulated potential energy
     real(PMFDP),pointer    :: inc_etotsum2(:)          ! accumulated square of potential energy
     real(PMFDP),pointer    :: inc_icfetotsum(:,:)      ! accumulated icfsum * etotsum
     real(PMFDP),pointer    :: inc_icfetotsum2(:,:)     ! accumulated square of icfsum * etotsum

     ! ABF force - block pre-sampling
     integer,pointer        :: block_nsamples(:)        ! number of hits into bins
     real(PMFDP),pointer    :: block_icfsum(:,:)        ! accumulated ABF force
     real(PMFDP),pointer    :: block_etotsum(:)         ! accumulated PotEne
     real(PMFDP),pointer    :: block_icfetotsum(:,:)    ! accumulated square of icfsum * etotsum
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
real(PMFDP),allocatable     :: zd1(:,:,:)     ! ZD1
real(PMFDP),allocatable     :: pxi0(:)        !
real(PMFDP),allocatable     :: pxi1(:)        !
real(PMFDP),allocatable     :: pxip(:)        !
real(PMFDP),allocatable     :: pxim(:)        !
real(PMFDP),allocatable     :: avg_values(:)  ! average values of coordinates at t - 3/2dt

real(PMFDP),allocatable     :: cvaluehist0(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist1(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist2(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist3(:)   ! history of coordinate values

real(PMFDP)                 :: epothist0   ! history of Epot
real(PMFDP)                 :: epothist1   ! history of Epot
real(PMFDP)                 :: epothist2   ! history of Epot
real(PMFDP)                 :: epothist3   ! history of Epot

real(PMFDP)                 :: ekinhist0   ! history of Ekin
real(PMFDP)                 :: ekinhist1   ! history of Ekin
real(PMFDP)                 :: ekinhist2   ! history of Ekin
real(PMFDP)                 :: ekinhist3   ! history of Ekin

real(PMFDP), allocatable    :: icf_cache(:,:)   ! icf_cache(2*ncvs,fnstlim)

! ------------------------------------------------------------------------------

end module abf_dat

