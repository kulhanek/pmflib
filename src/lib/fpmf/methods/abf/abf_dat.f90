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
integer     :: fmode        ! 0 - disable ABF, 1 - enabled ABF
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! 1 - restart job with previous data, 0 - otherwise not
integer     :: frstupdate   ! how often is restart file written
integer     :: feimode      ! extrapolation / interpolation mode
                            ! 1 - linear ramp I, 2 - linear ramp II
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
integer     :: fmask_mode   ! 0 - disable ABF mask, 1 - enable ABF mask
logical     :: fapply_abf   ! on - apply ABF, off - do not apply ABF
logical     :: fprint_ifc   ! T - print instanteous forces, F - do not print
logical     :: fcache_ifc   ! T - cache ifc into memory and dump them at the end,
                            ! F - write ifc immediatelly at each time step
logical     :: frawifc      ! T - use raw ifc data (in internal units), F - transform them to user req. units

! linear ramp mode I (feimode .eq. 1)
integer     :: fhramp       ! ramp size

! linear ramp mode II (feimode .eq. 2)
integer     :: fhramp_min   ! min value of linear ramp - ABF force is ignored below this value
integer     :: fhramp_max   ! max  value of linear ramp

! gaussian process (feimode .eq. 3)
integer     :: fgpmin_samples   ! minimum number of samples per bin
integer     :: fgpmodel_update  ! how often the model is updated
integer     :: fgpprint_period  ! how often the interpolated data are printed to fabfgpout

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
character(PMF_MAX_PATH) :: fpassword            ! abf-server password
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection

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
    real(PMFDP)             :: switch           ! switch to zero width
    ! gaussian process data
    real(PMFDP)             :: fgplen           ! characteristic lengh-scale
    real(PMFDP)             :: fgpsigmaoffset   ! force sigma offset
    real(PMFDP)             :: fgpsigmafac      ! force sigma factor
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
     type(CVInfoTypeABF), pointer  :: sizes(:)      ! accumulator informations
     integer                       :: tot_nbins     ! number of total bins

     ! ABF force
     real(PMFDP),pointer    :: mask(:)                  ! mask weights
     integer,pointer        :: nsamples(:)              ! number of hits into bins
     real(PMFDP),pointer    :: abfforce(:,:)            ! accumulated ABF force
     real(PMFDP),pointer    :: abfforce2(:,:)           ! accumulated square of ABF force
     real(PMFDP),pointer    :: tweights(:)              ! total weights for each bin

     ! ABF force - incremental part for ABF-server
     integer,pointer        :: nisamples(:)             ! number of hits into bins
     real(PMFDP),pointer    :: iabfforce(:,:)           ! accumulated ABF force
     real(PMFDP),pointer    :: iabfforce2(:,:)          ! accumulated square of ABF force
     real(PMFDP),pointer    :: iweights(:)              ! accumulated weights for each bin
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

real(PMFDP), allocatable    :: ifc_cache(:,:)   ! ifc_cache(2*ncvs,fnstlim)

! gaussian process ----
integer                     :: gpmaxsize
integer                     :: gpsize
real(PMFDP),allocatable     :: xvalues(:)
real(PMFDP),allocatable     :: yvalues(:)
real(PMFDP),allocatable     :: svalues(:)
real(PMFDP),allocatable     :: gpyvalues(:)
real(PMFDP),allocatable     :: xcov(:,:)
real(PMFDP),allocatable     :: kstar(:)
real(PMFDP),allocatable     :: alpha(:)
integer,allocatable         :: gpindx(:)

! ------------------------------------------------------------------------------

end module abf_dat

