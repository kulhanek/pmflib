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
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable ABF
                            ! 1 - enable ABF (simple ABF algorithm - 2p)
                            ! 2 - enable ABF (original ABF algorithm - 4p)
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! 1 - restart job with previous data, 0 - otherwise not
integer     :: frstupdate   ! how often is restart file written
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
logical     :: fapply_mask  ! off - disable ABF mask, on - enable ABF mask
logical     :: fapply_abf   ! on - apply ABF, off - do not apply ABF
integer     :: feimode      ! interpolation/extrapolation mode
                            ! 0 - disabled
                            ! 1 - linear ramp

logical     :: fenthalpy    ! collect data for enthalpy calculation
logical     :: fentropy     ! collect data for entropy calculation

logical     :: fsmoothetot  ! smooth etot prior covariance calculation

real(PMFDP) :: fepotaverage
real(PMFDP) :: fekinaverage

! linear ramp
integer     :: fhramp_min
integer     :: fhramp_max

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails
integer                 :: fmwamode             ! 0 - all is transferred
                                                ! 1 - only MICF is transfered

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
    real(PMFDP)             :: wfac             ! smoothing factor
end type CVTypeABF

! ----------------------
integer                     :: NumOfABFCVs      ! number of CVs in a group
type(CVTypeABF),allocatable :: ABFCVList(:)     ! definition of CVs

! ----------------------

type,extends(PMFAccuType) :: ABFAccuType

    ! biasing force manipulation
    real(PMFDP),pointer    :: weights(:)                ! mask weights

    ! MICF
    integer,pointer        :: nsamples(:)               ! number of hits into bins
    real(PMFDP),pointer    :: micf(:,:)                 ! mean ICF
    real(PMFDP),pointer    :: m2icf(:,:)                ! M2 of ICF
    real(PMFDP),pointer    :: mepot(:)                  ! mean of pot energy
    real(PMFDP),pointer    :: m2epot(:)                 ! M2 of pot energy
    real(PMFDP),pointer    :: metot(:)                  ! mean of total energy
    real(PMFDP),pointer    :: m2etot(:)                 ! M2 of total energy
    real(PMFDP),pointer    :: c11hh(:,:)                ! cov(H,H) for entropy

    ! applied ICF - this is not stored in PMF accumulator
    integer,pointer        :: bnsamples(:)              ! number of hits into bins
    real(PMFDP),pointer    :: bmicf(:,:)                ! applied MICF

    ! ABF force - incremental part for ABF-server
    integer,pointer        :: inc_nsamples(:)           ! number of hits into bins
    real(PMFDP),pointer    :: inc_micf(:,:)             ! accumulated mean ICF
    real(PMFDP),pointer    :: inc_m2icf(:,:)            ! accumulated M2 of ICF
    real(PMFDP),pointer    :: inc_mepot(:)              ! accumulated mean of pot energy
    real(PMFDP),pointer    :: inc_m2epot(:)             ! accumulated M2 of pot energy
    real(PMFDP),pointer    :: inc_metot(:)              ! accumulated mean of total energy
    real(PMFDP),pointer    :: inc_m2etot(:)             ! accumulated M2 of total energy
    real(PMFDP),pointer    :: inc_c11hh(:,:)            ! accumulated cov(H,H) for entropy
end type ABFAccuType

! ----------------------
type(ABFAccuType)           :: abfaccu                  ! accumulated forces
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
real(PMFDP),allocatable     :: v0(:,:)        ! velocity in previous step

real(PMFDP),allocatable     :: la(:)          ! ABF force in coordinate direction
real(PMFDP),allocatable     :: zd0(:,:,:)     ! ZD0
real(PMFDP),allocatable     :: zd1(:,:,:)     ! ZD1
real(PMFDP),allocatable     :: pxi0(:)        !
real(PMFDP),allocatable     :: pxi1(:)        !
real(PMFDP),allocatable     :: pxip(:)        !
real(PMFDP),allocatable     :: pxim(:)        !
real(PMFDP),allocatable     :: pdum(:)        !
real(PMFDP),allocatable     :: avg_values(:)  ! average values of coordinates at t - 3/2dt

integer                     :: hist_len
real(PMFDP),allocatable     :: cvhist(:,:)      ! history of CV values
real(PMFDP),allocatable     :: pcvhist(:,:)     ! history of CV momenta
real(PMFDP),allocatable     :: epothist(:)      ! history of Epot
real(PMFDP),allocatable     :: etothist(:)      ! history of Etot

integer                     :: gpr_len          ! MUST be odd number
real(PMFDP)                 :: gpr_width        ! SE width time steps
real(PMFDP)                 :: gpr_noise        !
integer                     :: gpr_kernel       ! 0 - MC(3/2)
                                                ! 1 - MC(5/2)
                                                ! 2 - ARDSE
real(PMFDP),allocatable     :: gpr_K(:,:)       ! covariance matrix
real(PMFDP),allocatable     :: gpr_model(:)     ! GPR model
real(PMFDP),allocatable     :: gpr_kff(:)
real(PMFDP),allocatable     :: gpr_kdf(:)
integer,allocatable         :: gpr_indx(:)
integer                     :: gpr_info

! ------------------------------------------------------------------------------

end module abf_dat

