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

module abp_dat

use pmf_sizes
use pmf_dat
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable ABP
                            ! 1 - enabled ABP (direct mode)
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! 1 - restart job with previous data, 0 - otherwise not
integer     :: frstupdate   ! how often is restart file written
integer     :: feimode      ! extrapolation / interpolation mode
                            ! 1 - linear ramp
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
real(PMFDP) :: fhbias       ! bias height


! linear ramp mode I
integer     :: fhramp       ! ramp size

character(PMF_MAX_PATH) :: fkernel  ! SE - squared exponential

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abp-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abp-server key file name
character(PMF_MAX_PATH) :: fserver              ! abp-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails

! abp server -----------------
integer                 :: client_id            ! abp walker client ID
integer                 :: failure_counter      ! current number of MWA failures


! item list --------------------------------------------------------------------
type CVTypeABP
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    integer                 :: set              ! coordinate set index
    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins
    real(PMFDP)             :: width            ! width factor
end type CVTypeABP

! ----------------------
integer                     :: NumOfABPCVs      ! number of CVs in a group
type(CVTypeABP),allocatable :: ABPCVList(:)     ! definition of CVs

! ----------------------

type,extends(PMFAccuType)   :: ABPAccuType
    ! ABP data
    real(PMFDP)                     :: M            ! max value of pop
    integer,pointer                 :: nsamples(:)  !
    real(PMFDP),pointer             :: dpop(:,:)    !
    real(PMFDP),pointer             :: pop(:)       ! without 1.0
    real(PMFDP),pointer             :: widths(:)    ! widths
    real(PMFDP),pointer             :: iwidths2(:)  ! inverse width squares

    ! helper data
    real(PMFDP),pointer             :: binpos(:,:)  ! CV values in bin centers
    real(PMFDP)                     :: dnorm

    ! ABP data - incremental part for ABP-server
    integer,pointer                 :: inc_nsamples(:)
    real(PMFDP),pointer             :: inc_dpop(:,:)
    real(PMFDP),pointer             :: inc_pop(:)
end type ABPAccuType

! ----------------------
type(ABPAccuType)           :: abpaccu          ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples

! ----------------------
real(PMFDP),allocatable     :: cvvalues(:)      ! CV values
real(PMFDP),allocatable     :: la(:)            ! ABP force in coordinate direction

real(PMFDP)                 :: cfac
real(PMFDP)                 :: kt
real(PMFDP),allocatable     :: diffvalues(:)    ! grid-current considering periodicity

! ------------------------------------------------------------------------------

end module abp_dat

