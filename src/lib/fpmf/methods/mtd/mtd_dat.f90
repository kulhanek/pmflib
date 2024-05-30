!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module mtd_dat

use pmf_sizes
use pmf_constants
use pmf_cvs
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode            ! mode of metadynamics:
                                    ! 0 - disable MTD
                                    ! 1 - grid MTD
integer     :: fsample          ! output sample period in steps
integer     :: fmetastep        ! step size of metadynamics
real(PMFDP) :: fheight          ! height of gaussian
real(PMFDP) :: fmetatemp        ! temperature for well-tempered mtd, if 0 K MTD-WT is disabled
logical     :: frestart         ! restart job with previous data
integer     :: frstupdate       ! how often is restart file written
integer     :: ftrjsample       ! how often save accumulator to "accumulator evolution"
logical     :: fwritehills      ! record deposition of gaussians
logical     :: fswitch2zero     ! switch to zero at discretionary boundary

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails

integer                 :: failure_counter
integer                 :: client_id

! CV list ----------------------------------------------------------------------

type CVTypeMTD
    integer                 :: cvindx           ! CV index
    class(CVType),pointer   :: cv               ! cv data

    real(PMFDP)             :: min_value        ! min value for recorded region
    real(PMFDP)             :: max_value        ! max value for recorded region

    real(PMFDP)             :: min_deposit      ! min value for deposited region
    real(PMFDP)             :: max_deposit      ! max value for deposited region

    real(PMFDP)             :: width            ! gaussian width
    integer                 :: nbins            ! number of bins (for mtd-energy)

    real(PMFDP)             :: buffer           ! switch biasing potential to zero at CV boundary
end type CVTypeMTD

integer                     :: NumOfMTDCVs              ! number of CVs
type(CVTypeMTD),allocatable :: MTDCVList(:)             ! list of CVs

! ------------------------------------------------------------------------------

type,extends(PMFAccuType) :: MTDGridType
    ! accumulated data
    real(PMFDP),pointer             :: binpos(:,:)              ! position of grids
    integer,pointer                 :: nsamples(:)              ! number of hits into bins
    real(PMFDP),pointer             :: mtdpot(:)                ! MTD potential
    real(PMFDP),pointer             :: mtdforce(:,:)            ! MTD forces
    real(PMFDP),pointer             :: aplforce(:,:)            ! MTD applied forces

    real(PMFDP),pointer             :: widths(:)
    real(PMFDP),pointer             :: iwidths2(:)

    ! incremental data
    integer,pointer                 :: inc_nsamples(:)          ! number of hits into bins
    real(PMFDP),pointer             :: inc_mtdpot(:)            ! MTD potential
    real(PMFDP),pointer             :: inc_mtdforce(:,:)        ! MTD forces
end type MTDGridType

! ----------------------
type(MTDGridType)           :: mtdaccu
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

real(PMFDP),allocatable     :: CVValues(:)          ! current CV values
real(PMFDP)                 :: TotalMTDEnergy       ! imposed MTD energy
real(PMFDP),allocatable     :: MTDForce(:)          ! imposed MTD force

! helpers
real(PMFDP),allocatable     :: sfac(:)              ! switching factors

! run-time information
integer                     :: meta_next_fstep  = 0
integer                     :: numofhills       = 0

!===============================================================================

end module mtd_dat

