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
logical     :: frecord      ! record time progress

real(PMFDP) :: fepotaverage
real(PMFDP) :: fekinaverage

! linear ramp
integer     :: fhramp_min
integer     :: fhramp_max

! kernel smoothing
integer     :: fsmooth_kernel
logical     :: fswitch2zero

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails
integer                 :: fmwamode             ! 0 - all is transferred
                                                ! 1 - only MICF is transferred

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
    real(PMFDP)             :: wfac             ! smoothing factor in number of bins
    real(PMFDP)             :: buffer           ! switch biasing potential to zero at CV boundary
end type CVTypeABF

! ----------------------
integer                     :: NumOfABFCVs          ! number of ALL CVs in a group
type(CVTypeABF),allocatable :: ABFCVList(:)         ! definition of CVs

! ----------------------

type,extends(PMFAccuType) :: ABFAccuType

    real(PMFDP),pointer    :: binpos(:,:)               ! position of grids
    real(PMFDP),pointer    :: weights(:)                ! mask weights
    real(PMFDP),pointer    :: nsamples(:)               ! number of hits into bins

! free energy
    real(PMFDP),pointer    :: micf(:,:)                 ! mean ICF
    real(PMFDP),pointer    :: m2icf(:,:)                ! M2 of ICF

! enthalpy
    real(PMFDP),pointer    :: mepot(:)                  ! mean of pot energy
    real(PMFDP),pointer    :: m2epot(:)                 ! M2 of pot energy
    real(PMFDP),pointer    :: merst(:)                  ! mean of rst energy
    real(PMFDP),pointer    :: m2erst(:)                 ! M2 of rst energy
    real(PMFDP),pointer    :: mekin(:)                  ! mean of kin energy
    real(PMFDP),pointer    :: m2ekin(:)                 ! M2 of kin energy

! entropy
    real(PMFDP),pointer    :: metot(:)                  ! mean of tot energy
    real(PMFDP),pointer    :: m2etot(:)                 ! M2 of tot energy
    real(PMFDP),pointer    :: mpp(:,:)                  ! mean of tot energy + icf
    real(PMFDP),pointer    :: m2pp(:,:)                 ! M2 of tot energy + icf
    real(PMFDP),pointer    :: mpn(:,:)                  ! mean of tot energy - icf
    real(PMFDP),pointer    :: m2pn(:,:)                 ! M2 of tot energy - icf

! time recording for post-processing
    real(PMFDP),pointer    :: tcvs(:,:)
    real(PMFDP),pointer    :: tzinv(:,:,:)
    real(PMFDP),pointer    :: tbicf(:,:)
    real(PMFDP),pointer    :: tepot(:)
    real(PMFDP),pointer    :: terst(:)
    real(PMFDP),pointer    :: tekin(:)

! applied ICF - this is stored in accu but ignored
    real(PMFDP),pointer    :: bnsamples(:)              ! number of hits into bins
    real(PMFDP),pointer    :: bmicf(:,:)                ! applied MICF

! ABF force - incremental part for ABF-server
    real(PMFDP),pointer    :: inc_nsamples(:)           ! number of hits into bins
    real(PMFDP),pointer    :: inc_micf(:,:)             ! accumulated mean ICF
    real(PMFDP),pointer    :: inc_m2icf(:,:)            ! accumulated M2 of ICF
end type ABFAccuType

! ----------------------
type(ABFAccuType)           :: abfaccu                  ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

! global variables for abf - results -------------------------------------------
real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix              in t
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix   in t
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)              ! for LU decomposition
real(PMFDP),allocatable     :: fzinv0(:,:)          ! inverse of Z matrix   in t-dt

! determinant of mass metric tensors
real(PMFDP)                 :: fzdet
real(PMFDP)                 :: fzdetall

! helper arrays -------
real(PMFDP),allocatable     :: la(:)
real(PMFDP),allocatable     :: cvave(:)
real(PMFDP),allocatable     :: pxi0(:)
real(PMFDP),allocatable     :: pxi1(:)
real(PMFDP),allocatable     :: pxip(:)
real(PMFDP),allocatable     :: pxim(:)
real(PMFDP),allocatable     :: sfac(:)              ! switching factors

! ------------------------------------------------------------------------------

integer                     :: hist_len
real(PMFDP),allocatable     :: cvhist(:,:)          ! history of CV values (nCVS,hist_len)
real(PMFDP),allocatable     :: fhist(:,:,:)         ! history of forces
real(PMFDP),allocatable     :: xhist(:,:,:)         ! history of velocities
real(PMFDP),allocatable     :: vhist(:,:,:)         ! history of velocities
real(PMFDP),allocatable     :: fzinvhist(:,:,:)     ! history of fzinv
real(PMFDP),allocatable     :: xvelhist(:,:)        ! history of CV velocities
real(PMFDP),allocatable     :: xphist(:,:)          ! history of CV momenta
real(PMFDP),allocatable     :: icfhist(:,:)         ! history of ICF forces
real(PMFDP),allocatable     :: zdhist(:,:,:,:)      ! history of ZD
real(PMFDP),allocatable     :: micfhist(:,:)        ! history of ABF bias
real(PMFDP),allocatable     :: epothist(:)          ! history of Epot
real(PMFDP),allocatable     :: ersthist(:)          ! history of Erst
real(PMFDP),allocatable     :: ekinhist(:)          ! history of Ekin

! GPR facility -----------------------------------------------------------------
integer                     :: gpr_len          ! MUST be odd number

! CVS GPR
real(PMFDP)                 :: gpr_cvs_width    ! kernel width in fs
integer                     :: gpr_cvs_kernel   ! 0 - Exponential
                                                ! 1 - MC(3/2)
                                                ! 2 - MC(5/2)
                                                ! 3 - ARDSE
                                                ! 4 - Epanechnikov (parabolic)
                                                ! 5 - Quartic (biweight)
                                                ! 6 - Triweight
real(PMFDP)                 :: gpr_cvs_noise    ! noise magnitude

! ICF GPR
logical                     :: gpr_icf_cdf      ! use central differences for the second differentiation in ICF calc.
real(PMFDP)                 :: gpr_icf_width    ! kernel width in fs
integer                     :: gpr_icf_kernel   ! kernel type
real(PMFDP)                 :: gpr_icf_noise    ! noise magnitude

! ENE GPR
integer                     :: gpr_ene_smooth   ! 0 - no smoothing
                                                ! 1 - smooth etot
                                                ! 2 - smooth ekin
real(PMFDP)                 :: gpr_ene_width    ! kernel width in fs
integer                     :: gpr_ene_kernel   ! kernel type
real(PMFDP)                 :: gpr_ene_noise    ! noise magnitude

integer                     :: gpr_boundary     ! skip analysis at boundaries of gpr_len
integer                     :: gpr_rank         ! rank for SVD inversion
real(PMFDP)                 :: gpr_rcond        ! rcond for automatic rank determination

real(PMFDP),allocatable     :: gpr_K_cvs(:,:)   ! co-variance matrix for CVS
real(PMFDP),allocatable     :: gpr_K_icf(:,:)   ! co-variance matrix for ICF
real(PMFDP),allocatable     :: gpr_K_ene(:,:)   ! co-variance matrix for ENE
real(PMFDP),allocatable     :: gpr_data(:)      ! GPR input data
real(PMFDP),allocatable     :: gpr_model(:)     ! GPR model
real(PMFDP),allocatable     :: gpr_kff_cvs(:,:) ! inference for cvs
real(PMFDP),allocatable     :: gpr_kfd_cvs(:,:) ! inference for cvs
real(PMFDP),allocatable     :: gpr_kff_icf(:,:) ! inference for icf
real(PMFDP),allocatable     :: gpr_kfd_icf(:,:) ! inference for icf
real(PMFDP),allocatable     :: gpr_kff_ene(:,:) ! inference for ene

! SG setup ---------------------------------------------------------------------
integer     :: fsgframelen
integer     :: fsgorder

real(PMFDP),allocatable     :: sg_c0(:)
real(PMFDP),allocatable     :: sg_c1(:)
real(PMFDP),allocatable     :: sg_c2(:)

! ------------------------------------------------------------------------------

! smoothing facility
integer                     :: max_snb_size
integer,allocatable         :: snb_list(:,:)
real(PMFDP),allocatable     :: sweights(:)

! ------------------------------------------------------------------------------

end module abf_dat

