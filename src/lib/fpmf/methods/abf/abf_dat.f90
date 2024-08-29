!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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
logical     :: fupdate_abf  ! on - update ABF bias, off - keep initial bias
integer     :: ficfsample   ! how often update ABF accumulator for ICF

! enthalpy/entropy calculations
logical     :: fenthalpy        ! collect data for enthalpy calculation
integer     :: fenthalpy_der    ! collect data for enthalpy derivative calculation
                                ! 0 - no
                                ! 1 - from forces
                                ! 2 - from velocities
logical     :: fentropy         ! collect data for entropy calculation
logical     :: fentdecomp       ! collect additional correlation terms
logical     :: ftds_add_bias    ! include ABF bias into TdS calculation
integer     :: ftds_ekin_src    ! source of kinetic energy, see abf_core_update_history_ene for supported values
real(PMFDP) :: fepotaverage
real(PMFDP) :: fekinaverage
integer     :: fenesample       ! how often update ABF accumulator for ENT and TDS
integer     :: finclude_pv      ! include pV term

integer     :: fepotsmooth
integer     :: ferstsmooth
integer     :: fekinsmooth

! US mode
logical     :: fusmode      ! enable US mode
logical     :: falignbias   ! move bottom of the biasing potential into the closest bin position

! interpolation modes
integer     :: feimode      ! interpolation/extrapolation mode
                            ! 0 - disabled
                            ! 1 - linear ramp
                            ! 2 - kernel smoother
                            ! 3 - linear smoother, only one CV

! linear ramp
integer     :: fhramp_min
integer     :: fhramp_max

! kernel smoothing
integer     :: fsmooth_kernel
logical     :: fswitch2zero

! 2p algorithm
integer     :: abf_p2_vx
integer     :: abf_p2_px

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is abf-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! abf-server key file name
character(PMF_MAX_PATH) :: fserver              ! abf-server name
integer                 :: fserverupdate        ! how often to communicate with server
integer                 :: fconrepeats          ! how many times to repeat connection
logical                 :: fabortonmwaerr       ! abort if communication with MWA fails
integer                 :: fmwamode             ! 0 - full update of MICF
                                                ! 1 - incremental update of MICF (default)

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

    ! US mode
    logical                 :: set_value        ! set target value to start value
    real(PMFDP)             :: target_value     ! required value of restraint
    real(PMFDP)             :: force_constant   ! sigma value
    real(PMFDP)             :: deviation        ! deviation between real and actual value
    real(PMFDP)             :: energy
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
    real(PMFDP),pointer    :: mgfx(:,:)                 ! mean -GFX
    real(PMFDP),pointer    :: m2gfx(:,:)                ! M2 of GFX

    real(PMFDP),pointer    :: msrdetz(:)                ! mean sqrt(det(Z))
    real(PMFDP),pointer    :: m2srdetz(:)               ! M2 of srdetz

    real(PMFDP),pointer    :: msrzii(:,:)                ! mean sqrt(Zii)
    real(PMFDP),pointer    :: m2srzii(:,:)               ! M2 of srzii

    real(PMFDP),pointer    :: mvol(:)                   ! mean volume
    real(PMFDP),pointer    :: m2vol(:)                  ! M2 of volume

! enthalpy & entropy
    real(PMFDP),pointer    :: ntds(:)                   ! number of hits into bins

    real(PMFDP),pointer    :: meint(:)                  ! mean of internal energy
    real(PMFDP),pointer    :: m2eint(:)                 ! M2 of internal energy
    real(PMFDP),pointer    :: mepot(:)                  ! mean of pot energy
    real(PMFDP),pointer    :: m2epot(:)                 ! M2 of pot energy
    real(PMFDP),pointer    :: merst(:)                  ! mean of rst energy
    real(PMFDP),pointer    :: m2erst(:)                 ! M2 of rst energy
    real(PMFDP),pointer    :: mekin(:)                  ! mean of kin energy
    real(PMFDP),pointer    :: m2ekin(:)                 ! M2 of kin energy
    real(PMFDP),pointer    :: mepv(:)                   ! mean of pV energy
    real(PMFDP),pointer    :: m2epv(:)                  ! M2 of pV energy

    real(PMFDP),pointer    :: metot(:)                  ! mean of tot energy
    real(PMFDP),pointer    :: m2etot(:)                 ! M2 of tot energy
    real(PMFDP),pointer    :: mpp(:,:)                  ! mean of tot energy + icf
    real(PMFDP),pointer    :: m2pp(:,:)                 ! M2 of tot energy + icf
    real(PMFDP),pointer    :: mpn(:,:)                  ! mean of tot energy - icf
    real(PMFDP),pointer    :: m2pn(:,:)                 ! M2 of tot energy - icf

    real(PMFDP),pointer    :: mpit(:,:)                 ! mean of ICF and tot energy product
    real(PMFDP),pointer    :: m2pit(:,:)                ! M2 of ICF and tot energy product

! enthalpy derivative
    real(PMFDP),pointer    :: micfp(:,:)                ! mean of ICF-P
    real(PMFDP),pointer    :: m2icfp(:,:)               ! M2 of internal energy

! entropy - decomposition
    real(PMFDP),pointer    :: mhicf(:,:)                ! mean of ICF - hamiltonian
    real(PMFDP),pointer    :: m2hicf(:,:)               ! M2 of ICF - hamiltonian
    real(PMFDP),pointer    :: mbicf(:,:)                ! mean of ICF - bias
    real(PMFDP),pointer    :: m2bicf(:,:)               ! M2 of ICF - bias

    real(PMFDP),pointer    :: c11hp(:,:)                ! co-variances
    real(PMFDP),pointer    :: c11hr(:,:)
    real(PMFDP),pointer    :: c11hk(:,:)
    real(PMFDP),pointer    :: c11hv(:,:)
    real(PMFDP),pointer    :: c11bp(:,:)
    real(PMFDP),pointer    :: c11br(:,:)
    real(PMFDP),pointer    :: c11bk(:,:)
    real(PMFDP),pointer    :: c11bv(:,:)

! applied ICF - this is stored in accu but ignored
    real(PMFDP),pointer    :: bnsamples(:)              ! number of hits into bins
    real(PMFDP),pointer    :: bmicf(:,:)                ! applied MICF

! ABF force - incremental part for ABF-server
    real(PMFDP),pointer    :: inc_nsamples(:)           ! number of hits into bins
    real(PMFDP),pointer    :: inc_micf(:,:)             ! accumulated mean ICF  - reverted indexing do to C/Fortran interoperability
    real(PMFDP),pointer    :: inc_m2icf(:,:)            ! accumulated M2 of ICF - reverted indexing do to C/Fortran interoperability
end type ABFAccuType

! ----------------------
type(ABFAccuType)           :: abfaccu                  ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

! global variables for abf - results -------------------------------------------
real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix              in t
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix   in t
real(PMFDP)                 :: fzdet
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)              ! for LU decomposition

! helper arrays -------
real(PMFDP),allocatable     :: la(:)
real(PMFDP),allocatable     :: pxia(:)
real(PMFDP),allocatable     :: pxif(:)
real(PMFDP),allocatable     :: picf(:)
real(PMFDP),allocatable     :: sfac(:)              ! switching factors
real(PMFDP),allocatable     :: vint(:,:)

! ------------------------------------------------------------------------------

integer                     :: hist_len
integer                     :: hist_fidx
integer                     :: fene_step


! it is not possible to buffer atom coordinates
! there are jumps due to nscm (removal of COM motion)

real(PMFDP),allocatable     :: cvhist(:,:)          ! history of CV values (nCVS,hist_len)
real(PMFDP),allocatable     :: vhist(:,:,:)         ! history of velocities
real(PMFDP),allocatable     :: cvderhist(:,:,:,:)   ! history of CV derivatives
real(PMFDP),allocatable     :: fzinvhist(:,:,:)     ! history of fzinv
real(PMFDP),allocatable     :: xphist(:,:)          ! history of CV momenta
real(PMFDP),allocatable     :: icfhist(:,:)         ! history of ABF ICF
real(PMFDP),allocatable     :: micfhist(:,:)        ! history of ABF bias

real(PMFDP),allocatable     :: fdetzhist(:)         ! history of fdetz
real(PMFDP),allocatable     :: fziihist(:,:)        ! history of fzii

real(PMFDP),allocatable     :: zdhist(:,:,:,:)      ! history of ZD
real(PMFDP),allocatable     :: fhist(:,:,:)         ! history of forces
real(PMFDP),allocatable     :: icfphist(:,:)        ! history of ABF ICF-P

real(PMFDP),allocatable     :: epothist(:)          ! history of Epot
real(PMFDP),allocatable     :: ersthist(:)          ! history of Erst
real(PMFDP),allocatable     :: ekinhist(:)          ! history of Ekin
real(PMFDP),allocatable     :: ekinlfhist(:)        ! history of EkinLF
real(PMFDP),allocatable     :: epvhist(:)           ! history of pV
real(PMFDP),allocatable     :: volhist(:)           ! history of volume
logical,allocatable         :: enevalidhist(:)      ! is energy valid?

! ------------------------------------------------------------------------------

! smoothing facility
integer                     :: max_snb_size
integer,allocatable         :: snb_list(:,:)
real(PMFDP),allocatable     :: sweights(:)

! ------------------------------------------------------------------------------

end module abf_dat

