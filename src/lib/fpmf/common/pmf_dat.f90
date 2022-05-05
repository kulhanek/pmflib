!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_dat

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

! MASTER =======================================================================

! [control] --------------------------------------------------------------------
character(PMF_MAX_PATH) :: ftopology            ! fake topology
logical                 :: fdebug               ! more verbose output
logical                 :: frepmpifrag          ! report MPI fragmentation (pmemd-new, only master)
logical                 :: fprint_inpcrds       ! print input coordinates including atom
                                                ! names and residue names
logical                 :: fprint_masks         ! print atom masks in CV definitions
logical                 :: fenable_pbc          ! enable PBC condition (partially)
logical                 :: fmonitor_paths       ! enable path monitoring

! parallel setup ---------------------------------------------------------------
logical     :: fmaster      = .true.    ! .true. for master process
integer     :: fmytaskid    = 0         ! process id
integer     :: fnumoftasks  = 0         ! number of tasks

! local copies of MD variables -------------------------------------------------
integer     :: fnatoms      ! total number of atoms
integer     :: fnstlim      ! length of simulation in steps
integer     :: fstep        ! current MD step
integer     :: fsystype     ! system type (NT,NTV,NTP)
real(PMFDP) :: fdt          ! dt of step in [fs]
real(PMFDP) :: ftime        ! actual time in [fs]
real(PMFDP) :: ftemp        ! simulation temperature in [K]
integer     :: fintalg      ! integration algorithm
logical     :: fshake       ! SHAKE or different constraint algorithm in use in main MD engine

real(PMFDP) :: fdtx         ! timestep in internal units
real(PMFDP) :: ifdtx        ! reciprocal timestep

real(PMFDP),allocatable     :: frmass(:)     ! atom masses of all atoms, should be initialized in pmf_xxxx_end_init

! status variables -------------------------------------------------------------
logical     :: fcanexmdloop ! client know how to exit md loop
integer     :: fexit_mdloop ! non zero value should terminate MD loop

! PBC support ------------------------------------------------------------------
integer     :: fbox_type    ! box type
real(PMFDP) :: fucell(3,3)  ! direct lattice vectors
real(PMFDP) :: frecip(3,3)  ! reciprocal lattice vectors
real(PMFDP) :: fbox_volume  ! box volume
real(PMFDP) :: fbox_sphere  ! maximum radius of sphere inscribed to cell

! MASTER/SLAVE =================================================================

! translation data -------------------------------------------------------------
integer                     :: NumOfLAtoms      ! number of local atoms in CVs
integer,allocatable         :: RIndexes(:)      ! translation from local (PMFLib)
                                                ! to global (external code) atom indexes

! methods ----------------------------------------------------------------------
logical                     :: pmf_enabled      ! if any method is on
logical                     :: cst_enabled
logical                     :: rst_enabled
logical                     :: mtd_enabled
logical                     :: abf_enabled
logical                     :: abp_enabled
logical                     :: mon_enabled
logical                     :: stm_enabled
logical                     :: pdrv_enabled
! testing
logical                     :: tabf_enabled
logical                     :: usabf_enabled

! requests
logical                     :: shake_force_required     ! we need force due to SHAKE
logical                     :: lng_force_required       ! collect Langevin forces

! MASTER =======================================================================

character(PMF_KEYLINE)      :: DriverName

! local atoms ------------------------------------------------------------------
real(PMFDP),allocatable     :: InitialCrd(:,:)  ! initial coordinates of local atoms
real(PMFDP),allocatable     :: Mass(:)          ! local atom masses
real(PMFDP),allocatable     :: MassInv(:)       ! mass inverse
real(PMFDP),allocatable     :: Crd(:,:)         ! current coordinates in t
real(PMFDP),allocatable     :: Frc(:,:)         ! current system forces in t due to potential energy
real(PMFDP),allocatable     :: Vel(:,:)         ! current system velocities in t-dt/2
! FIXME - better description
real(PMFDP)                 :: KinEne           ! current system kinetic energy in t-dt
real(PMFDP)                 :: PotEne           ! current system potential energy in t
real(PMFDP)                 :: PMFEne           ! current PMFLib potential energy in t (from RST, MTD, STM)
type(CVContextType)         :: CVContext        ! current CV context (values and derivatives) in t

! used by Blue moon
real(PMFDP),allocatable     :: CrdP(:,:)        ! coordinates in t+dt
real(PMFDP),allocatable     :: VelP(:,:)        ! velocities
type(CVContextType)         :: CVContextP

real(PMFDP),allocatable     :: CrdBar(:,:)      ! coordinates in t+dt without SHAKE
real(PMFDP),allocatable     :: SHAKEFrc(:,:)    ! current SHAKE forces in t(+dt) - after PotForce

real(PMFDP)                 :: LNG_c_implic     ! Langevin setup from AMBER
real(PMFDP)                 :: LNG_c_explic
real(PMFDP),allocatable     :: LNGFrc(:,:)      ! current Langevin forces in t(+dt) - after PotForce

! MASTER =======================================================================

! file names -------------------------------------------------------------------
character(PMF_MAX_PATH)     :: fcvsdef
character(PMF_MAX_PATH)     :: fpathsdef

! constraint dynamics -----------------------------
character(PMF_MAX_PATH)     :: fcstdef
character(PMF_MAX_PATH)     :: fcstout
character(PMF_MAX_PATH)     :: fcstrst
character(PMF_MAX_PATH)     :: fcstfrst
character(PMF_MAX_PATH)     :: fcstctr
character(PMF_MAX_PATH)     :: fcsttrj

! restraint dynamics -------------------------------
character(PMF_MAX_PATH)     :: frstdef
character(PMF_MAX_PATH)     :: frstout
character(PMF_MAX_PATH)     :: frsthist
character(PMF_MAX_PATH)     :: frstctr

! metadynamics -----------------------------------
character(PMF_MAX_PATH)     :: fmtddef
character(PMF_MAX_PATH)     :: fmtdout
character(PMF_MAX_PATH)     :: fmtdrst
character(PMF_MAX_PATH)     :: fmtdtrj
character(PMF_MAX_PATH)     :: fmtdhills

! adaptive biasing force method ------------------
character(PMF_MAX_PATH)     :: fabfdef
character(PMF_MAX_PATH)     :: fabfmask
character(PMF_MAX_PATH)     :: fabfout
character(PMF_MAX_PATH)     :: fabfrst
character(PMF_MAX_PATH)     :: fabftrj

! adaptive biasing force method (testing) --------
character(PMF_MAX_PATH)     :: ftabfdef
character(PMF_MAX_PATH)     :: ftabfout
character(PMF_MAX_PATH)     :: ftabficf
character(PMF_MAX_PATH)     :: ftabfrst
character(PMF_MAX_PATH)     :: ftabftrj

! umbrella sampling/adaptive biasing force method (testing) --------
character(PMF_MAX_PATH)     :: fusabfdef
character(PMF_MAX_PATH)     :: fusabfout
character(PMF_MAX_PATH)     :: fusabfrst
character(PMF_MAX_PATH)     :: fusabftrj

! adaptive biasing potential method ------------------
character(PMF_MAX_PATH)     :: fabpdef
character(PMF_MAX_PATH)     :: fabpout
character(PMF_MAX_PATH)     :: fabprst
character(PMF_MAX_PATH)     :: fabptrj

! string method ----------------------------------
character(PMF_MAX_PATH)     :: fstmdef
character(PMF_MAX_PATH)     :: fstmout

! monitoring -------------------------------------
character(PMF_MAX_PATH)     :: fmondef
character(PMF_MAX_PATH)     :: fmonout

! path driving -----------------------------------
character(PMF_MAX_PATH)     :: fpdrvdef
character(PMF_MAX_PATH)     :: fpdrvout

!-------------------------------------------------------------------------------

end module pmf_dat

