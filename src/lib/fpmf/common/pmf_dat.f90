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
logical                 :: fprint_inpcrds       ! print input coordinates including atom
                                                ! names and residue names
logical                 :: fprint_masks         ! print atom masks in CV definitions
logical                 :: fenable_pbc          ! enable PBC condition (partially)
logical                 :: fenable_hessian      ! enable the calculation of second derivatives
logical                 :: fmonitor_paths       ! enable path monitoring

! parallel setup ---------------------------------------------------------------
logical     :: fmaster      = .true.    ! .true. for master process
integer     :: fmytaskid    = 0         ! process id

! local copies of MD variables -------------------------------------------------
integer     :: fnatoms       ! total number of atoms
integer     :: fnstlim      ! length of simulation in steps
integer     :: fstep        ! current MD step
integer     :: fsystype     ! system type (NT,NTV,NTP)
real(PMFDP) :: fdt          ! dt of step in [fs]
real(PMFDP) :: ftime        ! actual time in [fs]
real(PMFDP) :: ftemp        ! simulation temperature in [K]
integer     :: fintalg      ! integration algorithm

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
logical                     :: con_enabled
logical                     :: rst_enabled
logical                     :: mtd_enabled
logical                     :: abf_enabled
logical                     :: mon_enabled
logical                     :: remd_enabled
logical                     :: stm_enabled
logical                     :: pdrv_enabled
logical                     :: gap_enabled

! MASTER =======================================================================

! local atoms ------------------------------------------------------------------
real(PMFDP),allocatable     :: InitialCrd(:,:)  ! initial coordinates of local atoms
real(PMFDP),allocatable     :: Mass(:)          ! local atom masses
real(PMFDP),allocatable     :: MassInv(:)       ! mass inverse
real(PMFDP),allocatable     :: Crd(:,:)         ! current coordinates
real(PMFDP),allocatable     :: Frc(:,:)         ! current system forces
real(PMFDP),allocatable     :: Vel(:,:)         ! current system velocities
real(PMFDP),allocatable     :: DelV(:,:)        ! current system -forces
real(PMFDP)                 :: PotEne           ! current system potential energy
type(CVContextType)         :: CVContext        ! current CV context (values and derivatives)

! used by Blue moon
real(PMFDP),allocatable     :: CrdP(:,:)        ! coordinates
real(PMFDP),allocatable     :: VelP(:,:)        ! velocities
type(CVContextType)         :: CVContextP

! MASTER =======================================================================

! file names -------------------------------------------------------------------
character(PMF_MAX_PATH)     :: fcvsdef      = '{CVS}'
character(PMF_MAX_PATH)     :: fpathsdef    = '{PATHS}'

! constraint dynamics -----------------------------
character(PMF_MAX_PATH)     :: fcondef      = '{CON}'
character(PMF_MAX_PATH)     :: fconout      = '_con.out'
character(PMF_MAX_PATH)     :: fconrst      = '_con.rst'
character(PMF_MAX_PATH)     :: fconctr

! restraint dynamics -------------------------------
character(PMF_MAX_PATH)     :: frstdef      = '{RST}'
character(PMF_MAX_PATH)     :: frstout      = '_rst.out'
character(PMF_MAX_PATH)     :: frsthist     = '_rst.hist'
character(PMF_MAX_PATH)     :: frstctr

! metadynamics -----------------------------------
character(PMF_MAX_PATH)     :: fmtddef     = '{MTD}'
character(PMF_MAX_PATH)     :: fmtdout     = '_mtd.out'
character(PMF_MAX_PATH)     :: fmtdrst     = '_mtd.rst'
character(PMF_MAX_PATH)     :: fmtdcvs     = '_mtd.cvs'
character(PMF_MAX_PATH)     :: fmtdhills   = '_mtd.hills'
character(PMF_MAX_PATH)     :: fmtdgpout   = '_mtd.gp'

! adaptive biasing force method ------------------
character(PMF_MAX_PATH)     :: fabfdef      = '{ABF}'
character(PMF_MAX_PATH)     :: fabfmask     = '_abf.mask'
character(PMF_MAX_PATH)     :: fabfout      = '_abf.out'
character(PMF_MAX_PATH)     :: fabfrst      = '_abf.rst'
character(PMF_MAX_PATH)     :: fabftrj      = '_abf.trj'
character(PMF_MAX_PATH)     :: fabfgpout    = '_abf.gpout'

! string method ----------------------------------
character(PMF_MAX_PATH)     :: fstmdef      = '{STM}'
character(PMF_MAX_PATH)     :: fstmout      = '_stm.out'

! monitoring -------------------------------------
character(PMF_MAX_PATH)     :: fmondef      = '{MON}'
character(PMF_MAX_PATH)     :: fmonout      = '_mon.out'

! path driving -----------------------------------
character(PMF_MAX_PATH)     :: fpdrvdef     = '{PDRV}'
character(PMF_MAX_PATH)     :: fpdrvout     = '_pdrv.out'

! remd -------------------------------------------
character(PMF_MAX_PATH)     :: fremdout     = '_remd.out'

! gap --------------------------------------------
character(PMF_MAX_PATH)     :: fgapdef      = '{GAP}'
character(PMF_MAX_PATH)     :: fgapout      = '_gap.out'

!-------------------------------------------------------------------------------

end module pmf_dat

