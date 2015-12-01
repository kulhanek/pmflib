! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

module pmfdyn_system_dat

use pmf_sizes
use pmf_constants
use prmfile
use smf_xyzfile_type
use smf_sizes

implicit none

! ------------------------------------------------------------------------------
! input files
character(len=PMF_MAX_PATH) :: ControlFile
type(PRMFILE_TYPE)          :: ControlPrmfile

character(len=PMF_MAX_PATH) :: InputRstFile
logical                     :: NeedRestart          = .false.
logical                     :: Restart              = .false.
logical                     :: OnlyInit             = .false.
character(len=PMF_MAX_PATH) :: OutputRstFile
character(len=PMF_MAX_PATH) :: OutputTrjFile

integer,parameter           :: IO_RST               = 324
integer,parameter           :: IO_TRJ               = 321
integer,parameter           :: IO_FRC               = 323

! [dynamics] ===================================================================
integer                     :: nsteps               = 1000
real(PMFDP)                 :: stepsize             = 1.0d0     ! in fs

! [intervals] ==================================================================
integer                     :: rst_freq             = 1000      ! how often actualize restart file
integer                     :: out_freq             = 100       ! output write freq
integer                     :: traj_freq            = 0         ! trajectory write freq
integer                     :: nscm_freq            = 0         ! removal of COM rot and trans

! ------------------------------------------------------------------------------
! files
type(XYZFILE_TYPE)          :: TrajectoryCrd

!-------------------------------------------------------------------------------

integer                                     :: natoms     = 0   ! number of atoms in the system
character(len=SMF_MAX_SYMBOL), allocatable  :: symbols(:)
real(PMFDP), allocatable                    :: mass(:)          ! system masses
real(PMFDP), allocatable                    :: winv(:)          ! system inverted masses
real(PMFDP), allocatable                    :: xtop(:,:)        ! initial system coordinates
logical, allocatable                        :: heavy(:)         ! heavy or light atom?

! system energies ----------------------
real(PMFDP)                 :: Etot
real(PMFDP)                 :: Ekin
real(PMFDP)                 :: Epot
real(PMFDP)                 :: Erst
real(PMFDP)                 :: Epmf

!-------------------------------------------------------------------------------

! general variables
real(PMFDP)                 :: dt                   = 0.0d0     ! in local vdt units
real(PMFDP)                 :: idt                  = 0.0d0
integer                     :: istep                = 0
real(PMFDP)                 :: time                 = 0.0d0     ! time from beginning in fs

! general arrays
real(PMFDP), allocatable    :: md_x(:,:)        ! system coordinates
real(PMFDP), allocatable    :: md_d(:,:)        ! system derivatives
real(PMFDP), allocatable    :: md_old_x(:,:)    ! for shake and for numerical forces
real(PMFDP), allocatable    :: md_v(:,:)
real(PMFDP), allocatable    :: md_old_v(:,:)    ! old velocities for thermostat

!-------------------------------------------------------------------------------

integer                     :: INITIALIZATION_TIMER     =  -2
integer                     :: CORE_TIMER               =  -3
integer                     :: FORCES_TIMER             =  -5
integer                     :: EXTERNAL_TIMER           =  -6
integer                     :: RE_TIMER                 =  -7
integer                     :: FINALIZATION_TIMER       =  -8

!-------------------------------------------------------------------------------

interface
	subroutine pot_ext_energy(x,d,ene)
    	real(8)            :: x(:,:)
    	real(8)            :: d(:,:)
    	real(8)            :: ene
	end subroutine
end interface

! ------------------------------------------------------------------------------

procedure(pot_ext_energy),pointer	::	pot_ext_energy_pts => NULL()

! ------------------------------------------------------------------------------

end module pmfdyn_system_dat
