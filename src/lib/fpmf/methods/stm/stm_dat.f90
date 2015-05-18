!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module stm_dat

use pmf_sizes
use pmf_dat

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer                 :: fmode            ! 0 - disable STM, 1 - enabled STM
integer                 :: fsample          ! output sample period in steps
integer                 :: ftensor          ! 0 - unity, 1 - normal, 2 - massweighted
character(PMF_MAX_PATH) :: fbeadidfile      ! name of file containing bead id

! server part ------------------------------------------------------------------
character(PMF_MAX_PATH) :: fserverkey           ! stm-server key file name
character(PMF_MAX_PATH) :: fserver              ! stm-server name
character(PMF_MAX_PATH) :: fpassword            ! stm-server password

! stm server -----------------
integer                 :: bead_id              ! bead ID
integer                 :: client_id            ! stm client ID
logical                 :: use_key              ! is stm-server enabled?

! item list --------------------------------------------------------------------
type CVTypeSTM
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data

    real(PMFDP)             :: startvalue       ! start value
    real(PMFDP)             :: stopvalue        ! stop value
    real(PMFDP)             :: target_value     ! required value of restraint

    real(PMFDP)             :: force_constant   ! force constant
    real(PMFDP)             :: deviation        ! deviation between real and actual value
    real(PMFDP)             :: energy           ! restraint energy
end type CVTypeSTM

! ----------------------
integer                     :: NumOfSTMCVs      ! number of CVs in a group
type(CVTypeSTM),allocatable :: STMCVList(:)     ! definition of CVs
real(PMFDP)                 :: TotalSTMEnergy   ! total restraints energy

! ----------------------
! Modes:
integer,parameter :: BMO_UNKNOWN           = 0
integer,parameter :: BMO_INITIALIZATION    = 1
integer,parameter :: BMO_ACCUMULATION      = 2
integer,parameter :: BMO_EQUILIBRATION     = 3
integer,parameter :: BMO_PRODUCTION        = 4
integer,parameter :: BMO_WAITFORRENDEZVOUS = 5
integer,parameter :: BMO_TERMINATE         = 6

integer                     :: stmmode          ! current STM mode
integer                     :: stmsteps         ! number of steps for current mode
integer                     :: curstep          ! current step
real(PMFDP),allocatable     :: beadpos(:)       ! bead position
real(PMFDP),allocatable     :: pmf(:)           ! mean force
real(PMFDP),allocatable     :: MTZ(:,:)         ! MTZ matrix

! ------------------------------------------------------------------------------

end module stm_dat

