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

module pmfdyn_thermostat_dat

use pmf_sizes
use pmf_constants

implicit none

! ------------------------------------------------------------------------------
! thermostat types
integer, parameter      :: THERMOSTAT_NONE          = 0
integer, parameter      :: THERMOSTAT_BERENDSEN     = 1
integer, parameter      :: THERMOSTAT_LANGEVIN      = 2
integer, parameter      :: THERMOSTAT_CHEATHAM      = 3

! [controls] ===================================================================
real(PMFDP)             :: Tmaxw                = -1.0d0
real(PMFDP)             :: Temp0                = 300.0d0
real(PMFDP)             :: Temp1                = -1.0d0
real(PMFDP)             :: Tau_T                = 400.0d0        ! in fs
real(PMFDP)             :: Gamma_T              = 0.001d0       ! in fs^-1
integer                 :: ThermostatType       = THERMOSTAT_BERENDSEN
integer                 :: Iseed                = 314159265
real(PMFDP)             :: HotAtomTreshold      = 1000.0d0
logical                 :: MonitorHotAtoms      = .true.
logical                 :: RemoveHotAtoms       = .false.

! ------------------------------------------------------------------------------

real(PMFDP)             :: NOF = 0.d0    ! total number of freedom
real(PMFDP)             :: TempT = 0.d0  ! target temperature
real(PMFDP)             :: TempA = 0.d0  ! temperature of whole system (Ndegf)

real(PMFDP)             :: HotAtomMaxEkin = 0.d0    ! maximal kinetic energy for "cold" atoms
real(PMFDP)             :: lgam_c = 0.d0 ! used by Langevin thermostat
real(PMFDP)             :: lgam_p = 0.d0 ! used by Langevin thermostat
real(PMFDP)             :: lgam_m = 0.d0 ! used by Langevin thermostat

! ------------------------------------------------------------------------------

end module pmfdyn_thermostat_dat
