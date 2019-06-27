!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module pmf_constants

use pmf_sizes

implicit none

!===============================================================================

! INTERNAL UNITS are:
! time      - fs
! length    - A
! velocity  - ?????
! energy    - kcal/mol
! force     - kcal/(mol.A)
! gradient  - kcal/(mol.A)
! mass      - g/mol

!===============================================================================

! math part ------------------------------------------------------------------
real(PMFDP),parameter   :: PMF_PI       = 3.1415926535897932384626433832795d0
real(PMFDP),parameter   :: PMF_HPI      = 0.5d0 * PMF_PI
real(PMFDP),parameter   :: PMF_R2D      = 180.0d0 / PMF_PI
real(PMFDP),parameter   :: PMF_D2R      = PMF_PI / 180.0d0

! gas constant 8.314 472(15) J mol-1 K-1
real(PMFDP), parameter  :: PMF_Rgas     = 0.0019872065d0     ! kcal mol-1 K-1 = 8.314 472 / 4184

! mass * velocity^2 -> kcal/mol
! g/mol * A^2 / fs^2 = 0.001 kg/mol * (10^-10)^2 m^2 / ((10^-15)^s^2)
! 0.001*10^-20/10^-30 * kg/mol*m^2/s^s = 10^-3*10^10 J/mol = 10^7 J/mol
! 10^7 J/mol = 10^7 / 4184 kcal/mol
! dt (fs) = vdt * sqrt( 4184 / 10^7) = 0.02045482828087295384
real(PMFDP), parameter  :: PMF_DT2VDT   = 0.02045482828087295384d0 ! sqrt(pc_cal/1d4)
real(PMFDP), parameter  :: PMF_VDT2DT   = 1.0d0 / PMF_DT2VDT

! hartree -> kcal/mol
real(PMFDP), parameter  :: PMF_HARTREE2KCL = 627.5059d0
real(PMFDP), parameter  :: PMF_KCL2HARTREE = 1.0d0 / PMF_HARTREE2KCL

! kJ/mol -> kcal/mol
real(PMFDP), parameter  :: PMF_KJ2KCL   = 1.0d0 / 4.184d0
real(PMFDP), parameter  :: PMF_KCL2KJ   = 1.0d0 / PMF_KJ2KCL

! eV -> kcal/mol
real(PMFDP), parameter  :: PMF_eV2KCL   = 23.06035d0
real(PMFDP), parameter  :: PMF_KCL2eV   = 1.0d0 / PMF_eV2KCL

! lambda (g A^2/(mol fs^2 CV) -> kcal/mol/CV
! 1000 * (10^10)^2 / (10^15)^2 = 1000 * 10^20 / 10^30 = 10^-7
real(PMFDP), parameter  :: PMF_L2CL     = 1e7 / 4184.d0
real(PMFDP), parameter  :: PMF_CL2L     = 1.0d0 / PMF_L2CL

! atomic unit time to fs (2.418 884 326 505(16)×10-17s)
real(PMFDP), parameter  :: PMF_AU2FS    = 2.418884326505d-2
real(PMFDP), parameter  :: PMF_FS2AU    = 1.0d0 / PMF_AU2FS

! ps to fs
real(PMFDP), parameter  :: PMF_PS2FS    = 1000.0d0
real(PMFDP), parameter  :: PMF_FS2PS    = 1.0d0 / PMF_PS2FS

! atomic unit length to A 5.291 772 108(18)×10-11m
real(PMFDP), parameter  :: PMF_AU2A     = 5.291772108d-1
real(PMFDP), parameter  :: PMF_A2AU     = 1.0d0 / PMF_AU2A

! pm to A
real(PMFDP), parameter  :: PMF_PM2A     = 0.01d0
real(PMFDP), parameter  :: PMF_A2PM     = 1.0d0 / PMF_PM2A

! nm to A
real(PMFDP), parameter  :: PMF_NM2A     = 10.0d0
real(PMFDP), parameter  :: PMF_A2NM     = 1.0d0 / PMF_NM2A

! libatoms mass to g/mol
! real(dp), parameter :: ELECTRONMASS_GPERMOL =  5.48579903e-4_dp !% grams/mol
! real(dp), parameter :: ELEM_CHARGE = 1.60217653e-19_dp !% coulombs
! real(dp), parameter :: HARTREE = 27.2113961_dp !% eV
! real(dp), parameter :: RYDBERG = 0.5_dp*HARTREE !% eV
! real(dp), parameter :: BOHR = 0.529177249_dp !% Angstrom
! real(dp), parameter :: HBAR_EVSEC = 6.5821220e-16_dp !% hbar in eV seconds
! real(dp), parameter :: HBAR_AU = 1.0_dp              !% hbar in a.u.
! real(dp), parameter :: HBAR = (HBAR_EVSEC*1e-15_dp)    !% hbar in eV fs
! real(dp), parameter :: ONESECOND = 1e15_dp           !% 1 second in fs
! real(dp), parameter :: ONESECOND_AU = (1.0_dp/(HBAR_EVSEC/(HBAR_AU*HARTREE))) !% 1 second in a.u.
! real(dp), parameter :: AU_FS = (1.0_dp/ONESECOND_AU*ONESECOND) !% a.u. time in fs
! real(dp), parameter :: MASSCONVERT = (1.0_dp/ELECTRONMASS_GPERMOL*HARTREE*AU_FS*AU_FS/(BOHR*BOHR))
! AU_FS= 1.0/1.0_dp/(6.5821220e-16/(1.0*27.2113961))*1e15
! AU_FS = (6.5821220e-16/(27.2113961))*1e15
! AU_FS = (6.5821220e-1/(27.2113961)) = 0.02418884343828283033 ----> PMF_AUT2FS
! 1.0/5.48579903e-4*27.2113961*0.02418884343828283033*0.02418884343828283033/(0.529177249*0.529177249)
! MASSCONVERT = 103.6427217590198535
real(PMFDP), parameter  :: PMF_AMU2LIBATOMSM = 103.6427217590198535d0
real(PMFDP), parameter  :: PMF_LIBATOMSM2AMU = 1.0d0 / PMF_AMU2LIBATOMSM

! a.u. to g/mol (amu)
real(PMFDP), parameter  :: PMF_AU2AMU = 1.0d0 / 1822.88842718d0
real(PMFDP), parameter  :: PMF_AMU2AU = 1.0d0 / PMF_AU2AMU

! common part ------------------------------------------------------------------
integer,parameter       :: PMF_INP      = 5
integer,parameter       :: PMF_OUT      = 6
integer,parameter       :: PMF_TEST     = 1342
integer,parameter       :: PMF_DEBUG    = 1000

! blue moon part ---------------------------------------------------------------
integer,parameter       :: CON_INP       = 150
integer,parameter       :: CON_OUT       = 151
integer,parameter       :: CON_RST       = 152
integer,parameter       :: CON_CTR       = 153

! umbrella part ----------------------------------------------------------------
integer,parameter       :: RST_INP       = 160
integer,parameter       :: RST_OUT       = 161
integer,parameter       :: RST_RST       = 162
integer,parameter       :: RST_CTR       = 163

! abf part ---------------------------------------------------------------------
integer,parameter       :: ABF_INP      = 170
integer,parameter       :: ABF_OUT      = 171
integer,parameter       :: ABF_RST      = 172
integer,parameter       :: ABF_TRJ      = 173
integer,parameter       :: ABF_GPOUT    = 174
integer,parameter       :: ABF_IFC      = 175

! metadyn part ----------------------------------------------------------------
integer,parameter       :: MTD_INP      = 180
integer,parameter       :: MTD_OUT      = 181
integer,parameter       :: MTD_RST      = 182
integer,parameter       :: MTD_CVS      = 183
integer,parameter       :: MTD_HILLS    = 184
integer,parameter       :: MTD_GPOUT    = 185

! string method ----------------------------------------------------------------
integer,parameter       :: STM_INP      = 310
integer,parameter       :: STM_OUT      = 311
integer,parameter       :: STM_BEADID   = 312

! monitoring part --------------------------------------------------------------
integer,parameter       :: MON_INP      = 190
integer,parameter       :: MON_OUT      = 191

! remd part --------------------------------------------------------------------
integer,parameter       :: REMD_OUT     = 211

! gap part ---------------------------------------------------------------------
integer,parameter       :: GAP_INP      = 221
integer,parameter       :: GAP_OUT      = 222

! path driving method ----------------------------------------------------------
integer,parameter       :: PDRV_INP     = 412
integer,parameter       :: PDRV_OUT     = 413

! xyz file ---------------------------------------------------------------------
integer,parameter       :: PMF_XYZ      = 231
integer,parameter       :: PMF_XYZ_S    = 232   ! xyz stream
integer,parameter       :: PMF_EVEC     = 233

!===============================================================================

! sys types
integer,parameter       :: SYS_UNK              = -1
integer,parameter       :: SYS_NT               = 0
integer,parameter       :: SYS_NTV              = 1
integer,parameter       :: SYS_NTP              = 2

! box types
integer,parameter       :: BOX_ISOLATED_SYSTEM  = 0
integer,parameter       :: BOX_ORTHOGONAL       = 1
integer,parameter       :: BOX_GENERAL          = 2

!===============================================================================

integer,parameter       :: IA_LEAP_FROG         = 0
integer,parameter       :: IA_VEL_VERLET        = 1

!===============================================================================

! conditions
integer,parameter       :: CND_GT         = 1
integer,parameter       :: CND_GE         = 2
integer,parameter       :: CND_LT         = 3
integer,parameter       :: CND_LE         = 4
integer,parameter       :: CND_EQ         = 5
integer,parameter       :: CND_NE         = 6

!===============================================================================

end module pmf_constants
