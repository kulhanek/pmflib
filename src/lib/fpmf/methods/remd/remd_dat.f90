!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module remd_dat

use pmf_sizes

implicit none

! This implementation of replica-exchange molecular dynamics uses Epot from previous
! time step to decide if replicas should be exchanged or not
! This approximation should not cause serious problems if MD time step
! is not too large

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer                 :: fmode        ! 0 - disable REMD, 1 - enabled REMD
integer                 :: fsample      ! sampling frequency
character(PMF_MAX_PATH) :: fserverkey   ! remd-server key file name
character(PMF_MAX_PATH) :: fserver      ! remd-server name
character(PMF_MAX_PATH) :: fpassword    ! remd-server password

! remd server -----------------
logical                 :: use_key          ! is remd-server key enabled?
integer                 :: ReplicaId        ! replica ID
integer                 :: BathId           ! bath ID

! internal data ----------------------------------------------------------------
logical         :: REMDFullSwap ! temp or full swap exchange mode?
integer         :: REMDPeriod   ! exchange data period
real(PMFDP)     :: OldBathTemp  ! old system temperature
real(PMFDP)     :: CurBathTemp  ! current system temperature

! ------------------------------------------------------------------------------

end module remd_dat

