!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2010,2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
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

module pmf_sander_dat_d01

use prmfile_dat
use pmf_sizes

implicit none

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

character(len=PMF_MAX_PATH)     :: ControlFileName
type(PRMFILE_TYPE)              :: ControlPrmfile

#ifdef MPI
integer,allocatable             :: atm_owner_map(:) ! atom map among processes
real(PMFDP),allocatable         :: tmp_a(:,:)   ! helper array
real(PMFDP),allocatable         :: tmp_b(:,:)   ! helper array
real(PMFDP),allocatable         :: tmp_c(:,:)   ! helper array
#endif

! interface binding check
integer, parameter              :: PMFLIB_CHECK_INT1 = 1089523658
real(PMFDP), parameter          :: PMFLIB_CHECK_R81  = 1.78493547
character(len=10), parameter    :: PMFLIB_CHECK_STR1 = 'PMFLib v06'
character(len=10), parameter    :: PMFLIB_CHECK_STR2 = 'DRVABI d01'

! energy array
integer, parameter              :: PMFLIB_EKIN_VV               = 1
integer, parameter              :: PMFLIB_EKIN_LF               = 2
integer, parameter              :: PMFLIB_EKIN_HA               = 3
integer, parameter              :: PMFLIB_EKIN_SIZE             = PMFLIB_EKIN_HA

! setup array
integer, parameter              :: PMFLIB_SETUP_ISCHEME         = 1
integer, parameter              :: PMFLIB_SETUP_SIZE            = PMFLIB_SETUP_ISCHEME

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_sander_dat_d01
