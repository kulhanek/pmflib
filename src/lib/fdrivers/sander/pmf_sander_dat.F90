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

module pmf_sander_dat

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

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_sander_dat
