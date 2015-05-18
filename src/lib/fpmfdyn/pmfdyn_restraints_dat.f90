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

module pmfdyn_restraints_dat

use pmf_sizes
use pmfdyn_restraints_types

implicit none

! ------------------------------------------------------------------------------
! types of sequential restraints
integer, parameter              :: SEQRES_FIX_TO_XYZ    = 0
integer, parameter              :: SEQRES_FIX_TO_COG    = 1
integer, parameter              :: SEQRES_FIX_TO_COM    = 2

! [sequence_restraints] ========================================================
integer                         ::  nrstr_seq   = 0
type(RSTRSEQ_TYPE), allocatable ::  rstseq(:)

! ------------------------------------------------------------------------------

end module pmfdyn_restraints_dat
