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

module pmfdyn_restraints_types

use pmf_sizes

implicit none
!-------------------------------------------------------------------------------

type RSTRSEQ_TYPE
    integer             ::  i,j
    real(PMFDP)         ::  fk
    integer             ::  ih
    integer             ::  to_centre ! flag for restraining to geometry or mass centre
    integer             ::  change
    real(PMFDP)         ::  fk_final
end type RSTRSEQ_TYPE

!-------------------------------------------------------------------------------

end module pmfdyn_restraints_types
