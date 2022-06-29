!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_kinene

use pmf_sizes

implicit none

!-------------------------------------------------------------------------------

! energy blob
type PMFKineticEnergy
    real(PMFDP)                 :: KinEneVV         ! velocity verlet kinetic energy in t-dt
    real(PMFDP)                 :: KinEneHA         ! velocity verlet kinetic energy in t-dt
    real(PMFDP)                 :: KinEneLF         ! leapfrog kinetic energy in t-dt/2
    real(PMFDP)                 :: KinEneV3         ! accurate kinetic energy in t-dt
    real(PMFDP)                 :: KinEneV4         ! accurate kinetic energy in t-2*dt
    real(PMFDP)                 :: KinEneV5         ! accurate kinetic energy in t-2*dt
    real(PMFDP)                 :: KinEneV6         ! accurate kinetic energy in t-3*dt
end type PMFKineticEnergy

!-------------------------------------------------------------------------------

end module pmf_kinene

