!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf2_finalize

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf2_finalize_method
!===============================================================================

subroutine abf2_finalize_method

    use abf2_client
    use abf2_trajectory
    use abf2_restart
    use abf2_output

    implicit none
    ! --------------------------------------------------------------------------

    call abf2_client_unregister
    call abf2_trajectory_close
    call abf2_restart_write
    call abf2_output_close

end subroutine abf2_finalize_method

!===============================================================================

end module abf2_finalize
