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

module pmf_finalize

implicit none
contains

!===============================================================================
! Subroutine: pmf_finalize_all
!===============================================================================

subroutine pmf_finalize_all(do_profiling)

    use pmf_timers

    implicit none
    logical :: do_profiling
    ! ----------------------------------------------------------------------------

    call pmf_finalize_methods
    call pmf_finalize_pmf

    if( do_profiling ) then
        call pmf_timers_finalize
    end if

end subroutine pmf_finalize_all

!===============================================================================
! Subroutine: pmf_finalize_methods
!===============================================================================

subroutine pmf_finalize_pmf

    implicit none
    ! ----------------------------------------------------------------------------

    ! nothing to be here

end subroutine pmf_finalize_pmf

!===============================================================================
! Subroutine: pmf_finalize_methods
!===============================================================================

subroutine pmf_finalize_methods

    use pmf_dat
    use mon_finalize
    use rst_finalize
    use abf_finalize
    use mtd_finalize
    use con_finalize
    use pdrv_finalize
    use stm_finalize

    implicit none
    ! --------------------------------------------------------------------------

    if( abf_enabled ) then
        call abf_finalize_method
    end if

    if( con_enabled ) then
        call con_finalize_method
    end if

    if( mtd_enabled ) then
        call mtd_finalize_method
    end if

    if( mon_enabled ) then
        call mon_finalize_method
    end if

    if( rst_enabled ) then
        call rst_finalize_method
    end if

    if( pdrv_enabled ) then
        call pdrv_finalize_method
    end if

    if( stm_enabled ) then
        call stm_finalize_method
    end if

    return

end subroutine pmf_finalize_methods

!===============================================================================

end module pmf_finalize
