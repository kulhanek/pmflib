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
    call pmf_timers_finalize(do_profiling)

end subroutine pmf_finalize_all

!===============================================================================
! Subroutine: pmf_finalize_pmf
!===============================================================================

subroutine pmf_finalize_pmf

    use pmf_dat
    use pmf_cvs

    implicit none
    integer     :: i
    ! ----------------------------------------------------------------------------

    if( allocated(RIndexes) ) then
        deallocate(RIndexes)
    end if
    if( allocated(InitialCrd) ) then
        deallocate(InitialCrd)
    end if
    if( allocated(Mass) ) then
        deallocate(Mass)
    end if
    if( allocated(MassInv) ) then
        deallocate(MassInv)
    end if
    if( allocated(Crd) ) then
        deallocate(Crd)
    end if
    if( allocated(Frc) ) then
        deallocate(Frc)
    end if
    if( allocated(Vel) ) then
        deallocate(Vel)
    end if
    if( allocated(DelV) ) then
        deallocate(DelV)
    end if
    if( allocated(CrdP) ) then
        deallocate(CrdP)
    end if
    if( allocated(VelP) ) then
        deallocate(VelP)
    end if

    ! destroy CVs
    do i=1,NumOfCVs
        call CVList(i)%cv%free_cv()
        deallocate(CVList(i)%cv)
    end do
    NumOfCVs = 0
    if( allocated(CVList) ) then
        deallocate(CVList)
    end if

end subroutine pmf_finalize_pmf

!===============================================================================
! Subroutine: pmf_finalize_methods
!===============================================================================

subroutine pmf_finalize_methods

    use pmf_dat
    use mon_finalize
    use rst_finalize
    use abf_finalize
    use abp_finalize
    use mtd_finalize
    use cst_finalize
    use pdrv_finalize
    use stm_finalize

    implicit none
    ! --------------------------------------------------------------------------

    if( abf_enabled ) then
        call abf_finalize_method
    end if

    if( abp_enabled ) then
        call abp_finalize_method
    end if

    if( cst_enabled ) then
        call cst_finalize_method
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
