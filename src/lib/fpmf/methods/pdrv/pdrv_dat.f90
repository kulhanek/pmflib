!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pdrv_dat

use pmf_sizes
use pmf_paths

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable PDRV, 1 - enabled PDRV
integer     :: fsample      ! output sample period in steps

! item list --------------------------------------------------------------------
type CVTypePDRV
    integer                 :: pathindx         ! path index
    class(PathType),pointer :: path             ! path data
    character(PMF_MAX_MODE) :: mode             ! mode - attach to value (A)
                                                !      - incremental (I)
                                                !      - change to value (P)
    logical                 :: set_attach_to_value
    logical                 :: initial_value_set
    real(PMFDP)             :: initial_alpha    ! initial alpha
    real(PMFDP)             :: final_alpha      ! final alpha
    real(PMFDP)             :: req_alpha        ! requested alpha
    real(PMFDP),allocatable :: ipos(:)          ! initial CV position for attach_to mode
    real(PMFDP),allocatable :: fpos(:)          ! final CV position for attach_to mode
end type CVTypePDRV

! ----------------------
integer                      :: NumOfPDRVItems  ! number of monitored CVs
type(CVTypePDRV),allocatable :: PDRVCVList(:)   ! monitored items

! ------------------------------------------------------------------------------

end module pdrv_dat

