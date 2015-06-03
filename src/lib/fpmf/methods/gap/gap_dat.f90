!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module gap_dat

use pmf_sizes
use pmf_cvs
use gp_dat_mod

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable GAP, 1 - enabled GAP
integer     :: fsample      ! output sample period in steps
integer     :: fplevel      ! print level

! item list --------------------------------------------------------------------
type CVTypeGAP
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    character(PMF_MAX_PATH) :: gpfile           ! gp file name
    character(PMF_MAX_PATH) :: groupname        ! group name
    integer                 :: indexid          ! index id
    integer                 :: groupid          ! group id
    real(PMFDP)             :: min_value        ! minimum value
    real(PMFDP)             :: max_value        ! maximum value
end type CVTypeGAP

! ----------------------

integer                     :: NumOfGAPCVs           ! number of CVs used in GAP
type(CVTypeGAP),allocatable :: GAPCVList(:)          ! GAP CV items

! ----------------------

integer                     :: NumOfGAPGroups        ! number of GAP groups
type GroupTypeGAP
    type(gp_basic)          :: gp                    ! gp of GAP
    character(PMF_MAX_PATH) :: groupname             ! group name
    integer                 :: nindexes              ! number of indexes of GAP type
    integer,pointer         :: gapcvindx(:)          ! local GAP cv index
    real(PMFDP),pointer     :: values(:)             ! cv values
    real(PMFDP)             :: energy                ! GAP energy
    real(PMFDP),pointer     :: forces(:)             ! cv forces
    real(PMFDP),pointer     :: min_values(:)         ! minimum values
    real(PMFDP),pointer     :: max_values(:)         ! maximum values
end type GroupTypeGAP
type(GroupTypeGAP),allocatable   :: GAPGroupList(:)  ! GAP group items

real(PMFDP)                 :: TotalGAPEnergy        ! total imposed energy

! ------------------------------------------------------------------------------

end module gap_dat
