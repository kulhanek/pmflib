!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module mtd_dat

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode            ! mode of metadynamics
integer     :: fsample          ! output print sampling
integer     :: fmetastep        ! step size of metadynamics
real(PMFDP) :: fheight          ! height of gaussian
logical     :: frestart         ! restart job with previous data
integer     :: fextout          ! control extended output (cvs,hills)
                                ! 0- none, 1 - new, 2 - append
integer     :: fbuffersize      ! buffer size

! server part ------------------------------------------------------------------
logical                 :: fserver_enabled      ! is metadyn-server enabled?
logical                 :: fserver_key_enabled  ! is metadyn-server enabled?
character(PMF_MAX_PATH) :: fserverkey           ! metadyn-server key file name
character(PMF_MAX_PATH) :: fserver              ! metadyn-server name
character(PMF_MAX_PATH) :: fpassword            ! metadyn-server password
integer                 :: fclient_id           ! client id

! item list --------------------------------------------------------------------

type CVTypeMTD
    integer                 :: cvindx           ! CV index
    class(CVType),pointer   :: cv               ! cv data

    real(PMFDP)             :: min_value        ! min value for wall rest.
    real(PMFDP)             :: max_value        ! max value for wall rest.
    real(PMFDP)             :: width            ! gaussian width
    integer                 :: nbins            ! number of bins (for mtd-energy)

    ! this applies to direct and extended version
    real(PMFDP)             :: meta_force       ! system forces
end type CVTypeMTD

integer                     :: NumOfMTDCVs              ! number of CVs
type(CVTypeMTD),allocatable :: MTDCVList(:)             ! list of CVs
real(PMFDP)                 :: TotalMTDEnergy           ! total imposed energy

! hills list --------------------------------------------------------------------

type MTDHistType
    type(MTDHistType),pointer   :: next_history_buffer
    integer                     :: length_of_buffer
    integer                     :: nrst_of_values
    real(PMFDP),pointer         :: values(:,:)
    real(PMFDP),pointer         :: widths(:,:)
    real(PMFDP),pointer         :: heights(:)
end type MTDHistType

type(MTDHistType),pointer   :: hill_history

! determine if data should be exchanged with multiple-walker server
logical                     :: fserverupdate = .false.

! how many times it was restarted or client ID in multiple-walker approach
integer                     :: frstlevel        = 0
integer                     :: cvs_starts       = 0
integer                     :: meta_step        = 0

!===============================================================================

end module mtd_dat

