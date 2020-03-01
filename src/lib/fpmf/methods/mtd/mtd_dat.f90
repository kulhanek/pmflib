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

module mtd_dat

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode            ! mode of metadynamics: 
                                    ! 0 - disable MTD
                                    ! 1 - conventional MTD
                                    ! 2 - grid MTD
integer     :: fplevel          ! output print level
integer     :: fpsample         ! output print sampling
integer     :: fmetastep        ! step size of metadynamics
real(PMFDP) :: fheight          ! height of gaussian
real(PMFDP) :: fmetatemp        ! temperature for well-tempered mtd
integer     :: fmetavary        ! what to vary when well-tempered mtd is on
                                    ! 0 - nothing
                                    ! 1 - height
                                    ! 2 - step
logical     :: frestart         ! restart job with previous data
integer     :: fextout          ! control extended output (cvs,hills)
                                    ! 0- none
                                    ! 1 - new
                                    ! 2 - append
integer     :: fbuffersize      ! buffer size
integer     :: fscaling         ! scaling type: 0 - none, 1 - step signal, 2 - ramp signal
integer     :: fdelaying        ! delaying time in md steps
integer     :: ftransition      ! transition time in md steps

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
    real(PMFDP)             :: max_dist         ! maximum distance to consider a hill

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
    integer,pointer             :: steps(:)
end type MTDHistType

type(MTDHistType),pointer   :: hill_history

! grid list ---------------------------------------------------------------------

type CVInfoTypeMTD
    integer                  :: nbins           ! number of grid bins
    real(PMFDP)              :: min_value       ! left boundary of coordinate
    real(PMFDP)              :: max_value       ! left boundary of coordinate
    real(PMFDP)              :: bin_width       ! (right-left)/numbins
    real(PMFDP)              :: width           ! right - left
end type CVInfoTypeMTD

type MTDGridType
    integer                      :: tot_cvs           ! total number of independent CVs
    type(CVInfoTypeMTD), pointer :: sizes(:)          ! grid informations
    integer                      :: tot_nbins         ! number of total grids
    real(PMFDP),pointer          :: binpositions(:,:) ! position of grids

    ! MTD potential
    real(PMFDP),pointer          :: mtdpotential(:)   ! MTD potential
    real(PMFDP),pointer          :: mtdforce(:,:)     ! MTD forces   
end type MTDGridType

type(MTDGridType)           :: grid

! determine if data should be exchanged with multiple-walker server
logical                     :: fserverupdate = .false.

! how many times it was restarted or client ID in multiple-walker approach
integer                     :: frstlevel        = 0
integer                     :: cvs_starts       = 0
integer                     :: meta_step        = 0

! when fmetastep is varied use this counter
integer                     :: meta_next_fstep  = 0

!===============================================================================

end module mtd_dat

