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

module rst_dat

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable RST, 1 - enabled RST, 1 - enabled RST + acumulator
integer     :: fsample      ! output sample period in steps
integer     :: fplevel      ! print level
logical     :: frestart     ! restart job with previous data
integer     :: fhistupdate  ! how often is restart file written
integer     :: fhistclear   ! after first 'fhistclear' steps histogram
                            ! will be reset (default 0)

! item list --------------------------------------------------------------------
type CVTypeUM
    integer                 :: cvindx           ! CV index
    class(CVType),pointer   :: cv               ! cv data

    character(PMF_MAX_MODE) :: mode             ! mode - constant (C)
                                                !      - incremental (I)
                                                !      - change to value (V)
                                                !      - wall restraints (W)
    real(PMFDP)             :: startvalue       ! start value
    real(PMFDP)             :: stopvalue        ! stop value
    real(PMFDP)             :: target_value     ! required value of restraint

    real(PMFDP)             :: left_value       ! left value for wall restraint
    real(PMFDP)             :: right_value      ! right value for wall restraint

    real(PMFDP)             :: force_constant   ! sigma value
    real(PMFDP)             :: deviation        ! deviation between real and actual value
    real(PMFDP)             :: energy           ! restraint energy

    real(PMFDP)             :: min_value        ! left range for this coordinate
    real(PMFDP)             :: max_value        ! right range for this coordinate
    integer                 :: nbins            ! numbins for this coordinate

    logical                 :: set_value        ! set value from initial value
end type CVTypeUM

! accu types -------------------------------------------------------------------
type CVInfoTypeUM
    integer                 :: nbins            ! number of accumulator bins
    real(PMFDP)             :: min_value        ! left boundary of coordinate
    real(PMFDP)             :: max_value        ! left boundary of coordinate
    real(PMFDP)             :: bin_width        ! (right-left)/numbins
    real(PMFDP)             :: width            ! right - left
end type CVInfoTypeUM

! ----------------------
type UMAccuType
     integer                        :: tot_cvs      ! total number of independent CVs
     type(CVInfoTypeUM), pointer    :: sizes(:)     ! accumulator informations
     integer                        :: tot_nbins    ! number of total bins

     integer,pointer                :: nsamples(:)  ! number of hits into bins
end type UMAccuType

! global data ------------------------------------------------------------------

! ----------------------
! regular umbrella sampling
integer                     :: NumOfRSTItems        ! number of restraints
type(CVTypeUM),allocatable  :: RSTCVList(:)        ! input definition of restraints
real(PMFDP)                 :: TotalRstEnergy       ! total restraints energy

! ----------------------
! umbrella integration
type(UMAccuType)            :: Accumulator          ! accumulated data for TI
integer                     :: insidesamples
integer                     :: outsidesamples

!===============================================================================

end module rst_dat

