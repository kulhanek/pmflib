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

module rst_dat

use pmf_sizes
use pmf_constants
use pmf_cvs
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable RST, 1 - enabled RST, 2 - enabled TI, 3 - enabled SM
integer     :: fsample      ! output sample period in steps
integer     :: fplevel      ! print level
logical     :: frestart     ! restart job with previous data
integer     :: fhistupdate  ! how often is restart file written
integer     :: fhistclear   ! after first 'fhistclear' steps RST and US-ABF accumulators will be erased
                            ! will be reset (default 0)
integer     :: fsamplefreq  ! how often take samples
real(PMFDP) :: fwarnlevel   ! warning level for restraint energy


! item list --------------------------------------------------------------------
type CVTypeUM
    integer                 :: cvindx           ! CV index
    class(CVType),pointer   :: cv               ! cv data

    character(PMF_MAX_MODE) :: mode             ! mode - constant (C)
                                                !      - incremental (I)
                                                !      - change to value (V)
                                                !      - wall restraints (W)
                                                !      - controlled steering (S)
    real(PMFDP)             :: startvalue       ! start value
    real(PMFDP)             :: stopvalue        ! stop value
    real(PMFDP)             :: target_value     ! required value of restraint

    real(PMFDP)             :: left_value       ! left value for wall restraint
    real(PMFDP)             :: right_value      ! right value for wall restraint

    real(PMFDP)             :: force_constant   ! sigma value
    real(PMFDP)             :: deviation        ! deviation between real and actual value
    real(PMFDP)             :: energy           ! restraint energy

    real(PMFDP)             :: force            ! restraint force

    real(PMFDP)             :: min_value        ! left range for this coordinate
    real(PMFDP)             :: max_value        ! right range for this coordinate
    integer                 :: nbins            ! numbins for this coordinate

    logical                 :: set_value         ! set value from initial value

    real(PMFDP),pointer     :: control_values(:) ! values for controlled streering
end type CVTypeUM

! RST accumulator --------------------------------------------------------------

type,extends(PMFAccuType) ::  UMAccuType

     integer,pointer                :: nsamples(:)  ! number of hits into bins
end type UMAccuType

! global data ------------------------------------------------------------------

! ----------------------
! regular umbrella sampling
integer                     :: NumOfRSTCVs          ! number of restraints
type(CVTypeUM),allocatable  :: RSTCVList(:)         ! input definition of restraints
real(PMFDP)                 :: TotalRstEnergy       ! total restraints energy

! ----------------------
! umbrella integration
type(UMAccuType)            :: rstaccu              ! accumulated data for TI
integer                     :: insidesamples
integer                     :: outsidesamples

! ----------------------
! value accumulation
integer                     :: faccumulation        ! total number of accumulated steps
real(PMFDP),allocatable     :: svalues(:)           ! accumulated values
real(PMFDP),allocatable     :: s2values(:)          ! accumulated square of values

!===============================================================================

end module rst_dat

