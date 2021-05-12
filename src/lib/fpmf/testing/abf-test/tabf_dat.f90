!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module tabf_dat

use pmf_sizes
use pmf_dat

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable ABF
                            ! 1 - enable ABF (standard algorithm - 4p)
                            ! 2 - enable ABF (numerical algorithm - 2p)
                            ! 3 - enable ABF (numerical algorithm - 2p, no SHAKE)
integer     :: fsample      ! output sample period in steps
integer     :: frstupdate   ! how often is restart file written
integer     :: feimode      ! extrapolation / interpolation mode
                            ! 1 - linear ramp I
integer     :: ftrjsample   ! how often save accumulator to "accumulator evolution"
logical     :: fapply_abf   ! on - apply ABF, off - do not apply ABF
logical     :: fprint_icf   ! T - print instantaneous collective forces (icf) and other data, F - do not print

logical     :: fenthalpy    ! collect data for enthalpy calculation
logical     :: fentropy     ! collect data for entropy calculation
real(PMFDP) :: fepotoffset
real(PMFDP) :: fekinoffset

! linear ramp mode II (feimode .eq. 1)
integer     :: fhramp_min   ! min value of linear ramp - ABF force is ignored below this value
integer     :: fhramp_max   ! max  value of linear ramp

! item list --------------------------------------------------------------------
type CVTypeABF
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data
    integer                 :: set              ! coordinate set index
    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins
end type CVTypeABF

! ----------------------
integer                     :: NumOfABFCVs      ! number of CVs in a group
type(CVTypeABF),allocatable :: ABFCVList(:)     ! definition of CVs

! accu types -------------------------------------------------------------------
type CVInfoTypeABF
    integer                  :: nbins           ! number of accumulator bins
    real(PMFDP)              :: min_value       ! left boundary of coordinate
    real(PMFDP)              :: max_value       ! left boundary of coordinate
    real(PMFDP)              :: bin_width       ! (right-left)/numbins
    real(PMFDP)              :: width           ! right - left
end type CVInfoTypeABF

! ----------------------

type ABFAccuType
    integer                       :: tot_cvs       ! total number of independent CVs
    type(CVInfoTypeABF), pointer  :: sizes(:)      ! CV information
    integer                       :: tot_nbins     ! number of total bins

    ! biasing force
    integer,pointer        :: nsamples(:)               ! number of hits into bins

    ! MICF
    real(PMFDP),pointer    :: micf(:,:)                 ! mean ICF - total
    real(PMFDP),pointer    :: m2icf(:,:)                ! M2 of ICF - total
    real(PMFDP),pointer    :: micf_pot(:,:)             ! mean ICF - potential energy part
    real(PMFDP),pointer    :: m2icf_pot(:,:)            ! M2 of ICF - potential energy part
    real(PMFDP),pointer    :: micf_kin(:,:)             ! mean ICF - kinetic energy part
    real(PMFDP),pointer    :: m2icf_kin(:,:)            ! M2 of ICF - kinetic energy part
    ! ENTHALPY
    real(PMFDP),pointer    :: metot(:)                  ! mean of total energy
    real(PMFDP),pointer    :: m2etot(:)                 ! M2 of total energy
    real(PMFDP),pointer    :: mepot(:)                  ! mean of potential energy
    real(PMFDP),pointer    :: m2epot(:)                 ! M2 of potential energy
    real(PMFDP),pointer    :: mekin(:)                  ! mean of kinetic energy
    real(PMFDP),pointer    :: m2ekin(:)                 ! M2 of kinetic energy
    ! ENTROPY
    real(PMFDP),pointer    :: cds_pp(:,:)               ! cds - potential/potential
    real(PMFDP),pointer    :: cds_pk(:,:)               ! cds - potential/kinetic
    real(PMFDP),pointer    :: cds_kp(:,:)               ! cds - kinetic/potential
    real(PMFDP),pointer    :: cds_kk(:,:)               ! cds - kinetic/kinetic
end type ABFAccuType

! ----------------------
type(ABFAccuType)           :: accumulator          ! accumulated forces
integer                     :: insidesamples
integer                     :: outsidesamples
! ----------------------

! global variables for force calculation ---------------------------------------
real(PMFDP)                 :: fdtx                 ! timestep

! global variables for abf - results -------------------------------------------
real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)

! helper arrays -------
real(PMFDP),allocatable     :: a0(:,:)        ! acceleration from previous step (t-dt)
real(PMFDP),allocatable     :: a1(:,:)        ! acceleration in current step (t)
real(PMFDP),allocatable     :: v0(:,:)        ! velocity in previous step

real(PMFDP),allocatable     :: la(:)          ! ABF force in coordinate direction
real(PMFDP),allocatable     :: zd0(:,:,:)     ! ZD0
real(PMFDP),allocatable     :: zd1(:,:,:)     ! ZD1
real(PMFDP),allocatable     :: pxi0(:)        !
real(PMFDP),allocatable     :: pxi1(:)        !
real(PMFDP),allocatable     :: pxip(:)        !
real(PMFDP),allocatable     :: pxim(:)        !
real(PMFDP),allocatable     :: pdum(:)        !
real(PMFDP),allocatable     :: avg_values(:)  ! average values of coordinates at t - 3/2dt

real(PMFDP),allocatable     :: cvaluehist0(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist1(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist2(:)   ! history of coordinate values
real(PMFDP),allocatable     :: cvaluehist3(:)   ! history of coordinate values

real(PMFDP)                 :: epothist0   ! history of Epot
real(PMFDP)                 :: epothist1   ! history of Epot
real(PMFDP)                 :: epothist2   ! history of Epot
real(PMFDP)                 :: epothist3   ! history of Epot

real(PMFDP)                 :: ekinhist0   ! history of Ekin
real(PMFDP)                 :: ekinhist1   ! history of Ekin
real(PMFDP)                 :: ekinhist2   ! history of Ekin
real(PMFDP)                 :: ekinhist3   ! history of Ekin

real(PMFDP), allocatable    :: icf_cache(:,:)   ! icf_cache(2*ncvs,fnstlim)

! ------------------------------------------------------------------------------

end module tabf_dat
