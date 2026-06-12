!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mta_dat

use pmf_sizes
use pmf_cvs
use pmf_accu

implicit none

! MASTER variables =============================================================

! control section --------------------------------------------------------------
integer     :: fmode        ! 0 - disable MTA, 1 - enabled MTA
integer     :: fsample      ! output sample period in steps
logical     :: frestart     ! restart accumulation from the previous restart file
integer     :: frstupdate   ! how often is restart file written

! item list --------------------------------------------------------------------
type CVTypeMTA
    integer                 :: cvindx           ! general description of coordinate
    class(CVType),pointer   :: cv               ! cv data

    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins
end type CVTypeMTA

! ----------------------

integer                     :: NumOfMTACVs          ! number of monitored CVs
type(CVTypeMTA),allocatable :: MTACVList(:)         ! monitored items

! ------------------------------------------------------------------------------

type,extends(PMFAccuType) :: MTAAccuType

    real(PMFDP),pointer    :: nsamples(:)               ! number of hits into bins

! MTA and iMTA=1.0/MTA
    real(PMFDP),pointer    :: mmta(:)                   ! mean MTA
    real(PMFDP),pointer    :: m2mta(:)                  ! M2 of MTA
    real(PMFDP),pointer    :: mimta(:)                  ! mean iMTA
    real(PMFDP),pointer    :: m2imta(:)                 ! M2 of iMTA
end type MTAAccuType

! ----------------------
type(MTAAccuType)           :: mtaaccu                  ! accumulated data
integer                     :: insidesamples
integer                     :: outsidesamples

! ----------------------

real(PMFDP),allocatable     :: fz(:,:)              ! Z matrix              in t
real(PMFDP),allocatable     :: fzinv(:,:)           ! inverse of Z matrix   in t
real(PMFDP)                 :: fzdet
real(PMFDP),allocatable     :: vv(:)                ! for LU decomposition
integer,allocatable         :: indx(:)              ! for LU decomposition

! ------------------------------------------------------------------------------

end module mta_dat

