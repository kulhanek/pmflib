!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module con_dat

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

! control section --------------------------------------------------------------
integer           :: fmode          ! 0 - disable BM, 1 - enabled BM
integer           :: fsample        ! output sample period in steps
integer           :: faccurst       ! number of steps for equilibration, it is ignored if job is restarted
integer           :: fplevel        ! print level
logical           :: frestart       ! 1 - restart job with previous data, 0 - otherwise not
integer           :: flambdasolver  ! 0 - Newton method, 1 - chord method
real(PMFDP)       :: flambdatol     ! tolerance for lambda optimization
integer           :: fmaxiter       ! maximum of iteration in lambda optimization

! item list --------------------------------------------------------------------
type CVTypeBM
    integer                 :: cvindx           ! index to PMF CV
    class(CVType),pointer   :: cv               ! cv data

    !mode - constant (C), incremental (I), change to value (V)
    character(PMF_MAX_MODE) :: mode

    real(PMFDP)             :: startvalue       ! start value
    real(PMFDP)             :: stopvalue        ! stop value
    real(PMFDP)             :: value            ! current value in time t

    real(PMFDP)             :: deviation        ! deviation between real and value
    real(PMFDP)             :: sdevtot          ! total sum of deviation squares
    logical                 :: value_set        ! initial value is user provided
end type CVTypeBM

! global variables for blue moon -----------------------------------------------
integer                    :: NumOfCONs         ! number of constraints including shakes
type(CVTypeBM),allocatable :: CONList(:)        ! constraint list

! shake in collisions with CVs -------------------------------------------------
type CVTypeSHAKE
    integer                 :: at1
    integer                 :: at2
    real(PMFDP)             :: value
end type CVTypeSHAKE

integer                         :: NumOfSHAKECONs  = 0       ! number of shake constraints in collision
type(CVTypeSHAKE),allocatable   :: SHAKECONList(:)           ! SHAKE definition of constraints

! serial/MPI variables ---------------------------------------------------------
integer             :: NumOfConAtoms            ! number of constrained atoms (unique list)
integer,allocatable :: ConAtoms(:)              ! constrained atoms with SHAKE

! constants --------------------------------------------------------------------
integer, parameter  :: CON_LS_NM         = 0    ! Newton method lambda solver
integer, parameter  :: CON_LS_CM         = 1    ! Chord method lambda solver

! global variables for lambda calculation -----------------------------------
integer                    :: fliter            ! number of iterations in lambda solver
real(PMFDP),allocatable    :: lambda(:)         ! list of Lagrange multipliers
real(PMFDP),allocatable    :: cv(:)             ! constraint value vector
real(PMFDP),allocatable    :: jac(:,:)          ! Jacobian matrix
logical                    :: has_lambdav       ! mu values (lambdav)

! global variables for LU decomposition  -----------------------------------
real(PMFDP),allocatable    :: vv(:)            ! for LU decomposition
integer,allocatable        :: indx(:)

! global variables for blue moon - results ---------------------------------
integer                    :: faccumulation    ! total number of accumulated steps
real(PMFDP),allocatable    :: fz(:,:)          ! Z matrix
real(PMFDP)                :: fzdet            ! current value of det(Z)
real(PMFDP)                :: isrztotal        ! accumulated value of inverse square root of fzdet
real(PMFDP)                :: isrztotals       ! accumulated value of square of inverse square root of fzdet
real(PMFDP),allocatable    :: lambdatotal(:)   ! accumulated value of lambda
real(PMFDP),allocatable    :: lambdatotals(:)  ! accumulated value of square of lambda

! global variables for velocity update ------------------------------
real(PMFDP),allocatable    :: matv(:,:)             ! left side matrix
real(PMFDP),allocatable    :: lambdav(:)            ! velocity lambdas
real(PMFDP),allocatable    :: lambdavtotal(:)       ! accumulated value of kappa
real(PMFDP),allocatable    :: lambdavtotals(:)      ! accumulated value of square of kappa

!===============================================================================

end module con_dat

