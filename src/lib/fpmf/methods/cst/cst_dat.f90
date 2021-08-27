!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module cst_dat

use pmf_sizes
use pmf_constants
use pmf_dat
use pmf_accu

implicit none

! control section --------------------------------------------------------------
integer         :: fmode            ! 0 - disable BM, 1 - enabled BM
integer         :: fsample          ! output sample period in steps
integer         :: faccurst         ! number of steps for equilibration, it is ignored if job is restarted
integer         :: fplevel          ! print level
logical         :: frestart         ! 1 - restart job with previous data, 0 - otherwise not
integer         :: frstupdate       ! how often is restart file written
integer         :: ftrjsample       ! how often save restart to "restart evolution"
integer         :: flambdasolver    ! 0 - Newton method, 1 - chord method, 2 - Newton + SVD
real(PMFDP)     :: flambdatol       ! tolerance for lambda optimization
integer         :: fmaxiter         ! maximum of iteration in lambda optimization
integer         :: fsamplefreq      ! how often take samples
logical         :: freadranges      ! read ranges for CVs

logical         :: fenthalpy        ! collect data for enthalpy calculation
logical         :: fentropy         ! collect data for entropy calculation
real(PMFDP)     :: fepotoffset
real(PMFDP)     :: fekinoffset

real(PMFDP)     :: frcond

! item list --------------------------------------------------------------------
type CVTypeBM
    integer                 :: cvindx           ! index to PMF CV
    class(CVType),pointer   :: cv               ! cv data

    character(PMF_MAX_MODE) :: mode             ! mode - constant (C)
                                                !      - incremental (I)
                                                !      - change to value (V)
                                                !      - controlled steering (S)

    real(PMFDP)             :: startvalue       ! start value
    real(PMFDP)             :: stopvalue        ! stop value
    real(PMFDP)             :: value            ! current value in time t
    integer                 :: ibin             ! bin corresponding to value

    real(PMFDP)             :: deviation        ! deviation between real and value
    real(PMFDP)             :: sdevtot          ! total sum of deviation squares
    logical                 :: value_set        ! initial value user provided
    real(PMFDP),pointer     :: control_values(:) ! values for controlled steering

    real(PMFDP)             :: min_value        ! left range
    real(PMFDP)             :: max_value        ! right range
    integer                 :: nbins            ! number of bins

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

! ------------------------------------------------------------------------------

type(PMFAccuType)   :: cstaccu
logical             :: fallconstant

integer,allocatable         :: ibuf_B(:)
real(PMFDP),allocatable     :: rbuf_B(:)
real(PMFDP),allocatable     :: rbuf_M(:,:)

! constants --------------------------------------------------------------------
integer, parameter  :: CON_LS_NM         = 0    ! Newton method lambda solver - LU
integer, parameter  :: CON_LS_CM         = 1    ! Chord method lambda solver
integer, parameter  :: CON_LS_NM_SVD     = 2    ! Newton method lambda solver - SVD

! global variables for lambda calculation -----------------------------------
integer                     :: fliter           ! number of iterations in lambda solver
real(PMFDP),allocatable     :: lambda(:)        ! list of Lagrange multipliers
real(PMFDP),allocatable     :: cv(:)            ! constraint value vector
real(PMFDP),allocatable     :: jac(:,:)         ! Jacobian matrix
logical                     :: has_lambdav      ! mu values (lambdav)

! global variables for LU decomposition  -----------------------------------
real(PMFDP),allocatable     :: vv(:)            ! for LU decomposition
integer                     :: lwork            ! for SVD decomposition
real(PMFDP),allocatable     :: work(:)          ! for SVD decomposition
integer,allocatable         :: indx(:)

! global variables for blue moon - results ---------------------------------
integer                     :: faccumulation    ! total number of accumulated steps
real(PMFDP),allocatable     :: fz(:,:)          ! Z matrix
real(PMFDP)                 :: fzdet            ! current value of det(Z)
real(PMFDP)                 :: misrz            ! mean of inverse square root of fzdet
real(PMFDP)                 :: m2isrz           ! M2 of inverse square root of fzdet
real(PMFDP),allocatable     :: mlambda(:)       ! mean of lambdas
real(PMFDP),allocatable     :: m2lambda(:)      ! M2 of lambdas

! global variables for velocity update ------------------------------
real(PMFDP),allocatable     :: matv(:,:)        ! left side matrix
real(PMFDP),allocatable     :: lambdav(:)       ! velocity lambdas - kappa
real(PMFDP),allocatable     :: mlambdav(:)      ! mean of kappa
real(PMFDP),allocatable     :: m2lambdav(:)     ! M2 of kappa

! enthalpy and entropy ----------------------------------------------
integer                     :: fentaccu         ! number of step for enthalpy and entropy calculations
! all at t+dt
real(PMFDP)                 :: metot            ! mean of total energy
real(PMFDP)                 :: m2etot           ! M2 of total energy
real(PMFDP)                 :: mepot            ! mean of potential energy
real(PMFDP)                 :: m2epot           ! M2 of potential energy
real(PMFDP)                 :: mekin            ! mean of kinetic energy
real(PMFDP)                 :: m2ekin           ! M2 of kinetic energy
real(PMFDP)                 :: merst            ! mean of restraint energy
real(PMFDP)                 :: m2erst           ! M2 of restraint energy

! data for entropy, co-moments between dH/dx a H (total energy), Pot, Kin, Rst energies
real(PMFDP),allocatable     :: c11hh(:)
real(PMFDP),allocatable     :: c11hp(:)
real(PMFDP),allocatable     :: c11hk(:)
real(PMFDP),allocatable     :: c11hr(:)

real(PMFDP),allocatable     :: lambda0(:)       ! list of Lagrange multipliers, t-dt
real(PMFDP),allocatable     :: lambda1(:)       ! list of Lagrange multipliers, t
real(PMFDP)                 :: epothist0        ! history of Epot, t-dt
real(PMFDP)                 :: epothist1        ! history of Epot, t
real(PMFDP)                 :: ersthist0        ! history of Erst, t-dt
real(PMFDP)                 :: ersthist1        ! history of Erst, t
real(PMFDP)                 :: isrz0            ! history of isrz, t-dt
real(PMFDP)                 :: isrz1            ! history of isrz, t

! call in cst_core_main_lf
! lambda(t), epot(t), ekin(t-dt)

!===============================================================================

end module cst_dat

