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
logical         :: freadranges      ! read ranges for CVs
logical         :: frestart         ! 1 - restart job with previous data, 0 - otherwise not

! output -----------------------------------------
integer         :: fsample          ! output sample period in steps
integer         :: fplevel          ! print level
integer         :: ftrjsample       ! how often save restart to "restart evolution"

! restart ----------------------------------------
integer         :: faccurst         ! number of steps for equilibration, it is ignored if job is restarted
integer         :: frstupdate       ! how often is restart file written


! constraints ------------------------------------
integer         :: fshakesolver     ! SHAKE solvers
                                    ! 0 - fixed SHAKE
                                    ! 1 - mixed SHAKE
                                    ! 2 - Newton-Raphson SHAKE
                                    ! 3 - diagonal SHAKE
                                    ! 4 - diagonal SHAKE with initial guess from the previous step
real(PMFDP)     :: flambdatol       ! tolerance for lambda optimization

integer         :: frattlesolver    ! RATTLE solvers
                                    ! 1 - matrix algebra RATTLE
real(PMFDP)     :: frveltol         ! residual for velocity in rattle/rattlev

integer         :: fmaxiter         ! maximum of iteration in lambda optimization
integer         :: flamsample       ! how often update lambda and metric tensor corrections

! enthalpy/entropy calculations
logical         :: fdhtds           ! collect data for enthalpy/entropy calculation
integer         :: fenesample      ! how often take samples

real(PMFDP)     :: fepotaverage
real(PMFDP)     :: fekinaverage


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
integer                    :: NumOfCONs         ! number of constraints
integer                    :: NumOfSHAKECONs    ! number of shake constraints in collision
integer                    :: NumOfAllCONs      ! number of constraints including shakes

type(CVTypeBM),allocatable :: CONList(:)        ! constraint list

! shake in collisions with CVs -------------------------------------------------
type CVTypeSHAKE
    integer                 :: at1
    integer                 :: at2
    real(PMFDP)             :: value
end type CVTypeSHAKE

type(CVTypeSHAKE),allocatable   :: SHAKECONList(:)           ! SHAKE definition of constraints

! serial/MPI variables ---------------------------------------------------------
integer                     :: NumOfCONAtoms            ! number of constrained atoms (unique list)
integer,allocatable         :: CONAtoms(:)              ! constrained atoms to test with SHAKE

! constants --------------------------------------------------------------------
integer, parameter  :: CON_SHAKESOL_FM      = 0     ! fixed shake: JAC(0,0)
integer, parameter  :: CON_SHAKESOL_MM      = 1     ! mixed shake: JAC(0,P)
integer, parameter  :: CON_SHAKESOL_NM      = 2     ! Newton-Raphson shake: JAC(P,P)
integer, parameter  :: CON_SHAKESOL_DI      = 3     ! diagonal JAC(0,P)
integer, parameter  :: CON_SHAKESOL_DIWG    = 4     ! diagonal JAC(0,P) with initial guess from the previous step

! global variables for lambda calculation --------------------------------------
real(PMFDP)                 :: isfdts           ! internal conversion factor
integer                     :: fsiter           ! number of iterations in shake solver
real(PMFDP),allocatable     :: lambdax(:)       ! list of Lagrange multipliers, internal units
real(PMFDP),allocatable     :: cv(:)            ! constraint value vector

real(PMFDP)                 :: nsupdates        ! number of shake updates
real(PMFDP)                 :: mfsiter          ! mean value of fsiter
real(PMFDP)                 :: m2fsiter         ! M2 moment of fsiter

! constants --------------------------------------------------------------------
integer, parameter  :: CON_RATTLESOL_MA     = 0     ! matrix algebra

! global variables for velocity update -----------------------------------------
real(PMFDP)                 :: isfdtr           ! internal conversion factor
integer                     :: friter           ! number of iterations in rattlev solver
real(PMFDP),allocatable     :: lambdav(:)       ! velocity lambdas - kappa, internal units

real(PMFDP)                 :: nrupdates        ! number of rattle updates
real(PMFDP)                 :: mfriter          ! mean value of friter
real(PMFDP)                 :: m2friter         ! M2 moment of friter

! metric tensor correction -----------------------------------------------------
real(PMFDP),allocatable     :: lambda(:)        ! total lambda with corrected units
real(PMFDP),allocatable     :: fwfac            ! current value of Fixman weight

! ICF
real(PMFDP),allocatable     :: CSTFrc(:,:)      ! forces after constraints are imposed
real(PMFDP),allocatable     :: icfp(:)          ! ICF - potential part
real(PMFDP),allocatable     :: icfk(:)          ! ICF - the other part
real(PMFDP),allocatable     :: icfk_vec(:,:)    ! helper array

! global variables for LU decomposition and other helper variable  -------------
real(PMFDP),allocatable     :: jac(:,:)         ! Jacobian matrix
real(PMFDP),allocatable     :: vv(:)            ! for LU decomposition
integer,allocatable         :: indx(:)
real(PMFDP),allocatable     :: zmata(:,:)       ! Z-matrix - all constraints
real(PMFDP),allocatable     :: zmats(:,:)       ! Z-matrix - SHAKE constraints

! history buffers ---------------------------------------------------------------
integer                     :: hist_len
integer                     :: hist_fidx

real(PMFDP),allocatable     :: lambdahist(:,:)
real(PMFDP),allocatable     :: epothist(:)
real(PMFDP),allocatable     :: ersthist(:)
real(PMFDP),allocatable     :: ekinhist(:)
real(PMFDP),allocatable     :: fwhist(:)
real(PMFDP),allocatable     :: icfphist(:,:)
real(PMFDP),allocatable     :: icfkhist(:,:)
logical,allocatable         :: enevalidhist(:)      ! is energy valid?

! ------------------------------------------------------------------------------
! ACCUMULATOR
! ------------------------------------------------------------------------------

type(PMFAccuType)           :: cstaccu
logical                     :: fallconstant     ! all CST CVs must be constant for PMFAccumulator
integer                     :: faccustep        ! number of integration steps for accumulator data sampling, TdS

real(PMFDP),allocatable     :: rbuf_B(:)        ! helper buffers
real(PMFDP),allocatable     :: rbuf_M(:,:)

! global variables for blue moon - results -------------------------------------
real(PMFDP)                 :: nsamples         ! total number of accumulated steps
real(PMFDP)                 :: mfw              ! mean of Fixman weights
real(PMFDP)                 :: m2fw             ! M2 of Fixman weights
real(PMFDP),allocatable     :: mlambda(:)       ! mean of lambdas
real(PMFDP),allocatable     :: m2lambda(:)      ! M2 of lambdas

! fdhtds  ----------------------------------------------------------------------
real(PMFDP)                 :: ntds             ! number of step for enthalpy and entropy calculations
real(PMFDP)                 :: fwsum            ! Fixman weights sum

real(PMFDP),allocatable     :: mlamtds(:)       ! mean of ICF - hamiltonian
real(PMFDP),allocatable     :: m2lamtds(:)      ! M2 of ICF - hamiltonian
real(PMFDP),allocatable     :: mlamtdsfw(:)     ! mean of ICF - hamiltonian  - Fixman weighted
real(PMFDP),allocatable     :: m2lamtdsfw(:)    ! M2 of ICF - hamiltonian

real(PMFDP)                 :: metot            ! mean of total energy
real(PMFDP)                 :: m2etot           ! M2 of total energy
real(PMFDP)                 :: meint            ! mean of internal energy
real(PMFDP)                 :: m2eint           ! M2 of internal energy
real(PMFDP)                 :: mepot            ! mean of potential energy
real(PMFDP)                 :: m2epot           ! M2 of potential energy
real(PMFDP)                 :: merst            ! mean of restraint energy
real(PMFDP)                 :: m2erst           ! M2 of restraint energy
real(PMFDP)                 :: mekin            ! mean of kinetic energy
real(PMFDP)                 :: m2ekin           ! M2 of kinetic energy

real(PMFDP)                 :: metotfw          ! mean of total energy - Fixman weighted
real(PMFDP)                 :: m2etotfw         ! M2 of total energy
real(PMFDP)                 :: meintfw          ! mean of internal energy
real(PMFDP)                 :: m2eintfw         ! M2 of internal energy
real(PMFDP)                 :: mepotfw          ! mean of potential energy
real(PMFDP)                 :: m2epotfw         ! M2 of potential energy
real(PMFDP)                 :: merstfw          ! mean of restraint energy
real(PMFDP)                 :: m2erstfw         ! M2 of restraint energy
real(PMFDP)                 :: mekinfw          ! mean of kinetic energy
real(PMFDP)                 :: m2ekinfw         ! M2 of kinetic energy

real(PMFDP),allocatable     :: micf(:)          ! mean of ICF
real(PMFDP),allocatable     :: m2icf(:)         ! M2 of ICF

real(PMFDP),allocatable     :: micffw(:)        ! mean of ICF - Fixman weighted
real(PMFDP),allocatable     :: m2icffw(:)       ! M2 of ICF

real(PMFDP),allocatable     :: micfpfw(:)       ! mean of ICF-P
real(PMFDP),allocatable     :: m2icfpfw(:)      ! M2 of ICF-P

real(PMFDP),allocatable     :: micfkfw(:)       ! mean of ICF-K
real(PMFDP),allocatable     :: m2icfkfw(:)      ! M2 of ICF-K

real(PMFDP),allocatable     :: c11ii(:)         ! co-variances covar(ICF,Eint)
real(PMFDP),allocatable     :: c11iifw(:)       ! co-variances covar(ICF,Eint)- Fixman weighted

real(PMFDP),allocatable     :: c11lt(:)         ! co-moments between lambda and total energy
real(PMFDP),allocatable     :: c11ltfw(:)       ! weighted co-moments between lambda and various energies
real(PMFDP),allocatable     :: c11lifw(:)
real(PMFDP),allocatable     :: c11lpfw(:)
real(PMFDP),allocatable     :: c11lrfw(:)
real(PMFDP),allocatable     :: c11lkfw(:)

!===============================================================================

end module cst_dat

