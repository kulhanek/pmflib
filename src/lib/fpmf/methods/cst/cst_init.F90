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

module cst_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_init_method
!===============================================================================

subroutine cst_init_method

    use cst_output
    use cst_restart
    use cst_constraints
    use cst_trajectory

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_init_all
    call cst_init_core

    ! we need first output then restart (writes info to output)
    call cst_output_open
    call cst_restart_read
    call cst_trajectory_open

    ! print header need value from restart file
    call cst_init_print_summary
    call cst_output_write_header

end subroutine cst_init_method

!===============================================================================
! Subroutine:  cst_init_dat
!===============================================================================

subroutine cst_init_dat

    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0             ! 0 - disable BM, 1 - enabled BM
    fsample         =  2500         ! output sample period in steps
    fplevel         = 0             ! print level

    frestart        = .false.       ! 1 - restart job with previous data, 0 - otherwise not
    faccurst        = 0             ! number of steps for equilibration, it is ignored if job is restarted
    frstupdate      = 10000
    ftrjsample      = 0             ! how often save accumulator to "accumulator evolution"

    flamsample      = 1
    fshakesolver    = CON_SHAKESOL_MM       ! mixed shake
    frattlesolver   = CON_RATTLESOL_MA      ! matrix algebra

    flambdatol      = 1.0d-7        ! tolerance for lambda optimization
    frveltol        = 1.0d-9        ! residual velocity in rattle/rattle-v
    fmaxiter        = 50            ! maximum of iteration in lambda optimization

    fenthalpy       = .false.       ! accumulate enthalpy
    fentropy        = .false.       ! accumulate entropy
    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0
    fenesample      = 1

    freadranges     = .false.        ! request full definitions of CVs

    NumOfCONs       = 0
    NumOfSHAKECONs  = 0
    NumOfConAtoms   = 0

    nsupdates       = 0.0d0
    mfsiter         = 0.0d0
    m2fsiter        = 0.0d0

    nrupdates       = 0.0d0
    mfriter         = 0.0d0
    m2friter        = 0.0d0

end subroutine cst_init_dat

!===============================================================================
! Subroutine:  cst_init_print_summary
!===============================================================================

subroutine cst_init_print_summary

    use prmfile
    use pmf_dat
    use cst_dat
    use cst_constraints

    implicit none
    integer :: i
    ! -----------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' -------------- FREE ENERGY CALCULATION BY CONSTRAINED DYNAMICS --------------- '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Cartesian Constraint Dynamics Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Constrained dynamics mode (fmode)       : ', fmode
    write(PMF_OUT,125)  ' Constraint definition file (fcstdef)    : ', trim(fcstdef)
    write(PMF_OUT,130)  ' Total number of constraints             : ', NumOfCONs
    write(PMF_OUT,130)  ' SHAKE constraints in collisions         : ', NumOfSHAKECONs
    write(PMF_OUT,130)  ' Num of constrained atoms (no SHAKE)     : ', NumOfConAtoms
    write(PMF_OUT,125)  ' Read CV ranges (freadranges)            : ', prmfile_onoff(freadranges)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Constraint optimization options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,140)  ' SHAKE solver (fshakesolver)             : ', fshakesolver, &
                                                                    trim(cst_init_get_shakesol_name(fshakesolver))
    write(PMF_OUT,135)  ' SHAKE lambda tolerance (flambdatol)     : ', flambdatol

    write(PMF_OUT,140)  ' RATTLE solver (frattlesolver)           : ', frattlesolver, &
                                                                    trim(cst_init_get_rattlesol_name(frattlesolver))
    write(PMF_OUT,135)  ' RATTLE velocity tolerance (frveltol)    : ', frveltol

    write(PMF_OUT,130)  ' Maximum of iteration (fmaxiter)         : ', fmaxiter

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Enthalpy options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Accumulate enthalpy (fenthalpy)         : ', prmfile_onoff(fenthalpy)
    write(PMF_OUT,125)  ' Accumulate enth. deriv. (fenthalpy_der) : ', prmfile_onoff(fenthalpy_der)
    write(PMF_OUT,145)  ' Potential energy offset (fepotaverage)  : ', pmf_unit_get_rvalue(EnergyUnit,fepotaverage),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,130)  ' Sampling for -TdS and dH (fenesample)   : ', fenesample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Entropy options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Accumulate entropy (fentropy)           : ', prmfile_onoff(fentropy)
    write(PMF_OUT,125)  ' Decompose entropy (fentdecomp)          : ', prmfile_onoff(fentdecomp)
    write(PMF_OUT,145)  ' Potential energy offset (fepotaverage)  : ', pmf_unit_get_rvalue(EnergyUnit,fepotaverage),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,145)  ' Kinetic energy offset (fekinaverage)    : ', pmf_unit_get_rvalue(EnergyUnit,fekinaverage), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,130)  ' Sampling for -TdS and dH (fenesample)   : ', fenesample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fcstrst)                  : ', trim(fcstrst)
    write(PMF_OUT,125)  ' Restart from previous run (frestart)    : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Accumulators reset (faccurst)           : ', faccurst
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fcstout)                   : ', trim(fcstout)
    write(PMF_OUT,130)  ' Sample period (fsample)                 : ', fsample
    write(PMF_OUT,130)  ' Print level (fplevel)                   : ', fplevel
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Trajectory file (fcsttrj)               : ', trim(fcsttrj)
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of constraints'
    write(PMF_OUT,120)  ' -------------------------------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfCONs-NumOfSHAKECONs
        write(PMF_OUT,150) i
        call cst_constraints_cst_info(CONList(i))
        write(PMF_OUT,120)
    end do

    if( NumOfSHAKECONs .gt. 0 ) then
        write(PMF_OUT,120)
        write(PMF_OUT,120)  ' List of SHAKE constraints in collision'
        write(PMF_OUT,120)  ' -------------------------------------------------------------------------------'
    end if

    write(PMF_OUT,120)
    do i=NumOfCONs-NumOfSHAKECONs+1,NumOfCONs
        write(PMF_OUT,150) i
        call cst_constraints_cst_info(CONList(i))
        write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)
135 format(A,E12.5)
140 format(A,I6,1X,A)
145 format(A,F10.1,1X,A)

150 format(' == Constrained collective variable #',I4.4)

end subroutine cst_init_print_summary

!===============================================================================
! Function:  cst_init_get_shakesol_name
!===============================================================================

character(80) function cst_init_get_shakesol_name(solver_id)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: solver_id
    ! --------------------------------------------------------------------------

    select case(solver_id)
        case(CON_SHAKESOL_FM)
            cst_init_get_shakesol_name = "Fixed SHAKE (LU)"
        case(CON_SHAKESOL_MM)
            cst_init_get_shakesol_name = "Mixed SHAKE (LU)"
        case(CON_SHAKESOL_NM)
            cst_init_get_shakesol_name = "Newton-Raphson SHAKE (LU)"
        case(CON_SHAKESOL_DI)
            cst_init_get_shakesol_name = "Mixed SHAKE (diagonal solver)"
        case(CON_SHAKESOL_DIWG)
            cst_init_get_shakesol_name = "Mixed SHAKE (diagonal solver) with initial guess"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented shake solver in cst_init_get_shakesol_name!')
    end select

    return

end function cst_init_get_shakesol_name

!===============================================================================
! Function:  cst_init_get_rattlesol_name
!===============================================================================

character(80) function cst_init_get_rattlesol_name(solver_id)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: solver_id
    ! --------------------------------------------------------------------------

    select case(solver_id)
        case(CON_SHAKESOL_FM)
            cst_init_get_rattlesol_name = "Matrix Algebra RATTLE"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented rattle solver in cst_init_get_rattlesol_name!')
    end select

    return

end function cst_init_get_rattlesol_name

!===============================================================================
! Subroutine:  cst_init_add_shake_csts
!===============================================================================

subroutine cst_init_add_shake_csts

    use pmf_utils
    use pmf_dat
    use cv_ds
    use cst_dat
    use cst_constraints
    use pmf_unit

    implicit none
    type(CVPointer),allocatable    :: CVList_backup(:)
    type(CVTypeBM),allocatable     :: CONList_backup(:)
    integer                        :: i,cvid,conid,alloc_failed
    ! ------------------------------------------------------------------------------

    if( NumOfSHAKECONs .eq. 0 .or. NumOfCONs .eq. 0 ) return

    ! backup old CVs
    allocate(CVList_backup(NumOfCVs),   &
          CONList_backup(NumOfCONs), &
          stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Unable to allocate memory for CVList_backup/CONList_backup array!')
    end if

    CVList_backup(:) = CVList(:)
    CONList_backup(:) = CONList(:)

    ! reallocate
    if( allocated(CVList) ) deallocate(CVList)
    if( allocated(CONList) ) deallocate(CONList)

    allocate(CVList(NumOfCVs + NumOfSHAKECONs),  &
          CONList(NumOfCONs + NumOfSHAKECONs), &
          stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CVList/CONList array!')
    end if

    do i=1,NumOfCVs
        CVList(i)  = CVList_backup(i)
    end do
    do i=1,NumOfCONs
        CONList(i) = CONList_backup(i)
    end do

    deallocate(CVList_backup)
    deallocate(CONList_backup)

    ! add SHAKE constraints
    do i= 1,NumOfSHAKECONs
        cvid  = NumOfCVs + i
        conid = NumOfCONs + i
        ! CV -----------------------------------------
        allocate(CVTypeDS::CVList(cvid)%cv, &
                 stat = alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for SHAKE constraint!')
        end if
        call CVList(cvid)%cv%reset_cv()
        CVList(cvid)%cv%ctype     = 'DS'
        CVList(cvid)%cv%unit      = pmf_unit_power_unit(LengthUnit,2)
        CVList(cvid)%cv%idx       = i
        CVList(cvid)%cv%name      = 'SHAKE'
        CVList(cvid)%cv%natoms    = 2
        CVList(cvid)%cv%ngrps     = 2
        allocate(CVList(cvid)%cv%grps(CVList(cvid)%cv%ngrps), &
                 CVList(cvid)%cv%rindexes(2), &
                 CVList(cvid)%cv%lindexes(2), &
                 stat = alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CVList(i)%grps array!')
        end if
        CVList(cvid)%cv%grps(1)   = 1
        CVList(cvid)%cv%grps(2)   = 2
        CVList(cvid)%cv%rindexes(1) = SHAKECONList(i)%at1
        CVList(cvid)%cv%rindexes(2) = SHAKECONList(i)%at2
        ! CST ----------------------------------------
        call cst_constraints_reset_con(CONList(conid))
        CONList(conid)%cvindx       = cvid
        CONList(conid)%cv           => CVList(cvid)%cv
        CONList(conid)%mode         = 'C'
        CONList(conid)%value        = SHAKECONList(i)%value
        CONList(conid)%value_set    = .true.
    end do

    ! correct numbers
    NumOfCVs = NumOfCVs + NumOfSHAKECONs
    NumOfCONs = NumOfCONs + NumOfSHAKECONs

end subroutine cst_init_add_shake_csts

!===============================================================================
! Subroutine:  cst_init_cst_atoms
!===============================================================================

subroutine cst_init_cst_atoms

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                 :: ci, na, i, j, k, alloc_failed
    logical                 :: found
    integer,allocatable     :: tmp_indexes(:)
    ! ------------------------------------------------------------------------------

    ! count involved atoms
    na = 0

    do i=1, NumOfCONs
        ci = CONList(i)%cvindx
        na = na + CVList(ci)%cv%natoms
    end do

    ! allocate index array
    allocate(tmp_indexes(na),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for tmp_indexes array!')
    end if

    tmp_indexes(:) = 0
    NumOfConAtoms = 0;

    ! add conflicting atoms
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do j=1,CVList(ci)%cv%natoms
            found = .false.
            do k=1,na
                if( tmp_indexes(k) .eq. CVList(ci)%cv%rindexes(j) ) then
                    found = .true.
                    exit
                end if
            end do
            if( .not. found ) then
                NumOfConAtoms = NumOfConAtoms + 1
                tmp_indexes(NumOfConAtoms) = CVList(ci)%cv%rindexes(j)
            end if
        end do
    end do

    if( NumOfConAtoms .eq. 0 ) then
        ! release temporary array
        deallocate(tmp_indexes)
        return
    end if

    ! create final array
    allocate(ConAtoms(NumOfConAtoms),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for ConAtoms array!')
    end if

    do i=1,NumOfConAtoms
        ConAtoms(i) = tmp_indexes(i)
    end do

    ! release temporary array
    deallocate(tmp_indexes)

end subroutine cst_init_cst_atoms

#ifdef MPI

!===============================================================================
! Subroutine:  cst_init_mpi_bcast_constraints
! this is required for cst_shake_checkatom
!===============================================================================

subroutine cst_init_mpi_bcast_constraints

    use mpi
    use cst_dat
    use pmf_utils
    use pmf_dat

    implicit none
    integer        :: alloc_failed,ierr
    ! -----------------------------------------------------------------------------

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> Broadcasting constrained atoms (only master is reporting)'
    end if

    ierr = MPI_SUCCESS

    ! integers --------------------------------------
    call mpi_bcast(NumOfConAtoms, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to broadcast the value of NumOfConAtoms!')
    end if

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A,I6)') '   Number of constrained atoms: ', NumOfConAtoms
        write(PMF_DEBUG+fmytaskid,*)
    end if

    if( NumOfConAtoms .eq. 0 ) return ! no atoms are constrained

    ! allocate arrays on slaves ---------------------
    if( .not. fmaster ) then
        allocate(ConAtoms(NumOfConAtoms),    &
                stat= alloc_failed )
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate ConAtoms!')
        end if
    end if

    ! transfer ConAtoms
    call mpi_bcast(ConAtoms, NumOfConAtoms, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to broadcast the ConAtoms array!')
    end if

end subroutine cst_init_mpi_bcast_constraints

#endif

!===============================================================================
! Subroutine:  cst_init_core
!===============================================================================

subroutine cst_init_core

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer      :: alloc_failed, i, tot_nbins
    ! ------------------------------------------------------------------------------

! setup conversion factors
    select case(fintalg)
        case(IA_LEAP_FROG)
            isfdts = 1.0d0/(fdt*fdt) * PMF_L2CL
            isfdtr = 0.0d0
        case(IA_LF_MIDDLE) ! FIXME
            isfdts = 1.0d0/(fdt*fdt) * PMF_L2CL
            isfdtr = 1.0d0/(fdt*PMF_VDT2DT) * PMF_L2CL
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_init_print_summary!')
    end select

! required always - det(Z) is calculate in core_analyse
! allocate arrays for LU decomposition
    allocate(vv(NumOfCONs),             &
             indx(NumOfCONs),           &
             jac(NumOfCONs,NumOfCONs),  &
             zmata(NumOfCONs,NumOfCONs), stat= alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in LU decomposition!')
    end if

    if( NumOfSHAKECONS .gt. 0 ) then
        allocate( zmats(NumOfSHAKECONS,NumOfSHAKECONS), stat= alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used in LU decomposition!')
        end if
    end if

! allocate arrays for lambda calculation
    select case(fintalg)
        case(IA_LEAP_FROG) ! FIXME
            allocate(lambdax(NumOfCONs), cv(NumOfCONs), stat= alloc_failed )
            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used in lambda calculation!')
            end if
            lambdax(:) = 0.0d0
            cv(:) = 0.0d0
        case(IA_VEL_VERLET)
            allocate(lambdax(NumOfCONs), lambdav(NumOfCONs),  &
                     cv(NumOfCONs), stat= alloc_failed )
            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used in lambda calculation!')
            end if
            lambdax(:) = 0.0d0
            lambdav(:) = 0.0d0
            cv(:) = 0.0d0
        case(IA_LF_MIDDLE)
            allocate(lambdax(NumOfCONs), lambdav(NumOfCONs),  &
                     cv(NumOfCONs), stat= alloc_failed )
            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used in lambda calculation!')
            end if
            lambdax(:) = 0.0d0
            lambdav(:) = 0.0d0
            cv(:) = 0.0d0
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_init_print_summary!')
    end select


! history buffers
    hist_len = 2
    hist_fidx = -1

    allocate( lambdahist(NumOfCONs,hist_len),   &
              epothist(hist_len),               &
              ersthist(hist_len),               &
              ekinhist(hist_len),               &
              ifwhist(hist_len),                &
              icfphist(NumOfCONs,hist_len),     &
              icfkhist(NumOfCONs,hist_len),     &
              enevalidhist(hist_len),           &
              stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used for history recording!')
    end if

    lambdahist(:,:) = 0.0d0
    epothist(:)     = 0.0d0
    ersthist(:)     = 0.0d0
    ekinhist(:)     = 0.0d0
    ifwhist(:)      = 0.0d0
    icfphist(:,:)   = 0.0d0
    icfkhist(:,:)   = 0.0d0
    enevalidhist(:) = .false.

! accumulator setup for free energy calculation
    allocate( lambda(NumOfCONs),    &
              mlambda(NumOfCONs),   &
              m2lambda(NumOfCONs),  &
              stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in lambda calculation!')
    end if

    nsamples    = 0.0d0
    misrz       = 0.0d0
    m2isrz      = 0.0d0
    lambda(:)   = 0.0d0
    mlambda(:)  = 0.0d0
    m2lambda(:) = 0.0d0

! accumulator setup for entropy and enthalpy
    fene_step = 0

    if( fenthalpy .or. fentropy ) then
        ntds    = 0.0d0
        mfw     = 0.0d0
        m2fw    = 0.0d0
    end if

    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
        meint       = 0.0d0
        m2eint      = 0.0d0
        mepot       = 0.0d0
        m2epot      = 0.0d0
        merst       = 0.0d0
        m2erst      = 0.0d0
        mekin       = 0.0d0
        m2ekin      = 0.0d0

        meintfw     = 0.0d0
        m2eintfw    = 0.0d0
        mepotfw     = 0.0d0
        m2epotfw    = 0.0d0
        merstfw     = 0.0d0
        m2erstfw    = 0.0d0
        mekinfw     = 0.0d0
        m2ekinfw    = 0.0d0
    end if

    if( fenthalpy .and. fenthalpy_der ) then
        allocate( CSTFrc(3,NumOfLAtoms),    &
                  icfp(NumOfCONs),          &
                  icfk(NumOfCONs),          &
                  icfk_vec(3,NumOfLAtoms),  &
                  micfp(NumOfCONs),         &
                  m2icfp(NumOfCONs),        &
                  micfk(NumOfCONs),         &
                  m2icfk(NumOfCONs),        &
                  micf(NumOfCONs),         &
                  m2icf(NumOfCONs),        &
                  c11pp(NumOfCONs),         &
                  micfpfw(NumOfCONs),       &
                  m2icfpfw(NumOfCONs),      &
                  micfkfw(NumOfCONs),       &
                  m2icfkfw(NumOfCONs),      &
                  micfpeintfw(NumOfCONs),       &
                  m2icfpeintfw(NumOfCONs),      &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
        end if
        CSTFrc(:,:) = 0.0d0
        icfp(:)     = 0.0d0
        icfk(:)     = 0.0d0
        icfk_vec(:,:) = 0.0d0
        micfp(:)    = 0.0d0
        m2icfp(:)   = 0.0d0
        micfk(:)    = 0.0d0
        m2icfk(:)   = 0.0d0
        micf(:)     = 0.0d0
        m2icf(:)    = 0.0d0
        c11pp(:)    = 0.0d0

        micfpfw(:)  = 0.0d0
        m2icfpfw(:) = 0.0d0
        micfkfw(:)  = 0.0d0
        m2icfkfw(:) = 0.0d0
        micfpeintfw(:)  = 0.0d0
        m2icfpeintfw(:) = 0.0d0
    end if

    if( fentropy ) then
        allocate( mpp(NumOfCONs),       &
                  m2pp(NumOfCONs),      &
                  mpn(NumOfCONs),       &
                  m2pn(NumOfCONs),      &
                  mhicf(NumOfCONs),     &
                  m2hicf(NumOfCONs),    &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
        end if

        metot       = 0.0d0
        m2etot      = 0.0d0
        mpp(:)      = 0.0d0
        m2pp(:)     = 0.0d0
        mpn(:)      = 0.0d0
        m2pn(:)     = 0.0d0
        mhicf(:)    = 0.0d0
        m2hicf(:)   = 0.0d0
    end if

    if( fentropy .and. fentdecomp ) then
        allocate( c11hp(NumOfCONs),     &
                  c11hr(NumOfCONs),     &
                  c11hk(NumOfCONs),     &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
        end if
        c11hp(:)    = 0.0d0
        c11hr(:)    = 0.0d0
        c11hk(:)    = 0.0d0
    end if

! init PMF accu
    cstaccu%tot_cvs = NumOfCONs - NumOfSHAKECONs

    allocate(cstaccu%sizes(cstaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CST accumulator!')
    endif

    tot_nbins       = 1
    fallconstant    = .true.

    do i=1,cstaccu%tot_cvs
        if( CONList(i)%mode .ne. 'C' ) then
            fallconstant = .false.
        end if
        if( CONList(i)%ibin .eq. 0 ) then
            freadranges = .false.   ! without all cv bins this cannot be enabled
        end if
    end do

    do i=1,cstaccu%tot_cvs
        if( freadranges ) then
            cstaccu%sizes(i)%min_value  = CONList(i)%min_value
            cstaccu%sizes(i)%max_value  = CONList(i)%max_value
            cstaccu%sizes(i)%nbins      = CONList(i)%nbins
            cstaccu%sizes(i)%width      = abs(cstaccu%sizes(i)%max_value - cstaccu%sizes(i)%min_value)
            cstaccu%sizes(i)%bin_width  = cstaccu%sizes(i)%width / cstaccu%sizes(i)%nbins
        else
            cstaccu%sizes(i)%min_value  = CONList(i)%value
            cstaccu%sizes(i)%max_value  = CONList(i)%value
            cstaccu%sizes(i)%nbins      = 1
            cstaccu%sizes(i)%width      = 0.0d0
            cstaccu%sizes(i)%bin_width  = 0.0d0
        end if
        cstaccu%sizes(i)%cv => CONList(i)%cv
        tot_nbins = tot_nbins * cstaccu%sizes(i)%nbins
    end do

    cstaccu%tot_nbins = tot_nbins

    allocate(   rbuf_B(cstaccu%tot_nbins),                    &
                rbuf_M(cstaccu%tot_cvs,cstaccu%tot_nbins),    &
                stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CST accumulator!')
    endif

    return

end subroutine cst_init_core

!===============================================================================

end module cst_init
