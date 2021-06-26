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
    call cst_init_print_header
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
    fsample         = 500           ! output sample pariod in steps
    fplevel         = 0             ! print level

    frestart        = .false.       ! 1 - restart job with previous data, 0 - otherwise not
    faccurst        = 0             ! number of steps for equilibration, it is ignored if job is restarted
    frstupdate      = 500
    ftrjsample      = 0             ! how often save accumulator to "accumulator evolution"

    flambdasolver   = CON_LS_NM     ! Newton-Rapson solver
    flambdatol      = 1.0d-4        ! tolerance for lambda optimization
    fmaxiter        = 50            ! maximum of iteration in lambda optimization

    fenthalpy       = .false.       ! accumulate enthalpy
    fentropy        = .false.       ! accumulate entropy
    fepotoffset     = 0.0d0
    fekinoffset     = 0.0d0

    freadranges     = .true.        ! request full definitions of CVs

    NumOfCONs       = 0
    NumOfSHAKECONs  = 0
    NumOfConAtoms   = 0

end subroutine cst_init_dat

!===============================================================================
! Subroutine:  cst_init_print_header
!===============================================================================

subroutine cst_init_print_header

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
    write(PMF_OUT,130)  ' Constrained dynamics mode (fmode)    : ', fmode
    write(PMF_OUT,125)  ' Constraint definition file (fcstdef) : ', trim(fcstdef)
    write(PMF_OUT,130)  ' Total number of constraints          : ', NumOfCONs
    write(PMF_OUT,130)  ' SHAKE constraints in collisions      : ', NumOfSHAKECONs
    write(PMF_OUT,130)  ' Num of constrained atoms (no SHAKE)  : ', NumOfConAtoms
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Constraint optimization options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,140)  ' Lambda solver (flambdasolver)        : ', flambdasolver, &
                                                                 trim(cst_init_get_lsolver_name(flambdasolver))
    write(PMF_OUT,135)  ' Lambda tolerance (flambdatol)        : ', flambdatol
    write(PMF_OUT,130)  ' Maximum of iteration (fmaxiter)      : ', fmaxiter
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fcstout)                : ', trim(fcstout)
    write(PMF_OUT,130)  ' Sample period (fsample)              : ', fsample
    write(PMF_OUT,130)  ' Print level (fplevel)                : ', fplevel
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fcstrst)               : ', trim(fcstrst)
    write(PMF_OUT,125)  ' Restart from previous run (frestart) : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Accumulators reset (faccurst)        : ', faccurst
    write(PMF_OUT,120)

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

150 format(' == Constrained collective variable #',I4.4)

end subroutine cst_init_print_header

!===============================================================================
! Function:  cst_init_get_lsolver_name
!===============================================================================

character(80) function cst_init_get_lsolver_name(solver_id)

    use cst_dat

    implicit none
    integer     :: solver_id
    ! --------------------------------------------------------------------------

    select case(solver_id)
        case(CON_LS_NM)
            cst_init_get_lsolver_name = "Newton method"
        case(CON_LS_CM)
            cst_init_get_lsolver_name = "Chord method"
        case default
            cst_init_get_lsolver_name = "unknown"
    end select

    return

end function cst_init_get_lsolver_name

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
    integer        :: i,alloc_failed,ierr
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

! allocate arrays for LU decomposition ---------------------------------------
    allocate(vv(NumOfCONs), indx(NumOfCONs), stat= alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in LU decomposition!')
    end if

! allocate arrays for lambda calculation ---------------------------------------
    allocate(lambda(NumOfCONs), &
          cv(NumOfCONs), &
          jac(NumOfCONs,NumOfCONs), stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in lambda calculation!')
    end if

! allocate arrays for dF/dx calculation ---------------------------------------
    allocate( fz(NumOfCONs,NumOfCONs), &
           mlambda(NumOfCONs), &
           m2lambda(NumOfCONs), &
           stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in lambda calculation!')
    end if

    misrz   = 0.0d0
    m2isrz  = 0.0d0

! in velocity verlet there is an additional rattle step
    select case(fintalg)
        case(IA_VEL_VERLET)
            faccumulation = -1
            has_lambdav = .true.
        case default
            has_lambdav = .false.
            faccumulation = 0
    end select

    lambda(:) = 0.0d0
    mlambda(:) = 0.0d0
    m2lambda(:) = 0.0d0

    ! -----------------------------------------------
    if( has_lambdav ) then
        allocate( lambdav(NumOfCONs), &
                  matv(NumOfCONs,NumOfCONs), &
                  mlambdav(NumOfCONs), &
                  m2lambdav(NumOfCONs), &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used in velocity corrections!')
        end if

        lambdav(:) = 0.0d0
        matv(:,:) = 0.0d0
        mlambdav(:) = 0.0d0
        m2lambdav(:) = 0.0d0
    end if

! history ----------------------------------------
    allocate( lambda0(NumOfCONs),   &
              lambda1(NumOfCONs),   &
              stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used for history recording!')
    end if

    lambda0(:)  = 0.0d0
    lambda1(:)  = 0.0d0

! -----------------------------------------------
    if( fentropy ) then
        allocate( m11hp(NumOfCONs),     &
                  m11hk(NumOfCONs),     &
                  m11hr(NumOfCONs),     &
                  m22hp(NumOfCONs),     &
                  m22hk(NumOfCONs),     &
                  m22hr(NumOfCONs),     &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
        end if

        m11hp(:)    = 0.0d0
        m11hk(:)    = 0.0d0
        m11hr(:)    = 0.0d0
        m22hp(:)    = 0.0d0
        m22hk(:)    = 0.0d0
        m22hr(:)    = 0.0d0
    end if

    fentaccu    = 0.0d0
    metot       = 0.0d0
    m2etot      = 0.0d0
    mepot       = 0.0d0
    m2epot      = 0.0d0
    mekin       = 0.0d0
    m2ekin      = 0.0d0
    merst       = 0.0d0
    m2erst      = 0.0d0
    epothist0   = 0.0d0
    epothist1   = 0.0d0
    isrz0       = 0.0d0
    isrz1       = 0.0d0

    ! init accu
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

    allocate(   ibuf_B(cstaccu%tot_nbins),                    &
                rbuf_B(cstaccu%tot_nbins),                    &
                rbuf_M(cstaccu%tot_cvs,cstaccu%tot_nbins),    &
                stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CST accumulator!')
    endif

    return

end subroutine cst_init_core

!===============================================================================

end module cst_init
