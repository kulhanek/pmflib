!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025-2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    freadranges     = .false.       ! request full definitions of CVs

    fsample         =  2500         ! output sample period in steps
    fplevel         = 0             ! print level

    frestart        = .false.       ! 1 - restart job with previous data, 0 - otherwise not
    faccurst        = 0             ! number of steps for equilibration, it is ignored if job is restarted
    frstupdate      = 10000
    ftrjsample      = 0             ! how often save accumulator to "accumulator evolution"

    flam_sample     = 1

    fmdconmode      = CON_MDCON_INCLUDE
    fmdcon_cvtype   = CON_CVTYPE_DS

    fshakesolver    = CON_SHAKESOL_MM       ! mixed shake
    frattlesolver   = CON_RATTLESOL_MA      ! matrix algebra

    fshake_fdamp    = 0.0d0                 ! FIXME
    frattle_fdamp   = 0.0d0                 ! FIXME
    flamsol_fdamp   = 0.0d0                 ! FIXME
    ficf_fdamp      = 0.0d0                 ! FIXME

    flambdatol      = 1.0d-7        ! tolerance for lambda optimization
    frveltol        = 1.0d-9        ! residual velocity in rattle/rattle-v
    fmaxiter        = 50            ! maximum of iteration in lambda optimization
    frcond          = 1e-7

    fintcalc        = .false.       ! accumulate internal energy
    fint_der        = .false.
    ftds_icfsol     = CON_ICFSOL_V1

    ftdscalc        = .false.
    ftds_decomp     = .false.
    ftds_lamsol     = CON_LAMSOL_MD
    ftds_ekinsrc    = CON_EKINSRC_V4

    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0
    ftds_sample     = 1

    NumOfCONs       = 0
    NumOfMDCONs     = 0
    NumOfAllCONs    = 0
    NumOfConAtoms   = 0
    NumOfExcMDCONs  = 0

    nsupdates       = 0.0d0
    mfsiter         = 0.0d0
    m2fsiter        = 0.0d0

    nrupdates       = 0.0d0
    mfriter         = 0.0d0
    m2friter        = 0.0d0

    fpmf_div_dh     = 1e-5

    fpmf_sdiv_dh    = 1e-5
    fpmf_sdiv_S     = 16
    fpmf_sdiv_qr    = .false.

    frmmdcon_zdet   = .false.

    fdump_data      = .false.

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
    write(PMF_OUT,120)  ' Cartesian Constrained Dynamics Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Constrained dynamics mode (fmode)       : ', fmode
    write(PMF_OUT,125)  ' Constraint definition file (fcstdef)    : ', trim(fcstdef)
    write(PMF_OUT,125)  ' Read CV ranges (freadranges)            : ', prmfile_onoff(freadranges)

    write(PMF_OUT,130)  ' Number of constraints                   : ', NumOfCONs
    write(PMF_OUT,130)  ' Num of constrained atoms (no MD cons)   : ', NumOfConAtoms
    write(PMF_OUT,130)  ' MD constraints in collisions            : ', NumOfMDCONs
    write(PMF_OUT,130)  ' Total number of constraints             : ', NumOfAllCONs
    write(PMF_OUT,130)  ' Excluded MD constraints                 : ', NumOfExcMDCONs
    write(PMF_OUT,125)  ' Remove MD con. from FW (frmmdcon_zdet)  : ', prmfile_onoff(frmmdcon_zdet)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Constraint optimization options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,140)  ' SHAKE solver (fshakesolver)             : ', fshakesolver, &
                                                                       trim(cst_init_get_shakesol_name(fshakesolver))
    write(PMF_OUT,135)  ' SHAKE diag. reg. (fshake_fdamp)         : ', fshake_fdamp
    write(PMF_OUT,135)  ' SHAKE lambda tolerance (flambdatol)     : ', flambdatol

    write(PMF_OUT,140)  ' RATTLE solver (frattlesolver)           : ', frattlesolver, &
                                                                       trim(cst_init_get_rattlesol_name(frattlesolver))
    write(PMF_OUT,135)  ' RATTLE diag. reg. (frattle_fdamp)       : ', frattle_fdamp
    write(PMF_OUT,135)  ' RATTLE velocity tolerance (frveltol)    : ', frveltol

    write(PMF_OUT,130)  ' Maximum of iteration (fmaxiter)         : ', fmaxiter

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' MD constraints:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,140)  ' MD constraints mode (fmdconmode)        : ', fmdconmode, &
                                                                       trim(cst_init_get_mdcon_mode_name(fmdconmode))
    write(PMF_OUT,140)  ' CV type for SHAKEn bonds (fmdcon_cvtype): ', fmdcon_cvtype, &
                                                                       trim(cst_init_get_mdcon_cvtype_name(fmdcon_cvtype))
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Internal energy/Entropy options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Accumulate internal energy (fintcalc)   : ', prmfile_onoff(fintcalc)
    write(PMF_OUT,125)  ' Accumulate intene deriv. (fint_der)     : ', prmfile_onoff(fint_der)
    write(PMF_OUT,140)  ' ICF solver (ftds_icfsol)                : ', ftds_icfsol, &
                                                                       trim(cst_init_get_icfsol_name(ftds_icfsol))
    write(PMF_OUT,135)  ' ICF diag. reg. (ficf_fdamp)             : ', ficf_fdamp
    write(PMF_OUT,130)  ' ICF sdiv num of trials (fpmf_sdiv_S)    : ', fpmf_sdiv_S
    write(PMF_OUT,125)  ' ICF sdiv QR orthog. (fpmf_sdiv_qr)      : ', prmfile_onoff(fpmf_sdiv_qr)
    write(PMF_OUT,135)  ' ICF sdiv dh factor (fpmf_sdiv_dh)       : ', fpmf_sdiv_dh

    write(PMF_OUT,125)  ' Accumulate entropy (ftdscalc)           : ', prmfile_onoff(ftdscalc)
    write(PMF_OUT,125)  ' Decompose entropy (ftds_decomp)         : ', prmfile_onoff(ftds_decomp)

    write(PMF_OUT,145)  ' Potential energy offset (fepotaverage)  : ', pmf_unit_get_rvalue(EnergyUnit,fepotaverage),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,145)  ' Kinetic energy offset (fekinaverage)    : ', pmf_unit_get_rvalue(EnergyUnit,fekinaverage), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'

    write(PMF_OUT,140)  ' Lambda solver (ftds_lamsol)             : ', ftds_lamsol, &
                                                                       trim(cst_init_get_lamsol_name(ftds_lamsol))
    write(PMF_OUT,135)  ' LAMSOL diag. reg. (flamsol_fdamp)       : ', flamsol_fdamp
    write(PMF_OUT,140)  ' Kinetic energy source (ftds_ekinsrc)    : ', ftds_ekinsrc, &
                                                                       trim(cst_init_get_ekinsrc_name(ftds_ekinsrc))

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fcstrst)                  : ', trim(fcstrst)
    write(PMF_OUT,125)  ' Restart from previous run (frestart)    : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Accumulators reset (faccurst)           : ', faccurst
    write(PMF_OUT,130)  ' Sampling for FEN (flam_sample)          : ', flam_sample
    write(PMF_OUT,130)  ' Sampling for TDS and INT (ftds_sample)  : ', ftds_sample
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fcstout)                   : ', trim(fcstout)
    write(PMF_OUT,130)  ' Sample period (fsample)                 : ', fsample
    write(PMF_OUT,130)  ' Print level (fplevel)                   : ', fplevel
    write(PMF_OUT,125)  ' Dump dU/mTdS data (fdump_data)          : ', prmfile_onoff(fdump_data)
    write(PMF_OUT,125)  ' Dump file (fcstdump)                    : ', trim(fcstdump)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Trajectory file (fcsttrj)               : ', trim(fcsttrj)
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of constraints'
    write(PMF_OUT,120)  ' -------------------------------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfAllCONs-NumOfMDCONs
        write(PMF_OUT,150) i
        call cst_constraints_cst_info(CONList(i),.true.)
        write(PMF_OUT,120)
    end do

    if( NumOfMDCONs .gt. 0 ) then
        write(PMF_OUT,120)
        write(PMF_OUT,120)  ' List of MD constraints in collision'
        write(PMF_OUT,120)  ' -------------------------------------------------------------------------------'
    end if

    write(PMF_OUT,120)
    do i=NumOfAllCONs-NumOfMDCONs+1,NumOfAllCONs
        write(PMF_OUT,150) i
        call cst_constraints_cst_info(CONList(i),.false.)
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
            cst_init_get_shakesol_name = "Newton-Raphson SHAKE (LL)"
        case(CON_SHAKESOL_DI)
            cst_init_get_shakesol_name = "Mixed SHAKE (diagonal solver)"
        case(CON_SHAKESOL_DIWG)
            cst_init_get_shakesol_name = "Mixed SHAKE (diagonal solver) with initial guess"
        case(CON_SHAKESOL_NMSVD)
            cst_init_get_shakesol_name = "Newton-Raphson SHAKE (SVD)"
        case(CON_SHAKESOL_NMLU)
            cst_init_get_shakesol_name = "Newton-Raphson SHAKE (LU)"
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
! Function:  cst_init_get_mdcon_mode_name
!===============================================================================

character(80) function cst_init_get_mdcon_mode_name(mode)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(CON_MDCON_EXCLUDE)
            cst_init_get_mdcon_mode_name = "exclude"
        case(CON_MDCON_INCLUDE)
            cst_init_get_mdcon_mode_name = "include"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented CV type in cst_init_get_mdcon_mode_name!')
    end select

    return

end function cst_init_get_mdcon_mode_name

!===============================================================================
! Function:  cst_init_get_mdcon_cvtype_name
!===============================================================================

character(80) function cst_init_get_mdcon_cvtype_name(cvname)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: cvname
    ! --------------------------------------------------------------------------

    select case(cvname)
        case(CON_CVTYPE_DS)
            cst_init_get_mdcon_cvtype_name = "DS"
        case(CON_CVTYPE_DIS)
            cst_init_get_mdcon_cvtype_name = "DIS"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented CV type in cst_init_get_mdcon_cvtype_name!')
    end select

    return

end function cst_init_get_mdcon_cvtype_name

!===============================================================================
! Function:  cst_init_get_lamsol_name
!===============================================================================

character(80) function cst_init_get_lamsol_name(lamsol)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: lamsol
    ! --------------------------------------------------------------------------

    select case(lamsol)
        case(CON_LAMSOL_MD)
            cst_init_get_lamsol_name = "MD Engine"
        case(CON_LAMSOL_V1)
            cst_init_get_lamsol_name = "Lambda Calculation (Explicit)"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented lambda solver in cst_init_get_lamsol_name!')
    end select

    return

end function cst_init_get_lamsol_name

!===============================================================================
! Function:  cst_init_get_icfsol_name
!===============================================================================

character(80) function cst_init_get_icfsol_name(icfsol)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: icfsol
    ! --------------------------------------------------------------------------

    select case(icfsol)
        case(CON_ICFSOL_V1)
            cst_init_get_icfsol_name = "V1 (numeric divergence)"
        case(CON_ICFSOL_V2)
            cst_init_get_icfsol_name = "V2 (analytic with analytic/numeric CV Hessian + symmetry)"
        case(CON_ICFSOL_V3)
            cst_init_get_icfsol_name = "V3 (stochastic divergence)"
        case(CON_ICFSOL_V4)
            cst_init_get_icfsol_name = "V4 (stochastic divergence - mass weighted)"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented ICF solver in cst_init_get_icfsol_name!')
    end select

    return

end function cst_init_get_icfsol_name

!===============================================================================
! Function:  cst_init_get_ekinsrc_name
!===============================================================================

character(80) function cst_init_get_ekinsrc_name(ekinsrc)

    use cst_dat
    use pmf_utils

    implicit none
    integer     :: ekinsrc
    ! --------------------------------------------------------------------------

    select case(ekinsrc)
        case(CON_EKINSRC_VV)
            cst_init_get_ekinsrc_name = "VV (velocity Verlet)"
        case(CON_EKINSRC_V4)
            cst_init_get_ekinsrc_name = "V4 (4th-order)"
        case(CON_EKINSRC_V6)
            cst_init_get_ekinsrc_name = "V6 (6th-order)"
        case default
            call pmf_utils_exit(PMF_OUT, 1, &
                        '[CST] Not implemented ekin srource in cst_init_get_ekinsrc_name!')
    end select

    return

end function cst_init_get_ekinsrc_name

!===============================================================================
! Subroutine:  cst_init_add_shake_csts
!===============================================================================

subroutine cst_init_add_shake_csts

    use pmf_utils
    use pmf_dat
    use cv_ds
    use cv_dis
    use cst_dat
    use cst_constraints
    use pmf_unit

    implicit none
    type(CVPointer),allocatable    :: CVList_backup(:)
    type(CVTypeBM),allocatable     :: CONList_backup(:)
    integer                        :: i,cvid,conid,alloc_failed
    ! ------------------------------------------------------------------------------

    NumOfCONs = NumOfAllCONs

    if( NumOfMDCONs .eq. 0 .or. NumOfAllCONs .eq. 0 ) return

    ! backup old CVs
    allocate(CVList_backup(NumOfCVs),   &
          CONList_backup(NumOfAllCONs), &
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

    allocate(CVList(NumOfCVs + NumOfMDCONs),  &
          CONList(NumOfAllCONs + NumOfMDCONs), &
          stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CVList/CONList array!')
    end if

    do i=1,NumOfCVs
        CVList(i)  = CVList_backup(i)
    end do
    do i=1,NumOfAllCONs
        CONList(i) = CONList_backup(i)
    end do

    deallocate(CVList_backup)
    deallocate(CONList_backup)

    ! add SHAKE constraints
    do i= 1,NumOfMDCONs
        cvid  = NumOfCVs + i
        conid = NumOfAllCONs + i
        ! CV -----------------------------------------

        select case(fmdcon_cvtype)
            case(CON_CVTYPE_DS)
                allocate(CVTypeDS::CVList(cvid)%cv, &
                         stat = alloc_failed)
                if( alloc_failed .ne. 0 ) then
                    call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for SHAKE constraint!')
                end if
                call CVList(cvid)%cv%reset_cv()
                CVList(cvid)%cv%ctype     = 'DS'
                CVList(cvid)%cv%unit      = pmf_unit_power_unit(LengthUnit,2)
                CVList(cvid)%cv%idx       = i
                CVList(cvid)%cv%name      = 'MDCON'
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
                CVList(cvid)%cv%rindexes(1) = MDCONList(i)%at1
                CVList(cvid)%cv%rindexes(2) = MDCONList(i)%at2
                ! CST ----------------------------------------
                call cst_constraints_reset_con(CONList(conid))
                CONList(conid)%cvindx       = cvid
                CONList(conid)%cv           => CVList(cvid)%cv
                CONList(conid)%mode         = 'C'
                CONList(conid)%value        = MDCONList(i)%value
                CONList(conid)%value_set    = .true.
            case(CON_CVTYPE_DIS)
                allocate(CVTypeDIS::CVList(cvid)%cv, &
                         stat = alloc_failed)
                if( alloc_failed .ne. 0 ) then
                    call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for SHAKE constraint!')
                end if
                call CVList(cvid)%cv%reset_cv()
                CVList(cvid)%cv%ctype     = 'DIS'
                CVList(cvid)%cv%unit      = LengthUnit
                CVList(cvid)%cv%idx       = i
                CVList(cvid)%cv%name      = 'MDCON'
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
                CVList(cvid)%cv%rindexes(1) = MDCONList(i)%at1
                CVList(cvid)%cv%rindexes(2) = MDCONList(i)%at2
                ! CST ----------------------------------------
                call cst_constraints_reset_con(CONList(conid))
                CONList(conid)%cvindx       = cvid
                CONList(conid)%cv           => CVList(cvid)%cv
                CONList(conid)%mode         = 'C'
                CONList(conid)%value        = sqrt(MDCONList(i)%value)
                CONList(conid)%value_set    = .true.
            case default
                call pmf_utils_exit(PMF_OUT, 1,'[CST] Unsupported fmdcon_cvtype in cst_init_add_shake_csts!')
        end select
    end do

    ! correct numbers
    NumOfCVs     = NumOfCVs + NumOfMDCONs
    NumOfAllCONs = NumOfCONs + NumOfMDCONs

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

    do i=1, NumOfAllCONs
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
    do i=1,NumOfAllCONs
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
    use cst_accu

    implicit none
    integer      :: alloc_failed, ntau
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
    allocate(lambda(NumOfAllCONs),              &
             cv(NumOfAllCONs),                  &
             vv(NumOfAllCONs),                  &
             indx(NumOfAllCONs),                &
             zmat(NumOfAllCONs,NumOfAllCONs),   &
             zmats(NumOfMDCONs,NumOfMDCONs), stat= alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in LU decomposition!')
    end if

    if( fshakesolver .eq. CON_SHAKESOL_NMSVD ) then
        ! allocate arrays for SVD decomposition
        lsvdwork = (3*NumOfAllCONs + max( 2*NumOfAllCONs, NumOfAllCONs, 1 ))*10
        allocate(svdwork(lsvdwork), stat= alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used in LU decomposition!')
        end if
    end if

    if( fint_der ) then
        ! allocate arrays for matrix inversion
        linvwork = NumOfAllCONs * 64
        allocate(invwork(linvwork), stat= alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used in matrix inversion!')
        end if
        ! for QR normalization
        lsdivwork = max(3*NumOfLAtoms,fpmf_sdiv_S) * 64
        ntau = min(3*NumOfLAtoms,fpmf_sdiv_S)
        allocate(sdivwork(lsdivwork),sdiv_z(3,NumOfLAtoms,fpmf_sdiv_S),sdivtau(ntau), stat= alloc_failed)
        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used in QR normalization!')
        end if
    end if

! allocate arrays for lambda calculation
    select case(fintalg)
        case(IA_LEAP_FROG) ! FIXME
            allocate(lambdax(NumOfAllCONs), stat= alloc_failed )
            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used in lambda calculation!')
            end if
            lambdax(:) = 0.0d0
        case(IA_VEL_VERLET)
            allocate(lambdax(NumOfAllCONs), lambdav(NumOfAllCONs), stat= alloc_failed )
            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used in lambda calculation!')
            end if
            lambdax(:) = 0.0d0
            lambdav(:) = 0.0d0
        case(IA_LF_MIDDLE)
            allocate(lambdax(NumOfAllCONs), lambdav(NumOfAllCONs), stat= alloc_failed )
            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used in lambda calculation!')
            end if
            lambdax(:) = 0.0d0
            lambdav(:) = 0.0d0
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_init_print_summary!')
    end select

! history buffers
    hist_len = 7       ! FIXED at 7
    hist_fidx = -3
    hist_fidx_tds = -3

    allocate( lambdaMhist(NumOfAllCONs,hist_len),   &
              lambdaEhist(NumOfAllCONs,hist_len),   &
              fwhist(hist_len),                     &
              epothist(hist_len),                   &
              ersthist(hist_len),                   &
              ekinhist(hist_len),                   &
              enevalidhist(hist_len),               &
              icfphist(NumOfAllCONs,hist_len),      &
              icfkhist(NumOfAllCONs,hist_len),      &
              crdhist(3,NumOfLAtoms,hist_len),               &
              cvderhist(3,NumOfLAtoms,NumOfCVs,hist_len),    &
              frchist(3,NumOfLAtoms,hist_len),               &
              velhist(3,NumOfLAtoms,hist_len),               &
              stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used for history recording!')
    end if

    lambdaMhist(:,:)    = 0.0d0
    lambdaEhist(:,:)    = 0.0d0
    fwhist(:)           = 0.0d0

    epothist(:)         = 0.0d0
    ersthist(:)         = 0.0d0
    ekinhist(:)         = 0.0d0
    enevalidhist(:)     = .false.

    icfphist(:,:)       = 0.0d0
    icfkhist(:,:)       = 0.0d0

    cvderhist(:,:,:,:)  = 0.0d0
    frchist(:,:,:)      = 0.0d0
    velhist(:,:,:)      = 0.0d0

! internal energy/entropy
    if( fintcalc .and. fint_der ) then
        allocate( icf_he(3,NumOfLAtoms),    &
                  icf_vi1(3,NumOfLAtoms),   &
                  icf_vi2(3,NumOfLAtoms),   &
                  icf_vin(3,NumOfLAtoms,NumOfAllCONs),    &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for internal energy/entropy calculation!')
        end if
    end if

! initialize PMFAccu
    call cst_accu_alloc

    return

end subroutine cst_init_core

!===============================================================================

end module cst_init
