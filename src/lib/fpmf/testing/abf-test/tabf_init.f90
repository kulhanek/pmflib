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

module tabf_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  tabf_init_method
!===============================================================================

subroutine tabf_init_method

    use tabf_output
    use tabf_restart
    use tabf_trajectory

    implicit none
    ! --------------------------------------------------------------------------

    call tabf_init_arrays
    call tabf_init_print_header
    call tabf_output_open
    call tabf_trajectory_open
    call tabf_output_write_header

end subroutine tabf_init_method

!===============================================================================
! Subroutine:  tabf_init_dat
!===============================================================================

subroutine tabf_init_dat

    use tabf_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0         ! 0 - disable ABF, 1 - enabled ABF
    fsample         = 5000      ! output sample period in steps
    frstupdate      = 5000      ! how often is restart file written
    feimode         = 1         ! extrapolation / interpolation mode
                                ! 1 - linear ramp I
    ftrjsample      = 0         ! how often save accumulator to "accumulator evolution"
    fapply_abf      = .true.    ! on - apply ABF, off - do not apply ABF
    fprint_icf      = .false.   ! T - print instantaneous collective forces, F - do not print

    fenthalpy       = .false.   ! accumulate enthalpy
    fentropy        = .false.   ! accumulate entropy
    fepotoffset     = 0.0d0
    fekinoffset     = 0.0d0

    fblock_size     = 0

    fhramp_min      = 100       ! definition of linear ramp
    fhramp_max      = 200       ! definition of linear ramp

    NumOfTABFCVs     = 0

    insidesamples   = 0
    outsidesamples  = 0

    fdtx            = 0.0d0     ! time step in internal units

end subroutine tabf_init_dat

!===============================================================================
! Subroutine:  tabf_print_header
!===============================================================================

subroutine tabf_init_print_header

    use prmfile
    use pmf_constants
    use pmf_dat
    use tabf_dat
    use tabf_cvs
    use pmf_utils

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ****************** ADAPTIVE BIASING FORCE METHOD (TESTING) ******************* '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' ABF mode (fmode)                        : ', fmode
    select case(fmode)
    case(1)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm'
    case(2)
    write(PMF_OUT,120)  '      |-> Original ABF algorithm'
    case(3)
    write(PMF_OUT,120)  '      |-> Analytical/Numerical ABF algorithm, SHAKE must be off'
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown fmode in tabf_init_print_header!')
    end select
    write(PMF_OUT,125)  ' Coordinate definition file (ftabfdef)   : ', trim(ftabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfTABFCVs
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Apply ABF force (fapply_abf)            : ', prmfile_onoff(fapply_abf)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Interpolation/Extrapolation '
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Data pre-blocking (fblock_size)         : ', fblock_size
    write(PMF_OUT,130)  ' Extra/interpolation mode (feimode)      : ', feimode
    select case(feimode)
    case(0)
    write(PMF_OUT,120)  '      |-> Direct application'
    case(1)
    write(PMF_OUT,120)  '      |-> Min/Max linear ramp'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode in tabf_init_print_header!')
    end select

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (ftabfrst)                 : ', trim(ftabfrst)
    write(PMF_OUT,130)  ' Final restart update (frstupdate)       : ', frstupdate
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (ftabfout)                  : ', trim(ftabfout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,125)  ' Print instantaneous forces (fprint_icf) : ', prmfile_onoff(fprint_icf)
    write(PMF_OUT,125)  ' ICF file (ftabficf)                     : ', trim(ftabficf)
    write(PMF_OUT,125)  ' Accumulate enthalpy (fenthalpy)         : ', prmfile_onoff(fenthalpy)
    write(PMF_OUT,125)  ' Accumulate entropy (fentropy)           : ', prmfile_onoff(fentropy)
    write(PMF_OUT,150)  ' Potential energy offset (fepotoffset)   : ', pmf_unit_get_rvalue(EnergyUnit,fepotoffset),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,150)  ' Kinetic energy offset (fekinoffset)     : ', pmf_unit_get_rvalue(EnergyUnit,fekinoffset), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample
    write(PMF_OUT,125)  ' Trajectory file (ftabftrj)              : ', trim(ftabftrj)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of TABF collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfTABFCVs
        write(PMF_OUT,140) i
        call tabf_cvs_cv_info(TABFCVList(i))
        write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

 return

120 format(A)
125 format(A,A)
130 format(A,I6)
150 format(A,F10.1,1X,A)

140 format(' == Collective variable #',I4.4)

end subroutine tabf_init_print_header

!===============================================================================
! Subroutine:  tabf_init_arrays
!===============================================================================

subroutine tabf_init_arrays

    use pmf_utils
    use pmf_dat
    use tabf_dat
    use tabf_accu

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

    fdtx = fdt*PMF_DT2VDT

    ! general arrays --------------------------------
    allocate(                               &
            a0(3,NumOfLAtoms),              &
            a1(3,NumOfLAtoms),              &
            v0(3,NumOfLAtoms),              &
            pxi0(NumOfTABFCVs),              &
            pxi1(NumOfTABFCVs),              &
            pxip(NumOfTABFCVs),              &
            pxim(NumOfTABFCVs),              &
            pdum(NumOfTABFCVs),              &
            avg_values(NumOfTABFCVs),        &
            la(NumOfTABFCVs),                &
            fz(NumOfTABFCVs,NumOfTABFCVs),    &
            fzinv(NumOfTABFCVs,NumOfTABFCVs), &
            zd0(3,NumOfLAtoms,NumOfTABFCVs), &
            zd1(3,NumOfLAtoms,NumOfTABFCVs), &
            cvaluehist0(NumOfTABFCVs),       &
            cvaluehist1(NumOfTABFCVs),       &
            cvaluehist2(NumOfTABFCVs),       &
            cvaluehist3(NumOfTABFCVs),       &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays used in ABF calculation!')
    end if

    a0(:,:) = 0.0d0
    a1(:,:) = 0.0d0
    v0(:,:) = 0.0d0
    pxi0(:) = 0.0d0
    pxi1(:) = 0.0d0
    pxip(:) = 0.0d0
    pxim(:) = 0.0d0
    pdum(:) = 0.0d0
    avg_values(:) = 0.0d0
    la(:) = 0.0d0
    fz(:,:) = 0.0d0
    fzinv(:,:) = 0.0d0
    zd0(:,:,:) = 0.0d0
    zd1(:,:,:) = 0.0d0
    cvaluehist0(:) = 0.0d0
    cvaluehist1(:) = 0.0d0
    cvaluehist2(:) = 0.0d0
    cvaluehist3(:) = 0.0d0

    ! for Z matrix inversion, only if fnitem > 1 ----
    if( NumOfTABFCVs .gt. 1 ) then
        allocate( vv(NumOfTABFCVs),               &
                  indx(NumOfTABFCVs),             &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[ABF] Unable to allocate memory for arrays used in Z matrix inversion!')
        end if
    end if

    ! init accumulator ------------------------------
    call tabf_accu_init

end subroutine tabf_init_arrays

!===============================================================================

end module tabf_init
