!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module usabf_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  usabf_init_method
!===============================================================================

subroutine usabf_init_method

    use usabf_output
    use usabf_restart
    use usabf_trajectory

    implicit none
    ! --------------------------------------------------------------------------

    call usabf_init_arrays
    call usabf_init_print_header
    call usabf_output_open
    call usabf_restart_read
    call usabf_trajectory_open
    call usabf_output_write_header

end subroutine usabf_init_method

!===============================================================================
! Subroutine:  usabf_init_dat
!===============================================================================

subroutine usabf_init_dat

    use usabf_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0         ! 0 - disable ABF, 1 - enabled ABF

    frestart        = .false.
    faccurst        = -1

    fsample         = 5000      ! output sample period in steps
    frstupdate      = 5000      ! how often is restart file written

    ftrjsample      = 0         ! how often save accumulator to "accumulator evolution"

    fenthalpy       = .false.   ! accumulate enthalpy
    fentropy        = .false.   ! accumulate entropy
    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0

    NumOfUSABFCVs     = 0

    insidesamples   = 0
    outsidesamples  = 0

    fdtx            = 0.0d0     ! time step in internal units

end subroutine usabf_init_dat

!===============================================================================
! Subroutine:  usabf_print_header
!===============================================================================

subroutine usabf_init_print_header

    use prmfile
    use pmf_constants
    use pmf_dat
    use usabf_dat
    use usabf_cvs
    use pmf_utils

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' *********************** ADAPTIVE BIASING FORCE METHOD ************************ '
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
    write(PMF_OUT,120)  '      |-> 7-points ABF'
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown fmode in usabf_init_print_header!')
    end select
    write(PMF_OUT,125)  ' Coordinate definition file (fusabfdef)  : ', trim(fusabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfUSABFCVs

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fusabfrst)                : ', trim(fusabfrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Final restart update (frstupdate)       : ', frstupdate
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fusabfout)                 : ', trim(fusabfout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,125)  ' Accumulate enthalpy (fenthalpy)         : ', prmfile_onoff(fenthalpy)
    write(PMF_OUT,125)  ' Accumulate entropy (fentropy)           : ', prmfile_onoff(fentropy)
    write(PMF_OUT,150)  ' Potential energy offset (fepotaverage)  : ', pmf_unit_get_rvalue(EnergyUnit,fepotaverage),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,150)  ' Kinetic energy offset (fekinaverage)    : ', pmf_unit_get_rvalue(EnergyUnit,fekinaverage), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample
    write(PMF_OUT,125)  ' Trajectory file (fusabftrj)             : ', trim(fusabftrj)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of US-ABF collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfUSABFCVs
        write(PMF_OUT,140) i
        call usabf_cvs_cv_info(USABFCVList(i))
        write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

 return

120 format(A)
125 format(A,A)
130 format(A,I6)
150 format(A,F10.1,1X,A)

140 format(' == Collective variable #',I4.4)

end subroutine usabf_init_print_header

!===============================================================================
! Subroutine:  usabf_init_arrays
!===============================================================================

subroutine usabf_init_arrays

    use pmf_utils
    use pmf_dat
    use usabf_dat
    use usabf_accu

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

    fdtx = fdt*PMF_DT2VDT

    ! general arrays --------------------------------
    allocate(                                   &
            a0(3,NumOfLAtoms),                  &
            a1(3,NumOfLAtoms),                  &
            v0(3,NumOfLAtoms),                  &
            pxi0(NumOfUSABFCVs),                &
            pxi1(NumOfUSABFCVs),                &
            pxip(NumOfUSABFCVs),                &
            pxim(NumOfUSABFCVs),                &
            avg_values(NumOfUSABFCVs),          &
            la(NumOfUSABFCVs),                  &
            fz(NumOfUSABFCVs,NumOfUSABFCVs),    &
            fzinv(NumOfUSABFCVs,NumOfUSABFCVs), &
            zd0(3,NumOfLAtoms,NumOfUSABFCVs),   &
            zd1(3,NumOfLAtoms,NumOfUSABFCVs),   &
            cvhist0(NumOfUSABFCVs),             &
            cvhist1(NumOfUSABFCVs),             &
            cvhist2(NumOfUSABFCVs),             &
            cvhist3(NumOfUSABFCVs),             &
            cvhist4(NumOfUSABFCVs),             &
            cvhist5(NumOfUSABFCVs),             &
            cvhist6(NumOfUSABFCVs),             &
            pcvhist0(NumOfUSABFCVs),            &
            pcvhist1(NumOfUSABFCVs),            &
            pcvhist2(NumOfUSABFCVs),            &
            pcvhist3(NumOfUSABFCVs),            &
            pcvhist4(NumOfUSABFCVs),            &
            icf2(NumOfUSABFCVs),                 &
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
    avg_values(:) = 0.0d0
    la(:) = 0.0d0
    fz(:,:) = 0.0d0
    fzinv(:,:) = 0.0d0
    zd0(:,:,:) = 0.0d0
    zd1(:,:,:) = 0.0d0

    cvhist0(:) = 0.0d0
    cvhist1(:) = 0.0d0
    cvhist2(:) = 0.0d0
    cvhist3(:) = 0.0d0
    cvhist4(:) = 0.0d0
    cvhist5(:) = 0.0d0
    cvhist6(:) = 0.0d0

    pcvhist0(:) = 0.0d0
    pcvhist1(:) = 0.0d0
    pcvhist2(:) = 0.0d0
    pcvhist3(:) = 0.0d0
    pcvhist4(:) = 0.0d0

    icf2(:) = 0.0d0

    epothist0 = 0.0d0
    epothist1 = 0.0d0
    epothist2 = 0.0d0
    epothist3 = 0.0d0
    epothist4 = 0.0d0
    epothist5 = 0.0d0
    epothist6 = 0.0d0

    etothist0 = 0.0d0
    etothist1 = 0.0d0
    etothist2 = 0.0d0
    etothist3 = 0.0d0
    etothist4 = 0.0d0
    etothist5 = 0.0d0
    etothist6 = 0.0d0

    ! for Z matrix inversion, only if fnitem > 1 ----
    if( NumOfUSABFCVs .gt. 1 ) then
        allocate( vv(NumOfUSABFCVs),               &
                  indx(NumOfUSABFCVs),             &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[ABF] Unable to allocate memory for arrays used in Z matrix inversion!')
        end if
    end if

    ! init accumulator ------------------------------
    call usabf_accu_init

end subroutine usabf_init_arrays

!===============================================================================

end module usabf_init
