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

    fsmoothetot     = .false.
    fcontbias       = .true.

    NumOfUSABFCVs     = 0

    insidesamples   = 0
    outsidesamples  = 0

    gpr_len         = 7
    gpr_width       = 6.0
    gpr_noise       = 0.05

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
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (2-points ABF)'
    case(2)
    write(PMF_OUT,120)  '      |-> 7-points ABF'
    case(3)
    write(PMF_OUT,120)  '      |-> 10-points ABF'
    case(4)
    write(PMF_OUT,120)  '      |-> GPR ABF'
    write(PMF_OUT,130)  '          gpr_len                        : ', gpr_len
    write(PMF_OUT,145)  '          gpr_width                      : ', gpr_width
    write(PMF_OUT,145)  '          gpr_noise                      : ', gpr_noise
    case(5)
    write(PMF_OUT,120)  '      |-> GPR ABF II'
    write(PMF_OUT,130)  '          gpr_len                        : ', gpr_len
    write(PMF_OUT,145)  '          gpr_width                      : ', gpr_width
    write(PMF_OUT,145)  '          gpr_noise                      : ', gpr_noise
    case default
    call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unknown fmode in usabf_init_print_header!')
    end select
    write(PMF_OUT,125)  ' Coordinate definition file (fusabfdef)  : ', trim(fusabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfUSABFCVs
    write(PMF_OUT,125)  ' Use continuous US bias (fcontbias)      : ', prmfile_onoff(fcontbias)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fusabfrst)                : ', trim(fusabfrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Final restart update (frstupdate)       : ', frstupdate

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Enthalpy/entropy options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Accumulate enthalpy (fenthalpy)         : ', prmfile_onoff(fenthalpy)
    write(PMF_OUT,125)  ' Accumulate entropy (fentropy)           : ', prmfile_onoff(fentropy)
    write(PMF_OUT,125)  ' Smooth Etot (fsmoothetot)               : ', prmfile_onoff(fsmoothetot)
    write(PMF_OUT,150)  ' Potential energy offset (fepotaverage)  : ', pmf_unit_get_rvalue(EnergyUnit,fepotaverage),  &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,150)  ' Kinetic energy offset (fekinaverage)    : ', pmf_unit_get_rvalue(EnergyUnit,fekinaverage), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fusabfout)                 : ', trim(fusabfout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample

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
145 format(A,F10.3)
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

    select case(fmode)
        case(1)
            hist_len = 2
        case(2)
            hist_len = 7
        case(3)
            hist_len = 10
        case(4,5)
            call usabf_init_gpr()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Not implemented fmode in usabf_init_arrays!')
    end select

    ! general arrays --------------------------------
    allocate(                                   &
            a0(3,NumOfLAtoms),                  &
            v0(3,NumOfLAtoms),                  &
            pxi0(NumOfUSABFCVs),                &
            pxi1(NumOfUSABFCVs),                &
            pxip(NumOfUSABFCVs),                &
            pxim(NumOfUSABFCVs),                &
            la(NumOfUSABFCVs),                  &
            fz(NumOfUSABFCVs,NumOfUSABFCVs),    &
            fzinv(NumOfUSABFCVs,NumOfUSABFCVs), &
            zd0(3,NumOfLAtoms,NumOfUSABFCVs),   &
            zd1(3,NumOfLAtoms,NumOfUSABFCVs),   &
            cvhist(NumOfUSABFCVs,hist_len),     &
            pcvhist(NumOfUSABFCVs,hist_len),    &
            epothist(hist_len),                 &
            etothist(hist_len),                 &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[US-ABF] Unable to allocate memory for arrays used in ABF calculation in usabf_init_arrays!')
    end if

    a0(:,:) = 0.0d0
    v0(:,:) = 0.0d0
    pxi0(:) = 0.0d0
    pxi1(:) = 0.0d0
    pxip(:) = 0.0d0
    pxim(:) = 0.0d0
    la(:) = 0.0d0
    fz(:,:) = 0.0d0
    fzinv(:,:) = 0.0d0
    zd0(:,:,:) = 0.0d0
    zd1(:,:,:) = 0.0d0

    cvhist(:,:) = 0.0d0
    pcvhist(:,:) = 0.0d0
    epothist(:) = 0.0d0
    etothist(:) = 0.0d0

    ! for Z matrix inversion, only if fnitem > 1 ----
    if( NumOfUSABFCVs .gt. 1 ) then
        allocate( vv(NumOfUSABFCVs),               &
                  indx(NumOfUSABFCVs),             &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[US-ABF] Unable to allocate memory for arrays used in Z matrix inversion in usabf_init_arrays!')
        end if
    end if

    ! init accumulator ------------------------------
    call usabf_accu_init

end subroutine usabf_init_arrays

!===============================================================================
! Subroutine:  usabf_init_gpr
!===============================================================================

subroutine usabf_init_gpr

    use pmf_utils
    use pmf_dat
    use usabf_dat
    use usabf_accu

    implicit none
    integer     :: i,j,alloc_failed
    real(PMFDP) :: dt,dt2,idw2
    ! --------------------------------------------------------------------------

    idw2 = 1.0d0/(2.0d0*(gpr_width * PMF_DT2VDT)**2)

! gpr_len must be an odd number
    if( mod(gpr_len,2) .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'[US-ABF] gpr_len must be an odd number in usabf_init_gpr!')
    end if

! allocate arrays
    allocate(                           &
            gpr_K(gpr_len,gpr_len),     &
            gpr_model(gpr_len),         &
            gpr_kff(gpr_len),           &
            gpr_kdf(gpr_len),           &
            gpr_indx(gpr_len),          &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[US-ABF] Unable to allocate memory for GPR arrays in usabf_init_gpr!')
    end if

! init covariance matrix
    do i=1,gpr_len
        do j=1,gpr_len
            dt2 = (real(i-j,PMFDP) * fdtx)**2
            gpr_K(i,j) = exp(- dt2 * idw2)
        end do
    end do

    do i=1,gpr_len
        gpr_K(i,i) = gpr_K(i,i) + gpr_noise
    end do

! run LU decomposition
    call dgetrf(gpr_len,gpr_len,gpr_K,gpr_len,gpr_indx,gpr_info)

    if( gpr_info .ne. 0 ) then
        ! throw error
        call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to run LU decomposition in usabf_init_gpr!')
    end if

! construct kff
    j = gpr_len / 2 + 1
    do i=1,gpr_len
        dt = real(i-j,PMFDP) * fdtx
        gpr_kdf(i) = - 2.0d0 * exp(- dt**2 * idw2) * dt * idw2
        gpr_kff(i) = exp(- dt**2 * idw2)
    end do

    hist_len = 0

    if( fmode .eq. 4 ) then
        hist_len = gpr_len + gpr_len/2
    end if

    if( fmode .eq. 5 ) then
        hist_len = gpr_len + 1

        allocate(                                           &
                cvcontex0%CVsValues(NumOfCVs),              &
                cvcontex0%CVsDrvs(3,NumOfLAtoms,NumOfCVs),  &
                stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                '[US-ABF] Unable to allocate memory for GPR arrays in usabf_init_gpr!')
        end if

    end if

end subroutine usabf_init_gpr

!===============================================================================

end module usabf_init
