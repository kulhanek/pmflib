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

module abf_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_init_method
!===============================================================================

subroutine abf_init_method

    use abf_output
    use abf_restart
    use abf_trajectory
    use abf_client

    implicit none
    ! --------------------------------------------------------------------------

    call abf_init_arrays
    call abf_init_print_header
    call abf_output_open
    call abf_restart_read
    call abf_trajectory_open
    call abf_client_register
    call abf_client_get_initial_data
    call abf_output_write_header

end subroutine abf_init_method

!===============================================================================
! Subroutine:  abf_init_dat
!===============================================================================

subroutine abf_init_dat

    use abf_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0
    fsample         = 5000
    frestart        = .false.
    frstupdate      = 5000
    ftrjsample      = 0
    fapply_mask     = .false.
    fapply_abf      = .true.

    fenthalpy       = .false.
    fentropy        = .false.
    fsmoothetot     = .false.

    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0

    feimode         = 1
    fhramp_min      = 100
    fhramp_max      = 200

    NumOfABFCVs     = 0

    fserver_enabled = .false.
    fserverkey      = ''
    fserver         = ''
    fserverupdate   = 20000
    fconrepeats     = 0
    fabortonmwaerr  = .true.
    fmwamode        = 0

    client_id       = -1
    failure_counter = 0

    insidesamples   = 0
    outsidesamples  = 0

    gpr_len         = 7
    gpr_width       = 4.0
    gpr_noise       = 0.01
    gpr_kernel      = 2

    fdtx            = 0.0d0

end subroutine abf_init_dat

!===============================================================================
! Subroutine:  abf_print_header
!===============================================================================

subroutine abf_init_print_header

    use prmfile
    use pmf_constants
    use pmf_dat
    use abf_dat
    use abf_cvs
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
    write(PMF_OUT,120)  '      |-> GPR ABF algorithm'
    write(PMF_OUT,130)  '          gpr_len                        : ', gpr_len
    write(PMF_OUT,145)  '          gpr_width                      : ', gpr_width
    write(PMF_OUT,145)  '          gpr_noise                      : ', gpr_noise

    select case(gpr_kernel)
        case(0)
    write(PMF_OUT,125)  '          gpr_kernel                      : ', '0 - Matern class v=3/2'
        case(1)
    write(PMF_OUT,125)  '          gpr_kernel                      : ', '1 - Matern class v=5/2'
        case(2)
    write(PMF_OUT,125)  '          gpr_kernel                      : ', '2 - Squared exponential'

        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_kernel in abf_init_print_header!')
    end select
    case(4)
    write(PMF_OUT,120)  '      |-> GPR ABF algorithm (CV momenta)'
    write(PMF_OUT,130)  '          gpr_len                        : ', gpr_len
    write(PMF_OUT,145)  '          gpr_width                      : ', gpr_width
    write(PMF_OUT,145)  '          gpr_noise                      : ', gpr_noise

    select case(gpr_kernel)
        case(0)
    write(PMF_OUT,125)  '          gpr_kernel                      : ', '0 - Matern class v=3/2'
        case(1)
    write(PMF_OUT,125)  '          gpr_kernel                      : ', '1 - Matern class v=5/2'
        case(2)
    write(PMF_OUT,125)  '          gpr_kernel                      : ', '2 - Squared exponential'

        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_kernel in abf_init_print_header!')
    end select
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown fmode in abf_init_print_header!')
    end select
    write(PMF_OUT,125)  ' Coordinate definition file (fabfdef)    : ', trim(fabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfABFCVs
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Apply ABF force (fapply_abf)            : ', prmfile_onoff(fapply_abf)
    write(PMF_OUT,125)  ' ABF mask mode (fapply_mask)             : ', prmfile_onoff(fapply_mask)
    write(PMF_OUT,125)  ' ABF mask file (fabfmask)                : ', trim(fabfmask)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Interpolation/Extrapolation '
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Extra/interpolation mode (feimode)      : ', feimode
    select case(feimode)
    case(0)
    write(PMF_OUT,120)  '      |-> Disabled'
    case(1)
    write(PMF_OUT,120)  '      |-> Min/Max linear ramp'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode in abf_init_print_header!')
    end select

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
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fabfrst)                  : ', trim(fabfrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Restart file update (frstupdate)        : ', frstupdate

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fabfout)                   : ', trim(fabfout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample
    write(PMF_OUT,125)  ' Trajectory file (fabftrj)               : ', trim(fabftrj)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' MWA server options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Server communication is                      : ', prmfile_onoff(fserver_enabled)
    if( fserver_enabled ) then
    write(PMF_OUT,125)  ' Server key file name (fserverkey)            : ', trim(fserverkey)
    else
    write(PMF_OUT,125)  ' Server key file name (fserverkey)            : ', 'none'
    end if
    write(PMF_OUT,130)  ' Server update interval (fserverupdate)       : ', fserverupdate
    write(PMF_OUT,130)  ' Number of connection repeats (fconrepeats)   : ', fconrepeats
    write(PMF_OUT,125)  ' Abort on MWA failure (fabortonmwaerr)        : ', prmfile_onoff(fabortonmwaerr)
    write(PMF_OUT,130)  ' MWA mode (fmwamode)                          : ', fmwamode
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of ABF collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfABFCVs
    write(PMF_OUT,140) i
    call abf_cvs_cv_info(ABFCVList(i))
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

end subroutine abf_init_print_header

!===============================================================================
! Subroutine:  abf_init_arrays
!===============================================================================

subroutine abf_init_arrays

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

    fdtx = fdt*PMF_DT2VDT

    select case(fmode)
        case(1)
            hist_len = 2
        case(2)
            hist_len = 4
        case(3)
            call abf_init_gpr()
        case(4)
            call abf_init_gpr()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_init_arrays!')
    end select

    ! general arrays --------------------------------
    allocate(                               &
            a1(3,NumOfLAtoms),              &
            a0(3,NumOfLAtoms),              &
            v0(3,NumOfLAtoms),              &
            cvcontex0%CVsValues(NumOfCVs),              &
            cvcontex0%CVsDrvs(3,NumOfLAtoms,NumOfCVs),  &
            la(NumOfABFCVs),                &
            zd0(3,NumOfLAtoms,NumOfABFCVs), &
            zd1(3,NumOfLAtoms,NumOfABFCVs), &
            pxi0(NumOfABFCVs),              &
            pxi1(NumOfABFCVs),              &
            pxip(NumOfABFCVs),              &
            pxim(NumOfABFCVs),              &
            cvave(NumOfABFCVs),             &
            fz(NumOfABFCVs,NumOfABFCVs),    &
            fzinv(NumOfABFCVs,NumOfABFCVs), &
            cvhist(NumOfABFCVs,hist_len),   &
            pchist(NumOfABFCVs,hist_len),   &
            epothist(hist_len),             &
            etothist(hist_len),             &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays used in ABF calculation!')
    end if

    a1(:,:)     = 0.0d0
    a0(:,:)     = 0.0d0
    v0(:,:)     = 0.0d0

    la(:)       = 0.0d0
    zd0(:,:,:)  = 0.0d0
    zd1(:,:,:)  = 0.0d0
    pxi0(:)     = 0.0d0
    pxi1(:)     = 0.0d0
    pxip(:)     = 0.0d0
    pxim(:)     = 0.0d0

    cvave(:)    = 0.0d0

    fz(:,:)     = 0.0d0
    fzinv(:,:)  = 0.0d0

    cvhist(:,:) = 0.0d0
    pchist(:,:) = 0.0d0
    epothist(:) = 0.0d0
    etothist(:) = 0.0d0

    ! for Z matrix inversion, only if fnitem > 1 ----
    if( NumOfABFCVs .gt. 1 ) then
        allocate( vv(NumOfABFCVs),               &
                  indx(NumOfABFCVs),             &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[ABF] Unable to allocate memory for arrays used in Z matrix inversion!')
        end if
    end if

    ! init accumulator ------------------------------
    call abf_accu_init

end subroutine abf_init_arrays

!===============================================================================
! Subroutine:  abf_init_gpr
!===============================================================================

subroutine abf_init_gpr

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer     :: i,j,alloc_failed
    real(PMFDP) :: r
    ! --------------------------------------------------------------------------

! gpr_len must be an odd number
    if( mod(gpr_len,2) .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'[US-ABF] gpr_len must be an odd number in abf_init_gpr!')
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
            '[US-ABF] Unable to allocate memory for GPR arrays in abf_init_gpr!')
    end if

! init covariance matrix
    do i=1,gpr_len
        do j=1,gpr_len
            r = abs(real(i-j,PMFDP)/gpr_width)
            select case(gpr_kernel)
                case(0)
                    gpr_K(i,j) = (1.0d0 + sqrt(3.0) * r) * exp(- sqrt(3.0d0) * r)
                case(1)
                    gpr_K(i,j) = (1.0d0 + sqrt(5.0) * r + 5.0d0/3.0d0 * r**2) * exp(- sqrt(5.0d0) * r)
                case(2)
                    gpr_K(i,j) = exp(- 0.5d0 * r**2)
                case default
                    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_kernel in abf_init_gpr I!')
            end select
        end do
    end do

    do i=1,gpr_len
        gpr_K(i,i) = gpr_K(i,i) + gpr_noise
    end do

! run LU decomposition
    call dgetrf(gpr_len,gpr_len,gpr_K,gpr_len,gpr_indx,gpr_info)

    if( gpr_info .ne. 0 ) then
        ! throw error
        call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to run LU decomposition in abf_init_gpr!')
    end if

! construct kff
    j = gpr_len / 2 + 1
    do i=1,gpr_len
        r = abs(real(i-j,PMFDP)/gpr_width)
        select case(gpr_kernel)
            case(0)
                gpr_kdf(i) = 3.0d0 * r * exp(- sqrt(3.0d0) * r) / (gpr_width * fdtx) * sign(1.0d0,real(i-j,PMFDP))
                gpr_kff(i) = (1.0d0 + sqrt(3.0) * r) * exp(- sqrt(3.0d0) * r)
            case(1)
                gpr_kdf(i) = 5.0d0/3.0d0 * r * (1.0d0 + sqrt(5.0d0) * r) * exp(- sqrt(5.0d0) * r ) / (gpr_width * fdtx) &
                           * sign(1.0d0,real(i-j,PMFDP))
                gpr_kff(i) = (1.0d0 + sqrt(5.0d0) * r + 5.0d0/3.0d0 * r**2) * exp(- sqrt(5.0d0) * r)
            case(2)
                gpr_kdf(i) = exp(- 0.5d0 * r**2) * r / (gpr_width * fdtx) * sign(1.0d0,real(i-j,PMFDP))
                gpr_kff(i) = exp(- 0.5d0 * r**2)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_kernel in abf_init_gpr II!')
        end select
    end do

    hist_len = 0

    if( fmode .eq. 3 ) then
        hist_len = gpr_len + gpr_len/2
    end if

    if( fmode .eq. 4) then
        ! with CV momenta
        hist_len = gpr_len + 1
    end if

end subroutine abf_init_gpr

!===============================================================================

end module abf_init
