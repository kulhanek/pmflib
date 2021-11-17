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
    frecord         = .false.

    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0

    feimode         = 1
    fhramp_min      = 100
    fhramp_max      = 500

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

    fsmooth_kernel  = 0

    fsgframelen     = 5
    fsgorder        = 3
    fsgsmoothall    = .false.

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
    write(PMF_OUT,120)  '      |-> Savitzky-Golay differentiation ABF algorithm'
    write(PMF_OUT,130)  '          Frame length (fsgframelen)     : ', fsgframelen
    write(PMF_OUT,130)  '          Polynomial order (fsgorder)    : ', fsgorder
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
    case(2)
    write(PMF_OUT,120)  '      |-> Kernel smoother'
    write(PMF_OUT,130)  '          Kernel type (fsmooth_kernel)   : ', fsmooth_kernel
    select case(fsmooth_kernel)
    case(0)
    write(PMF_OUT,120)  '          |-> Epanechnikov (parabolic)'
    case(1)
    write(PMF_OUT,120)  '          |-> Triweight'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown kernel in abf_init_print_header!')
    end select
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case(3)
    write(PMF_OUT,120)  '      |-> Linear interpolation'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode in abf_init_print_header!')
    end select

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Enthalpy/entropy options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Record time progress (frecord)          : ', prmfile_onoff(frecord)
    write(PMF_OUT,125)  ' Accumulate enthalpy (fenthalpy)         : ', prmfile_onoff(fenthalpy)
    write(PMF_OUT,125)  ' Accumulate entropy (fentropy)           : ', prmfile_onoff(fentropy)
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
    integer     :: alloc_failed, i
    ! --------------------------------------------------------------------------

    fdtx = fdt*PMF_DT2VDT

    NumOfBiasedABFCVs = 0
    do i=1,NumOfABFCVs
        if( ABFCVList(i)%nobias ) cycle
        NumOfBiasedABFCVs = NumOfBiasedABFCVs + 1
    end do

! general arrays --------------------------------
    allocate(                                   &
            a1(3,NumOfLAtoms),                  &
            a0(3,NumOfLAtoms),                  &
            v0(3,NumOfLAtoms),                  &
            la(NumOfABFCVs),                    &
            zd0(3,NumOfLAtoms,NumOfABFCVs),     &
            zd1(3,NumOfLAtoms,NumOfABFCVs),     &
            pxi0(NumOfABFCVs),                  &
            pxi1(NumOfABFCVs),                  &
            pxip(NumOfABFCVs),                  &
            pxim(NumOfABFCVs),                  &
            cvave(NumOfABFCVs),                 &
            cvcur(NumOfABFCVs),                 &
            cv1dr(NumOfABFCVs),                 &
            cv2dr(NumOfABFCVs),                 &
            fz(NumOfABFCVs,NumOfABFCVs),        &
            fzinv(NumOfABFCVs,NumOfABFCVs),     &
            fzinv0(NumOfABFCVs,NumOfABFCVs),    &
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
    cvcur(:)    = 0.0d0
    cv1dr(:)    = 0.0d0
    cv2dr(:)    = 0.0d0

    fz(:,:)     = 0.0d0
    fzinv(:,:)  = 0.0d0
    fzinv0(:,:) = 0.0d0

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

! history buffers ------------------------------------------

    select case(fmode)
        case(1)
            hist_len = 2
        case(2)
            hist_len = 4
        case(3)
            call abf_init_filter_sg
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_init_arrays!')
    end select

    allocate(                                           &
            cvhist(NumOfABFCVs,hist_len),               &
            xihist(hist_len,NumOfABFCVs),               &
            micfhist(hist_len,NumOfABFCVs),             &
            epothist(hist_len),                         &
            ersthist(hist_len),                         &
            ekinhist(hist_len),                         &
            zinvhist(hist_len,NumOfABFCVs,NumOfABFCVs), &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for buffers used in ABF calculation!')
    end if

    cvhist(:,:)     = 0.0d0
    xihist(:,:)     = 0.0d0
    micfhist(:,:)   = 0.0d0
    epothist(:)     = 0.0d0
    ersthist(:)     = 0.0d0
    ekinhist(:)     = 0.0d0
    zinvhist(:,:,:) = 0.0d0

! other setup ----------------------------------------------

    ! init accumulator
    call abf_accu_init

    if( feimode .eq. 2 ) then
        call abf_init_snb_list
    end if

    if( (feimode .eq. 3) .and. (NumOfABFCVs .gt. 1) ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] feimode == 3 can be used only with one CV!')
    end if

end subroutine abf_init_arrays

!===============================================================================
! Subroutine:  abf_init_arrays
!===============================================================================

subroutine abf_init_snb_list

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer         :: i,j,k,dcv,idx,alloc_failed,fac
    real(PMFDP)     :: dx,u2
    ! --------------------------------------------------------------------------

    ! use bigger distance buffer, fac is square of this buffer
    fac = 2**2 ! 2^2 = 4

    max_snb_size = 1
    do i=1,NumOfABFCVs
        dcv = fac*ceiling(ABFCVList(i)%wfac) + 1
        max_snb_size = max_snb_size * dcv
    end do

    allocate(                                                       &
            snb_list(max_snb_size,abfaccu%PMFAccuType%tot_nbins),   &
            sweights(max_snb_size),                                 &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays for kernel smoothing in abf_init_snb_list!')
    end if

    sweights(:)   = 0.0d0
    snb_list(:,:) = 0

    do i=1,abfaccu%PMFAccuType%tot_nbins
        idx = 0
        do j=1,abfaccu%PMFAccuType%tot_nbins
            u2 = 0.0d0
            do k=1,abfaccu%PMFAccuType%tot_cvs
                dx = (abfaccu%binpos(k,i) - abfaccu%binpos(k,j)) / (ABFCVList(k)%wfac * abfaccu%PMFAccuType%sizes(k)%bin_width)
                u2 = u2 + dx**2
            end do
            if( u2 .le. real(fac,PMFDP) ) then
                idx = idx + 1
                if( idx .gt. max_snb_size ) then
                    call pmf_utils_exit(PMF_OUT,1, &
                                '[ABF] Max index into snb_list overflow in abf_init_snb_list!')
                end if
                snb_list(idx,i) = j
            end if
        end do
    end do

end subroutine abf_init_snb_list

!===============================================================================
! Subroutine:  abf_init_filter_sg
!===============================================================================

subroutine abf_init_filter_sg

    use pmf_utils
    use pmf_dat
    use abf_dat

    implicit none
    integer                     :: m, np, nr, nl, ipj, k, imj, mm, info, ld, j, kk
    integer                     :: alloc_failed
    integer, allocatable        :: lindx(:)
    real(PMFDP)                 :: s, fac
    real(PMFDP),allocatable     :: a(:,:), b(:)
    ! --------------------------------------------------------------------------

    m = fsgorder
    if( m .le. 1 ) then
        call pmf_utils_exit(PMF_OUT,1, '[ABF] fsgorder too low in abf_init_sg!')
    end if

    np = fsgframelen

    if( mod(np,2) .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1, '[ABF] fsgframelen must be an odd number in abf_init_sg!')
    end if

    if( np .le. 2 ) then
        call pmf_utils_exit(PMF_OUT,1, '[ABF] fsgframelen too low in abf_init_sg!')
    end if

    nr = (np-1)/2
    nl = nr

    allocate( a(m+1,m+1), b(m+1), lindx(m+1), sg_c0(np), sg_c1(np), sg_c2(np), stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, '[ABF] Unable to allocate memory in abf_init_sg!')
    end if

    a(:,:) = 0.0d0

    do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit.
        s = 0.0d0
        if( ipj .eq. 0 ) s = 1.0d0

        do k=1,nr
            s = s + real(k,PMFDP)**ipj
        end do

        do k=1,nl
            s = s + real(-k,PMFDP)**ipj
        end do

        mm = min(ipj,2*m-ipj)
        do imj=-mm,mm,2
            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = s
        end do
    end do

    call dgetrf(m+1,m+1,a,m+1,lindx,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_init_sg!')
    end if

    sg_c0(:) = 0.0d0
    sg_c1(:) = 0.0d0
    sg_c2(:) = 0.0d0

    do ld=0,2
        do j=1,m+1
            b(j) = 0.0d0
        end do
        b(ld+1) = 1.0d0      !Right-hand side vector is unit vector, depending on which derivative we want.

        call dgetrs('N',m+1,1,a,m+1,lindx,b,m+1,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_init_sg!')
        end if

        do k=-nl,nr                        ! Each Savitzky-Golay coefficient is the dot product
            s   = b(1)                     ! of powers of an integer with the inverse matrix row.
            fac = 1.0d0
            do mm=1,m
                fac = fac*k
                s   = s + b(mm+1)*fac
            end do
            kk = k+nl+1
            select case(ld)
                case(0)
                    sg_c0(kk) = s
                case(1)
                    sg_c1(kk) = s
                case(2)
                    sg_c2(kk) = s * 2.0d0
            end select

        end do
    end do

    invfdtx = 1.0d0 / (fdt * PMF_DT2VDT)

    sg_c1(:) = sg_c1(:) * invfdtx
    sg_c2(:) = sg_c2(:) * invfdtx * invfdtx

    ! define history length
    hist_len = np + 1 ! +1 for delayed velocity

end subroutine abf_init_filter_sg

!===============================================================================

end module abf_init
