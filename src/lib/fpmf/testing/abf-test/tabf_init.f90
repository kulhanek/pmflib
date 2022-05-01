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

!    gpr_len         = 1001
!    gpr_boundary    = 50
!
!    gpr_cvs_kernel  = 1
!    gpr_cvs_width   = 30.0
!    gpr_cvs_noise   = 0.0
!
!    gpr_icf_cdf     = .true.
!    gpr_icf_kernel  = 1
!    gpr_icf_width   = 30.0
!    gpr_icf_noise   = 0.0
!
!    gpr_ene_smooth  = 0
!    gpr_ene_kernel  = 1
!    gpr_ene_width   = 30.0
!    gpr_ene_noise   = 0.0
!
!    gpr_rank        = -1
!    gpr_rcond       = 1e-16
!
!    fsgframelen     = 11
!    fsgorder        = 5

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

!!===============================================================================
!! Subroutine:  abf_init_gpr
!!===============================================================================
!
!subroutine abf_init_gpr
!
!    use pmf_utils
!    use pmf_dat
!    use abf_dat
!    use abf_accu
!
!    implicit none
!    integer             :: i,j,alloc_failed
!    ! --------------------------------------------------------------------------
!
!! gpr_len must be long enough
!    if( gpr_len  .le. 3 ) then
!        call pmf_utils_exit(PMF_OUT,1,'[ABF] gpr_len must be >= 3 in abf_init_gpr!')
!    end if
!
!! allocate arrays
!    allocate(                                   &
!            gpr_K_cvs(gpr_len,gpr_len),         &
!            gpr_K_icf(gpr_len,gpr_len),         &
!            gpr_K_ene(gpr_len,gpr_len),         &
!            gpr_data(gpr_len),                  &
!            gpr_model(gpr_len),                 &
!            gpr_kff_cvs(gpr_len,gpr_len),       &
!            gpr_kfd_cvs(gpr_len,gpr_len),       &
!            gpr_kff_icf(gpr_len,gpr_len),       &
!            gpr_kfd_icf(gpr_len,gpr_len),       &
!            gpr_kff_ene(gpr_len,gpr_len),       &
!            stat= alloc_failed )
!
!    if( alloc_failed .ne. 0 ) then
!        call pmf_utils_exit(PMF_OUT,1, &
!            '[ABF] Unable to allocate memory for GPR arrays in abf_init_gpr!')
!    end if
!
!! init covariance matrix
!    do i=1,gpr_len
!        do j=1,gpr_len
!            gpr_K_cvs(i,j) = abf_init_gpr_kernel(gpr_cvs_kernel,i,j,gpr_cvs_width)
!            gpr_K_icf(i,j) = abf_init_gpr_kernel(gpr_icf_kernel,i,j,gpr_icf_width)
!            gpr_K_ene(i,j) = abf_init_gpr_kernel(gpr_ene_kernel,i,j,gpr_ene_width)
!        end do
!    end do
!
!! add noise
!    do i=1,gpr_len
!        gpr_K_cvs(i,i) = gpr_K_cvs(i,i) + gpr_cvs_noise
!        gpr_K_icf(i,i) = gpr_K_icf(i,i) + gpr_icf_noise
!        gpr_K_ene(i,i) = gpr_K_ene(i,i) + gpr_ene_noise
!    end do
!
!    if( fdebug ) then
!        j = gpr_len/2
!        open(unit=7812,file='gpr-kernel.cvs',status='UNKNOWN')
!        open(unit=7813,file='gpr-kernel.icf',status='UNKNOWN')
!        open(unit=7814,file='gpr-kernel.ene',status='UNKNOWN')
!        do i=1,gpr_len
!            write(7812,*) i, gpr_K_cvs(i,j)
!            write(7813,*) i, gpr_K_icf(i,j)
!            write(7814,*) i, gpr_K_ene(i,j)
!        end do
!        close(7812)
!        close(7813)
!        close(7814)
!    end if
!
!    write(PMF_OUT,10)
!    call abf_init_gpr_invK(gpr_K_cvs)
!
!    write(PMF_OUT,20)
!    call abf_init_gpr_invK(gpr_K_icf)
!
!    write(PMF_OUT,30)
!    call abf_init_gpr_invK(gpr_K_ene)
!
!! init kff and kfd
!    do i=1,gpr_len
!        do j=1,gpr_len
!            gpr_kff_cvs(j,i) = abf_init_gpr_kernel(gpr_cvs_kernel,i,j,gpr_cvs_width)
!            gpr_kfd_cvs(j,i) = abf_init_gpr_kernel_1d(gpr_cvs_kernel,i,j,gpr_cvs_width)
!
!            gpr_kff_icf(j,i) = abf_init_gpr_kernel(gpr_icf_kernel,i,j,gpr_icf_width)
!            gpr_kfd_icf(j,i) = abf_init_gpr_kernel_1d(gpr_icf_kernel,i,j,gpr_icf_width)
!
!            gpr_kff_ene(j,i) = abf_init_gpr_kernel(gpr_ene_kernel,i,j,gpr_ene_width)
!        end do
!    end do
!
! 10 format('>>> ABF GPR - CVS: (K+Sigma)^-1')
! 20 format('>>> ABF GPR - ICF: (K+Sigma)^-1')
! 30 format('>>> ABF GPR - ENE: (K+Sigma)^-1')
!
!end subroutine abf_init_gpr
!
!!===============================================================================
!! Subroutine:  abf_init_gpr_invK
!!===============================================================================
!
!subroutine abf_init_gpr_invK(mat)
!
!    use pmf_dat
!    use abf_dat
!    use pmf_utils
!
!    implicit none
!    real(PMFDP)                 :: mat(:,:)
!    integer                     :: i,alloc_failed,info,lwork,irank
!    real(PMFDP)                 :: minv, maxv, a_rcond
!    real(PMFDP),allocatable     :: sig(:)
!    real(PMFDP),allocatable     :: u(:,:)
!    real(PMFDP),allocatable     :: vt(:,:)
!    real(PMFDP),allocatable     :: sig_plus(:,:)
!    real(PMFDP),allocatable     :: temp_mat(:,:)
!    real(PMFDP),allocatable     :: twork(:)
!    integer,allocatable         :: iwork(:)
!    ! --------------------------------------------------------------------------
!
!! allocate arrays
!    allocate(                           &
!            sig(gpr_len),               &
!            u(gpr_len,gpr_len),         &
!            vt(gpr_len,gpr_len),        &
!            sig_plus(gpr_len,gpr_len),  &
!            temp_mat(gpr_len,gpr_len),  &
!            iwork(8*gpr_len),           &
!            twork(1),                   &
!            stat= alloc_failed )
!
!    if( alloc_failed .ne. 0 ) then
!        call pmf_utils_exit(PMF_OUT,1, &
!            '[ABF] Unable to allocate memory I for GPR arrays in abf_init_gpr_invK!')
!    end if
!
!    sig(:)          = 0.0d0
!    u(:,:)          = 0.0d0
!    vt(:,:)         = 0.0d0
!    sig_plus(:,:)   = 0.0d0
!    temp_mat(:,:)   = 0.0d0
!
!! query work size
!    lwork = -1
!    call dgesdd('A', gpr_len, gpr_len, mat, gpr_len, sig, u, gpr_len, vt, gpr_len, twork, lwork, iwork, info)
!
!    if( info .ne. 0 ) then
!        call pmf_utils_exit(PMF_OUT,1, &
!            '[ABF] Unable to run SVD I in abf_init_gpr_invK!')
!    end if
!
!
!    lwork = int(twork(1)) + 1
!    if( lwork < 0 ) then
!        call pmf_utils_exit(PMF_OUT,1, &
!            '[ABF] Unable to illegal work size in abf_init_gpr_invK!')
!    end if
!
!    deallocate(twork)
!    allocate(   twork(lwork),   &
!                stat= alloc_failed )
!
!    if( alloc_failed .ne. 0 ) then
!        call pmf_utils_exit(PMF_OUT,1, &
!            '[ABF] Unable to allocate memory II for GPR arrays in abf_init_gpr_invK!')
!    end if
!
!! run SVD
!    call dgesdd('A', gpr_len, gpr_len, mat, gpr_len, sig, u, gpr_len, vt, gpr_len, twork, lwork, iwork, info)
!
!    if( info .ne. 0 )  then
!        call pmf_utils_exit(PMF_OUT,1, &
!            '[ABF] Unable to run SVD II in abf_init_gpr_invK!')
!    end if
!
!    deallocate(iwork,twork)
!
!! invert singular numbers
!    maxv = sig(1)
!    minv = sig(1)
!    do i=1,gpr_len
!        if( maxv .lt. sig(i) ) then
!            maxv = sig(i)
!        end if
!        if( minv .gt. sig(i) ) then
!            minv = sig(i)
!        end if
!    end do
!
!    a_rcond = minv/maxv
!
!    irank = 0
!    if( gpr_rank .gt. 1 ) then
!        do i=1,gpr_len
!            if( i .le. gpr_rank ) then
!               sig_plus(i,i) = 1.0d0/sig(i)
!               irank = irank + 1
!            else
!               sig_plus(i,i) = 0.0d0
!            end if
!        end do
!    else
!        do i=1,gpr_len
!            if( sig(i) .gt. gpr_rcond*maxv ) then
!               sig_plus(i,i) = 1.0d0/sig(i)
!               irank = irank + 1
!            else
!               sig_plus(i,i) = 0.0d0
!            end if
!        end do
!    end if
!
!! build pseudoinverse: V*sig_plus*UT
!    call dgemm('N', 'T', gpr_len, gpr_len, gpr_len, 1.0d0, sig_plus, gpr_len, u, gpr_len, 0.0d0, temp_mat, gpr_len)
!    call dgemm('T', 'N', gpr_len, gpr_len, gpr_len, 1.0d0, vt, gpr_len, temp_mat, gpr_len, 0.0d0, mat, gpr_len)
!
!    deallocate(sig,sig_plus,u,vt,temp_mat)
!
!    write(PMF_OUT,10) gpr_len,irank,a_rcond
!
! 10 format('    Size = ',I5,'; Rank = ',I5,'; Real rcond = ',E10.5)
!
!end subroutine abf_init_gpr_invK
!
!!===============================================================================
!! Subroutine:  abf_init_filter_sg
!!===============================================================================
!
!subroutine abf_init_filter_sg
!
!    use pmf_utils
!    use pmf_dat
!    use abf_dat
!
!    implicit none
!    integer                     :: m, np, nr, nl, ipj, k, imj, mm, info, ld, j, kk
!    integer                     :: alloc_failed
!    integer, allocatable        :: lindx(:)
!    real(PMFDP)                 :: s, fac
!    real(PMFDP),allocatable     :: a(:,:), b(:)
!    ! --------------------------------------------------------------------------
!
!    m = fsgorder
!    if( m .le. 1 ) then
!        call pmf_utils_exit(PMF_OUT,1, '[ABF] fsgorder too low in abf_init_sg!')
!    end if
!
!    np = fsgframelen
!
!    if( mod(np,2) .ne. 1 ) then
!        call pmf_utils_exit(PMF_OUT,1, '[ABF] fsgframelen must be an odd number in abf_init_sg!')
!    end if
!
!    if( np .le. 2 ) then
!        call pmf_utils_exit(PMF_OUT,1, '[ABF] fsgframelen too low in abf_init_sg!')
!    end if
!
!    nr = (np-1)/2
!    nl = nr
!
!    allocate( a(m+1,m+1), b(m+1), lindx(m+1), sg_c0(np), sg_c1(np), sg_c2(np), stat= alloc_failed )
!
!    if( alloc_failed .ne. 0 ) then
!        call pmf_utils_exit(PMF_OUT,1, '[ABF] Unable to allocate memory in abf_init_sg!')
!    end if
!
!    a(:,:) = 0.0d0
!
!    do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit.
!        s = 0.0d0
!        if( ipj .eq. 0 ) s = 1.0d0
!
!        do k=1,nr
!            s = s + real(k,PMFDP)**ipj
!        end do
!
!        do k=1,nl
!            s = s + real(-k,PMFDP)**ipj
!        end do
!
!        mm = min(ipj,2*m-ipj)
!        do imj=-mm,mm,2
!            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = s
!        end do
!    end do
!
!    call dgetrf(m+1,m+1,a,m+1,lindx,info)
!    if( info .ne. 0 ) then
!        call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_init_sg!')
!    end if
!
!    sg_c0(:) = 0.0d0
!    sg_c1(:) = 0.0d0
!    sg_c2(:) = 0.0d0
!
!    do ld=0,2
!        do j=1,m+1
!            b(j) = 0.0d0
!        end do
!        b(ld+1) = 1.0d0      !Right-hand side vector is unit vector, depending on which derivative we want.
!
!        call dgetrs('N',m+1,1,a,m+1,lindx,b,m+1,info)
!        if( info .ne. 0 ) then
!            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_init_sg!')
!        end if
!
!        do k=-nl,nr                        ! Each Savitzky-Golay coefficient is the dot product
!            s   = b(1)                     ! of powers of an integer with the inverse matrix row.
!            fac = 1.0d0
!            do mm=1,m
!                fac = fac*k
!                s   = s + b(mm+1)*fac
!            end do
!            kk = k+nl+1
!            select case(ld)
!                case(0)
!                    sg_c0(kk) = s
!                case(1)
!                    sg_c1(kk) = s
!                case(2)
!                    sg_c2(kk) = s * 2.0d0
!            end select
!
!        end do
!    end do
!
!    sg_c1(:) = sg_c1(:) * ifdtx
!    sg_c2(:) = sg_c2(:) * ifdtx * ifdtx
!
!end subroutine abf_init_filter_sg

!===============================================================================

end module tabf_init
