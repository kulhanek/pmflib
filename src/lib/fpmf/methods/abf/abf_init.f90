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

    call abf_init_print_header
    call abf_init_arrays
    call abf_init_print_summary
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

    NumOfABFCVs         = 0
    NumOfABFCVs   = 0

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
    fswitch2zero    = .false.

    gpr_len         = 1001
    gpr_boundary    = 50

    gpr_cvs_kernel  = 1
    gpr_cvs_width   = 30.0
    gpr_cvs_noise   = 0.0

    gpr_icf_cdf     = .true.
    gpr_icf_kernel  = 1
    gpr_icf_width   = 30.0
    gpr_icf_noise   = 0.0

    gpr_ene_smooth  = 0
    gpr_ene_kernel  = 1
    gpr_ene_width   = 30.0
    gpr_ene_noise   = 0.0

    gpr_rank        = -1
    gpr_rcond       = 1e-16

end subroutine abf_init_dat

!===============================================================================
! Subroutine:  abf_print_header
!===============================================================================

subroutine abf_init_print_header

    use pmf_constants
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' *********************** ADAPTIVE BIASING FORCE METHOD ************************ '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)

120 format(A)

end subroutine abf_init_print_header

!===============================================================================
! Subroutine:  abf_init_print_summary
!===============================================================================

subroutine abf_init_print_summary

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
    write(PMF_OUT,120)  ' ABF Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' ABF mode (fmode)                        : ', fmode
    select case(fmode)

    case(1)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (V)'

    case(2)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (GPR)'
    write(PMF_OUT,130)  '          gpr_len                        : ', gpr_len
    write(PMF_OUT,130)  '          gpr_boundary                   : ', gpr_boundary

    write(PMF_OUT,120)  '          === CVS GPR'
    write(PMF_OUT,130)  '              gpr_cvs_kernel             : ', gpr_cvs_kernel
    select case(gpr_cvs_kernel)
    case(0)
    write(PMF_OUT,120)  '              \-> EXP (exponential kernel)'
    case(1)
    write(PMF_OUT,120)  '              \-> MC(3/2) (Matern kernel v=3/2)'
    case(2)
    write(PMF_OUT,120)  '              \-> MC(5/2) (Matern kernel v=5/2)'
    case(3)
    write(PMF_OUT,120)  '              \-> SE (squared exponential kernel)'
    case(4)
    write(PMF_OUT,120)  '              \-> Epanechnikov (parabolic) kernel'
    case(5)
    write(PMF_OUT,120)  '              \-> Quartic (biweight) kernel'
    case(6)
    write(PMF_OUT,120)  '              \-> Triweigh kernel'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_cvs_kernel in abf_init_print_summary!')
    end select
    write(PMF_OUT,152)  '              gpr_cvs_width              : ', gpr_cvs_width,  &
                                       '['//trim(pmf_unit_label(TimeUnit))//']'
    write(PMF_OUT,160)  '              gpr_cvs_noise              : ', gpr_cvs_noise

    write(PMF_OUT,120)  '          === ICF GPR'
    write(PMF_OUT,125)  '              gpr_icf_cdf                : ', prmfile_onoff(gpr_icf_cdf)
    write(PMF_OUT,130)  '              gpr_icf_kernel             : ', gpr_icf_kernel
    select case(gpr_icf_kernel)
    case(0)
    write(PMF_OUT,120)  '              \-> EXP (exponential kernel)'
    case(1)
    write(PMF_OUT,120)  '              \-> MC(3/2) (Matern kernel v=3/2)'
    case(2)
    write(PMF_OUT,120)  '              \-> MC(5/2) (Matern kernel v=5/2)'
    case(3)
    write(PMF_OUT,120)  '              \-> SE (squared exponential kernel)'
    case(4)
    write(PMF_OUT,120)  '              \-> Epanechnikov (parabolic) kernel'
    case(5)
    write(PMF_OUT,120)  '              \-> Quartic (biweight) kernel'
    case(6)
    write(PMF_OUT,120)  '              \-> Triweigh kernel'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_icf_kernel in abf_init_print_summary!')
    end select
    write(PMF_OUT,152)  '              gpr_icf_width              : ', gpr_icf_width,  &
                                       '['//trim(pmf_unit_label(TimeUnit))//']'
    write(PMF_OUT,160)  '              gpr_icf_noise              : ', gpr_icf_noise

    write(PMF_OUT,120)  '          === ENE GPR'
    write(PMF_OUT,130)  '              gpr_ene_smooth             : ', gpr_ene_smooth
    select case(gpr_ene_smooth)
    case(0)
    write(PMF_OUT,120)  '              \-> no energy smoothing'
    case(1)
    write(PMF_OUT,120)  '              \-> smooth Etot'
    case(2)
    write(PMF_OUT,120)  '              \-> smooth Ekin'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_ene_smooth in abf_init_print_summary!')
    end select
    write(PMF_OUT,130)  '              gpr_ene_kernel             : ', gpr_ene_kernel
    select case(gpr_ene_kernel)
    case(0)
    write(PMF_OUT,120)  '              \-> EXP (exponential kernel)'
    case(1)
    write(PMF_OUT,120)  '              \-> MC(3/2) (Matern kernel v=3/2)'
    case(2)
    write(PMF_OUT,120)  '              \-> MC(5/2) (Matern kernel v=5/2)'
    case(3)
    write(PMF_OUT,120)  '              \-> SE (squared exponential kernel)'
    case(4)
    write(PMF_OUT,120)  '              \-> Epanechnikov (parabolic) kernel'
    case(5)
    write(PMF_OUT,120)  '              \-> Quartic (biweight) kernel'
    case(6)
    write(PMF_OUT,120)  '              \-> Triweigh kernel'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_ene_kernel in abf_init_print_summary!')
    end select
    write(PMF_OUT,152)  '              gpr_ene_width              : ', gpr_ene_width,  &
                                       '['//trim(pmf_unit_label(TimeUnit))//']'
    write(PMF_OUT,160)  '              gpr_ene_noise              : ', gpr_ene_noise

    write(PMF_OUT,120)  '          === K+Sigma inversion'
    write(PMF_OUT,130)  '              gpr_rank                   : ', gpr_rank
    write(PMF_OUT,160)  '              gpr_rcond                  : ', gpr_rcond

    case(10)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (F+V)'

    case(11)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (X)'

    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown fmode in abf_init_print_summary!')
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
    write(PMF_OUT,125)  ' Switch ICF to zero (fswitch2zero)       : ', prmfile_onoff(fswitch2zero)
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
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown kernel in abf_init_print_summary!')
    end select
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case(3)
    write(PMF_OUT,120)  '      |-> Linear interpolation'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode in abf_init_print_summary!')
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
152 format(A,F10.5,1X,A)
160 format(A,E10.4)

140 format(' == Collective variable #',I4.4)

end subroutine abf_init_print_summary

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

! general arrays --------------------------------
    allocate(                                   &
            la(NumOfABFCVs),                    &
            pxi0(NumOfABFCVs),                  &
            pxi1(NumOfABFCVs),                  &
            pxip(NumOfABFCVs),                  &
            pxim(NumOfABFCVs),                  &
            sfac(NumOfABFCVs),                  &
            cvave(NumOfABFCVs),                 &
            fz(NumOfABFCVs,NumOfABFCVs),        &
            fzinv(NumOfABFCVs,NumOfABFCVs),     &
            indx(NumOfABFCVs),                  &
            vv(NumOfABFCVs),                    &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays used in ABF calculation!')
    end if

    la(:)       = 0.0d0
    pxi0(:)     = 0.0d0
    pxi1(:)     = 0.0d0
    pxip(:)     = 0.0d0
    pxim(:)     = 0.0d0
    cvave(:)    = 0.0d0
    fz(:,:)     = 0.0d0
    fzinv(:,:)  = 0.0d0

    sfac(:)     = 1.0d0

! history buffers ------------------------------------------

    select case(fmode)
        ! standard
        case(1)
            hist_len = 3
        case(2)
            hist_len = gpr_len + 1 ! for ekin delay
            call abf_init_gpr
        ! experimental
        case(10)
            hist_len = 3
        case(11)
            hist_len = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_init_arrays!')
    end select

    allocate(                                           &
            cvhist(NumOfABFCVs,hist_len),               &
            micfhist(NumOfABFCVs,hist_len),             &
            fhist(3,NumOfLAtoms,hist_len),              &
            xhist(3,NumOfLAtoms,hist_len),              &
            vhist(3,NumOfLAtoms,hist_len),              &
            zdhist(3,NumOfLAtoms,NumOfABFCVs,hist_len), &
            fzinvhist(NumOfABFCVs,NumOfABFCVs,hist_len),&
            xvelhist(NumOfABFCVs,hist_len),             &
            xphist(NumOfABFCVs,hist_len),               &
            icfhist(NumOfABFCVs,hist_len),              &
            epothist(hist_len),                         &
            ersthist(hist_len),                         &
            ekinhist(hist_len),                         &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for buffers used in ABF calculation!')
    end if

    cvhist(:,:)     = 0.0d0
    micfhist(:,:)   = 0.0d0
    fhist(:,:,:)    = 0.0d0
    xhist(:,:,:)    = 0.0d0
    vhist(:,:,:)    = 0.0d0
    zdhist(:,:,:,:) = 0.0d0
    epothist(:)     = 0.0d0
    ersthist(:)     = 0.0d0
    ekinhist(:)     = 0.0d0
    fzinvhist(:,:,:)= 0.0d0
    xvelhist(:,:)   = 0.0d0
    xphist(:,:)     = 0.0d0
    icfhist(:,:)    = 0.0d0

! other setup ----------------------------------------------

! init accumulator
    call abf_accu_init

! sanity checks
    if( feimode .eq. 2 ) then
        call abf_init_snb_list
    end if

    if( (feimode .eq. 3) .and. (NumOfABFCVs .gt. 1) ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] feimode == 3 can be used only with one CV!')
    end if

end subroutine abf_init_arrays

!===============================================================================
! Function:  abf_init_gpr_kernel
!===============================================================================

real(PMFDP) function abf_init_gpr_kernel(ktype,t1,t2,width)

    use pmf_utils
    use pmf_dat
    use abf_dat

    implicit none
    integer     :: ktype,t1,t2
    real(PMFDP) :: width
    real(PMFDP) :: r,wfac
    ! --------------------------------------------------------------------------

    ! convert width to time steps
    wfac = width / fdt

    if( wfac .eq. 0.0d0 ) then
        abf_init_gpr_kernel = 0.0d0
        return
    end if

    r = abs(real(t1-t2,PMFDP)/wfac)

    select case(ktype)
        case(1)
            abf_init_gpr_kernel = (1.0d0 + sqrt(3.0) * r) * exp(- sqrt(3.0d0) * r)
        case(2)
            abf_init_gpr_kernel = (1.0d0 + sqrt(5.0d0) * r + 5.0d0/3.0d0 * r**2) * exp(- sqrt(5.0d0) * r)
        case(3)
            abf_init_gpr_kernel = exp(- 0.5d0 * r**2)
        case(4)
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel = (1.0d0 - r**2)
            else
                abf_init_gpr_kernel = 0.0d0
            end if
        case(5)
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel = (1.0d0 - r**2)**2
            else
                abf_init_gpr_kernel = 0.0d0
            end if
        case(6)
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel = (1.0d0 - r**2)**3
            else
                abf_init_gpr_kernel = 0.0d0
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown ktype in abf_init_gpr_kernel!')
    end select

end function abf_init_gpr_kernel

!===============================================================================
! Function:  abf_init_gpr_main_kernel_1d
! first derivative
!===============================================================================

real(PMFDP) function abf_init_gpr_kernel_1d(ktype,t1,t2,width)

    use pmf_utils
    use pmf_dat
    use abf_dat

    implicit none
    integer     :: ktype,t1,t2
    real(PMFDP) :: width
    real(PMFDP) :: r,wfac
    ! --------------------------------------------------------------------------

    ! convert width to time steps
    wfac = width / fdt

    if( wfac .eq. 0.0d0 ) then
        abf_init_gpr_kernel_1d = 0.0d0
        return
    end if

    r = abs(real(t1-t2,PMFDP)/wfac)

    select case(ktype)
        case(1)
            abf_init_gpr_kernel_1d = - 3.0d0 * r * exp(- sqrt(3.0d0) * r) / (wfac * fdtx) * sign(1.0d0,real(t1-t2,PMFDP))
        case(2)
            abf_init_gpr_kernel_1d = - 5.0d0/3.0d0 * r * (1.0d0 + sqrt(5.0d0) * r) * exp(- sqrt(5.0d0) * r ) / (wfac * fdtx) &
                       * sign(1.0d0,real(t1-t2,PMFDP))
        case(3)
            abf_init_gpr_kernel_1d = - exp(- 0.5d0 * r**2) * r / (wfac * fdtx) * sign(1.0d0,real(t1-t2,PMFDP))
        case(4)
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel_1d = -2.0d0 * r / (wfac * fdtx) * sign(1.0d0,real(t1-t2,PMFDP))
            else
                abf_init_gpr_kernel_1d = 0.0d0
            end if
        case(5)
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel_1d = -4.0d0 * r * (1.0d0 - r**2) / (wfac * fdtx) * sign(1.0d0,real(t1-t2,PMFDP))
            else
                abf_init_gpr_kernel_1d = 0.0d0
            end if
        case(6)
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel_1d = -6.0d0 * r * (1.0d0 - r**2)**2 / (wfac * fdtx) * sign(1.0d0,real(t1-t2,PMFDP))
            else
                abf_init_gpr_kernel_1d = 0.0d0
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown ktype in abf_init_gpr_kernel_1d!')
    end select

end function abf_init_gpr_kernel_1d

!===============================================================================
! Subroutine:  abf_init_gpr
!===============================================================================

subroutine abf_init_gpr

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer             :: i,j,alloc_failed
    ! --------------------------------------------------------------------------

! gpr_len must be long enough
    if( gpr_len  .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] gpr_len must be >= 3 in abf_init_gpr!')
    end if

! allocate arrays
    allocate(                                   &
            gpr_K_cvs(gpr_len,gpr_len),         &
            gpr_K_icf(gpr_len,gpr_len),         &
            gpr_K_ene(gpr_len,gpr_len),         &
            gpr_data(gpr_len),                  &
            gpr_model(gpr_len),                 &
            gpr_kff_cvs(gpr_len,gpr_len),       &
            gpr_kfd_cvs(gpr_len,gpr_len),       &
            gpr_kff_icf(gpr_len,gpr_len),       &
            gpr_kfd_icf(gpr_len,gpr_len),       &
            gpr_kff_ene(gpr_len,gpr_len),       &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for GPR arrays in abf_init_gpr!')
    end if

! init covariance matrix
    do i=1,gpr_len
        do j=1,gpr_len
            gpr_K_cvs(i,j) = abf_init_gpr_kernel(gpr_cvs_kernel,i,j,gpr_cvs_width)
            gpr_K_icf(i,j) = abf_init_gpr_kernel(gpr_icf_kernel,i,j,gpr_icf_width)
            gpr_K_ene(i,j) = abf_init_gpr_kernel(gpr_ene_kernel,i,j,gpr_ene_width)
        end do
    end do

! add noise
    do i=1,gpr_len
        gpr_K_cvs(i,i) = gpr_K_cvs(i,i) + gpr_cvs_noise
        gpr_K_icf(i,i) = gpr_K_icf(i,i) + gpr_icf_noise
        gpr_K_ene(i,i) = gpr_K_ene(i,i) + gpr_ene_noise
    end do

    if( fdebug ) then
        j = gpr_len/2
        open(unit=7812,file='gpr-kernel.cvs',status='UNKNOWN')
        open(unit=7813,file='gpr-kernel.icf',status='UNKNOWN')
        open(unit=7814,file='gpr-kernel.ene',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, gpr_K_cvs(i,j)
            write(7813,*) i, gpr_K_icf(i,j)
            write(7814,*) i, gpr_K_ene(i,j)
        end do
        close(7812)
        close(7813)
        close(7814)
    end if

    write(PMF_OUT,10)
    call abf_init_gpr_invK(gpr_K_cvs)

    write(PMF_OUT,20)
    call abf_init_gpr_invK(gpr_K_icf)

    write(PMF_OUT,30)
    call abf_init_gpr_invK(gpr_K_ene)

! init kff and kfd
    do i=1,gpr_len
        do j=1,gpr_len
            gpr_kff_cvs(j,i) = abf_init_gpr_kernel(gpr_cvs_kernel,i,j,gpr_cvs_width)
            gpr_kfd_cvs(j,i) = abf_init_gpr_kernel_1d(gpr_cvs_kernel,i,j,gpr_cvs_width)

            gpr_kff_icf(j,i) = abf_init_gpr_kernel(gpr_icf_kernel,i,j,gpr_icf_width)
            gpr_kfd_icf(j,i) = abf_init_gpr_kernel_1d(gpr_icf_kernel,i,j,gpr_icf_width)

            gpr_kff_ene(j,i) = abf_init_gpr_kernel(gpr_ene_kernel,i,j,gpr_ene_width)
        end do
    end do

 10 format('>>> ABF GPR - CVS: (K+Sigma)^-1')
 20 format('>>> ABF GPR - ICF: (K+Sigma)^-1')
 30 format('>>> ABF GPR - ENE: (K+Sigma)^-1')

end subroutine abf_init_gpr

!===============================================================================
! Subroutine:  abf_init_gpr_invK
!===============================================================================

subroutine abf_init_gpr_invK(mat)

    use pmf_dat
    use abf_dat
    use pmf_utils

    implicit none
    real(PMFDP)                 :: mat(:,:)
    integer                     :: i,alloc_failed,info,lwork,irank
    real(PMFDP)                 :: minv, maxv, a_rcond
    real(PMFDP),allocatable     :: sig(:)
    real(PMFDP),allocatable     :: u(:,:)
    real(PMFDP),allocatable     :: vt(:,:)
    real(PMFDP),allocatable     :: sig_plus(:,:)
    real(PMFDP),allocatable     :: temp_mat(:,:)
    real(PMFDP),allocatable     :: twork(:)
    integer,allocatable         :: iwork(:)
    ! --------------------------------------------------------------------------

! allocate arrays
    allocate(                           &
            sig(gpr_len),               &
            u(gpr_len,gpr_len),         &
            vt(gpr_len,gpr_len),        &
            sig_plus(gpr_len,gpr_len),  &
            temp_mat(gpr_len,gpr_len),  &
            iwork(8*gpr_len),           &
            twork(1),                   &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory I for GPR arrays in abf_init_gpr_invK!')
    end if

    sig(:)          = 0.0d0
    u(:,:)          = 0.0d0
    vt(:,:)         = 0.0d0
    sig_plus(:,:)   = 0.0d0
    temp_mat(:,:)   = 0.0d0

! query work size
    lwork = -1
    call dgesdd('A', gpr_len, gpr_len, mat, gpr_len, sig, u, gpr_len, vt, gpr_len, twork, lwork, iwork, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to run SVD I in abf_init_gpr_invK!')
    end if


    lwork = int(twork(1)) + 1
    if( lwork < 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to illegal work size in abf_init_gpr_invK!')
    end if

    deallocate(twork)
    allocate(   twork(lwork),   &
                stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory II for GPR arrays in abf_init_gpr_invK!')
    end if

! run SVD
    call dgesdd('A', gpr_len, gpr_len, mat, gpr_len, sig, u, gpr_len, vt, gpr_len, twork, lwork, iwork, info)

    if( info .ne. 0 )  then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to run SVD II in abf_init_gpr_invK!')
    end if

    deallocate(iwork,twork)

! invert singular numbers
    maxv = sig(1)
    minv = sig(1)
    do i=1,gpr_len
        if( maxv .lt. sig(i) ) then
            maxv = sig(i)
        end if
        if( minv .gt. sig(i) ) then
            minv = sig(i)
        end if
    end do

    a_rcond = minv/maxv

    irank = 0
    if( gpr_rank .gt. 1 ) then
        do i=1,gpr_len
            if( i .le. gpr_rank ) then
               sig_plus(i,i) = 1.0d0/sig(i)
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    else
        do i=1,gpr_len
            if( sig(i) .gt. gpr_rcond*maxv ) then
               sig_plus(i,i) = 1.0d0/sig(i)
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    end if

! build pseudoinverse: V*sig_plus*UT
    call dgemm('N', 'T', gpr_len, gpr_len, gpr_len, 1.0d0, sig_plus, gpr_len, u, gpr_len, 0.0d0, temp_mat, gpr_len)
    call dgemm('T', 'N', gpr_len, gpr_len, gpr_len, 1.0d0, vt, gpr_len, temp_mat, gpr_len, 0.0d0, mat, gpr_len)

    deallocate(sig,sig_plus,u,vt,temp_mat)

    write(PMF_OUT,10) gpr_len,irank,a_rcond

 10 format('    Size = ',I5,'; Rank = ',I5,'; Real rcond = ',E10.5)

end subroutine abf_init_gpr_invK

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

    max_snb_size = max_snb_size + 1 ! terminating null

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

end module abf_init
