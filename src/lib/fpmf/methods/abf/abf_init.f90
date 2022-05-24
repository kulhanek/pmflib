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
    fentdecomp      = .false.
    frecord         = .false.

    ftds_ekin_src   = 1
    ftds_add_bias   = .false.

    fepotaverage    = 0.0d0
    fekinaverage    = 0.0d0

    flpfilter       = 0

    feimode         = 1
    fhramp_min      = 20000
    fhramp_max      = 30000

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

    gpr_len         = 51

    gpr_filter_aicf = .false.
    gpr_filter_etot = .false.
    gpr_filter_ficf = .false.
    gpr_filter_epot = .false.
    gpr_filter_ekin = .false.

    gpr_kernel      = 3
    gpr_width       = 12.0
    gpr_delay       = 0d0
    gpr_noise       = 1e-6

    gpr_rank        = -1
    gpr_rcond       = 1e-16
    gpr_rsigma      = 1e-5

    cbuff_pos       = 0

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
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (3pV)'

    case(2)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (3pF)'

    case(3)
    write(PMF_OUT,120)  '      |-> Simplified ABF algorithm (3pV4)'

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
    write(PMF_OUT,130)  ' Kinetic energy source (ftds_ekin_src)   : ', ftds_ekin_src
    select case(ftds_ekin_src)
    case(1)
    write(PMF_OUT,120)  '      |-> VV (velocity Verlet)'
    case(2)
    write(PMF_OUT,120)  '      |-> LF (leap-frog)'
    case(3)
    write(PMF_OUT,120)  '      |-> PE (potential energy corrected)'
    case(4)
    write(PMF_OUT,120)  '      |-> LF (corrected)'
    case(5)
    write(PMF_OUT,120)  '      |-> KIE12'
    case(6)
    write(PMF_OUT,120)  '      |-> KIE11'
    case(7)
    write(PMF_OUT,120)  '      |-> EkinH'
    case(8)
    write(PMF_OUT,120)  '      |-> LF (corrected, II)'
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown kinetic energy source in abf_init_print_summary!')
    end select
    write(PMF_OUT,125)  ' Incl. ABF bias into -TdS (ftds_add_bias): ', prmfile_onoff(ftds_add_bias)


! LP filter
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Low-pass filters:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    select case(flpfilter)
    case(0)
    write(PMF_OUT,120)  '      |-> No filter'
    case(1)
    write(PMF_OUT,120)  '      |-> GPR (Gaussian Process Regression)'
    write(PMF_OUT,130)  '          gpr_len                        : ', gpr_len
    write(PMF_OUT,125)  '          gpr_filter_aicf                : ', prmfile_onoff(gpr_filter_aicf)
    write(PMF_OUT,125)  '          gpr_filter_etot                : ', prmfile_onoff(gpr_filter_etot)
    write(PMF_OUT,125)  '          gpr_filter_ficf                : ', prmfile_onoff(gpr_filter_ficf)
    write(PMF_OUT,125)  '          gpr_filter_epot                : ', prmfile_onoff(gpr_filter_epot)
    write(PMF_OUT,125)  '          gpr_filter_ekin                : ', prmfile_onoff(gpr_filter_ekin)

    write(PMF_OUT,130)  '          gpr_kernel                     : ', gpr_kernel
    select case(gpr_kernel)
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
    case(7)
    write(PMF_OUT,120)  '              \-> sinc kernel'
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown gpr_kernel in abf_init_print_summary!')
    end select
    write(PMF_OUT,152)  '              gpr_width                  : ', gpr_width,  &
                                       '['//trim(pmf_unit_label(TimeUnit))//']'
    write(PMF_OUT,160)  '              gpr_delay                  : ', gpr_delay
    write(PMF_OUT,160)  '              gpr_noise                  : ', gpr_noise

    write(PMF_OUT,120)  '          === K+Sigma inversion'
    write(PMF_OUT,130)  '              gpr_rank                   : ', gpr_rank
    write(PMF_OUT,160)  '              gpr_rcond                  : ', gpr_rcond
    write(PMF_OUT,160)  '              gpr_rsigma                 : ', gpr_rsigma
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown flpfilter mode in abf_init_print_summary!')
    end select

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
            pxi2(NumOfABFCVs),                  &
            pxi3(NumOfABFCVs),                  &
            pxip(NumOfABFCVs),                  &
            pxif(NumOfABFCVs),                  &
            sfac(NumOfABFCVs),                  &
            cvave(NumOfABFCVs),                 &
            fz(NumOfABFCVs,NumOfABFCVs),        &
            fzinv(NumOfABFCVs,NumOfABFCVs),     &
            fzshake(NumOfABFSHAKECVs,NumOfABFSHAKECVs),     &
            indxshake(NumOfABFSHAKECVs),                    &
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
    pxi2(:)     = 0.0d0
    pxi3(:)     = 0.0d0
    pxip(:)     = 0.0d0
    pxif(:)     = 0.0d0
    cvave(:)    = 0.0d0
    fz(:,:)     = 0.0d0
    fzinv(:,:)  = 0.0d0

    fzshake(:,:)    = 0.0d0
    sfac(:)         = 1.0d0

! history buffers ------------------------------------------

    select case(fmode)
        ! standard - 3pV
        case(1)
            hist_len = 5
        ! standard - 3pF
        case(2)
            hist_len = 3
        ! standard - 3pV4
        case(3)
            hist_len = 5
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_init_arrays!')
    end select

    allocate(                                           &
            cvhist(NumOfABFCVs,hist_len),               &
            micfhist(NumOfABFCVs,hist_len),             &
            fhist(3,NumOfLAtoms,hist_len),              &
            shist(3,NumOfLAtoms,hist_len),              &
            lhist(3,NumOfLAtoms,hist_len),              &
            vhist(3,NumOfLAtoms,hist_len),              &
            zdhist(3,NumOfLAtoms,NumOfABFCVs,hist_len), &
            fzinvhist(NumOfABFCVs,NumOfABFCVs,hist_len),&
            xphist(NumOfABFCVs,hist_len),               &
            icfhist(NumOfABFCVs,hist_len),              &
            epothist(hist_len),                         &
            ersthist(hist_len),                         &
            ekinhist(hist_len),                         &
            ekinvvhist(hist_len),                       &
            ekinlfhist(hist_len),                       &
            mtchist(hist_len),                          &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for buffers used in ABF calculation!')
    end if

    cvhist(:,:)     = 0.0d0
    micfhist(:,:)   = 0.0d0
    fhist(:,:,:)    = 0.0d0
    shist(:,:,:)    = 0.0d0
    lhist(:,:,:)    = 0.0d0
    vhist(:,:,:)    = 0.0d0
    zdhist(:,:,:,:) = 0.0d0
    epothist(:)     = 0.0d0
    ersthist(:)     = 0.0d0
    ekinhist(:)     = 0.0d0
    fzinvhist(:,:,:)= 0.0d0
    xphist(:,:)     = 0.0d0
    icfhist(:,:)    = 0.0d0

! filter setup
    select case(flpfilter)
        case(0)
            ! nothing
        case(1)
            ! GPR LP filter
            lp_len = gpr_len
            allocate(      &
                cvlph(NumOfABFCVs,lp_len),            &
                flph(NumOfABFCVs,lp_len),             &
                slph(NumOfABFCVs,lp_len),             &
                vlph(NumOfABFCVs,lp_len),             &
                llph(NumOfABFCVs,lp_len),             &
                blph(NumOfABFCVs,lp_len),             &
                mtclph(lp_len),                       &
                epotlph(lp_len),                      &
                erstlph(lp_len),                      &
                ekinlph(lp_len),                      &
                stat= alloc_failed )

                if( alloc_failed .ne. 0 ) then
                    call pmf_utils_exit(PMF_OUT,1, &
                        '[ABF] Unable to allocate memory for LP buffers used in ABF calculation!')
                end if

                cvlph(:,:)  = 0.0d0
                flph(:,:)   = 0.0d0
                slph(:,:)   = 0.0d0
                vlph(:,:)   = 0.0d0
                llph(:,:)   = 0.0d0
                blph(:,:)   = 0.0d0
                mtclph(:)   = 0.0d0
                epotlph(:)  = 0.0d0
                erstlph(:)  = 0.0d0
                ekinlph(:)  = 0.0d0

                call abf_init_gpr()

        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented flpfilter mode in abf_init_arrays!')
    end select

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
    integer     :: ktype
    real(PMFDP) :: t1,t2
    real(PMFDP) :: width
    real(PMFDP) :: r,wfac
    ! --------------------------------------------------------------------------

    ! convert width to time steps
    wfac = width / fdt

    if( wfac .eq. 0.0d0 ) then
        abf_init_gpr_kernel = 0.0d0
        return
    end if

    r = abs((t1-t2)/wfac)

    select case(ktype)
        case(1)
            ! MC(3/2)
            abf_init_gpr_kernel = (1.0d0 + sqrt(3.0) * r) * exp(- sqrt(3.0d0) * r)
        case(2)
            ! MC(5/2)
            abf_init_gpr_kernel = (1.0d0 + sqrt(5.0d0) * r + 5.0d0/3.0d0 * r**2) * exp(- sqrt(5.0d0) * r)
        case(3)
            ! SE
            abf_init_gpr_kernel = exp(- 0.5d0 * r**2)
        case(4)
            ! Epanechnikov (parabolic) kernel
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel = (1.0d0 - r**2)
            else
                abf_init_gpr_kernel = 0.0d0
            end if
        case(5)
            ! Quartic (biweight) kernel
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel = (1.0d0 - r**2)**2
            else
                abf_init_gpr_kernel = 0.0d0
            end if
        case(6)
            ! Triweigh kernel
            if( r .le. 1.0d0 ) then
                abf_init_gpr_kernel = (1.0d0 - r**2)**3
            else
                abf_init_gpr_kernel = 0.0d0
            end if
        case(7)
            ! sinc(x)
            if( r .eq. 0.0d0 ) then
                abf_init_gpr_kernel = 1.0d0
            else
                abf_init_gpr_kernel = 2.0d0/wfac*sin(2.0d0*PMF_PI*r)/(2.0d0*PMF_PI*r)
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown ktype in abf_init_gpr_kernel!')
    end select

end function abf_init_gpr_kernel

!===============================================================================
! Subroutine:  abf_init_gpr
!===============================================================================

subroutine abf_init_gpr

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accu

    implicit none
    integer     :: i,j,alloc_failed,gpr_mid
    real(PMFDP) :: v
    ! --------------------------------------------------------------------------

! gpr_len must be long enough
    if( gpr_len  .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] gpr_len must be >= 3 in abf_init_gpr!')
    end if

    if( mod(gpr_len,2) .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1, '[ABF] gpr_len must be an odd number in abf_init_gpr!')
    end if

! allocate arrays
    allocate(                                   &
            gpr_K(gpr_len,gpr_len),             &
            gpr_data(gpr_len),                  &
            gpr_model(gpr_len),                 &
            gpr_kff1(gpr_len),                  &
            gpr_kff2(gpr_len),                  &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for GPR arrays in abf_init_gpr!')
    end if

! init co-variance matrix
    do i=1,gpr_len
        do j=1,gpr_len
            gpr_K(i,j) = abf_init_gpr_kernel(gpr_kernel,real(i,PMFDP),real(j,PMFDP),gpr_width)
        end do
    end do

! add noise
    do i=1,gpr_len
        gpr_K(i,i) = gpr_K(i,i) + gpr_noise
    end do

    if( fdebug ) then
        j = gpr_len/2
        open(unit=7812,file='gpr-kernel.K',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, gpr_K(i,j)
        end do
        close(7812)
    end if

    write(PMF_OUT,10)
    call abf_init_gpr_invK(gpr_K,gpr_K_logdet)

    gpr_mid = gpr_len/2 + 1

! init kff1
    do j=1,gpr_len
        gpr_data(j) = abf_init_gpr_kernel(gpr_kernel,real(gpr_mid,PMFDP),real(j,PMFDP),gpr_width)
    end do

    if( fdebug ) then
        open(unit=7812,file='gpr-kernel.kff1',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, gpr_data(i)
        end do
        close(7812)
    end if

! multiply by Kinv from right
    do i=1,gpr_len
        v = 0.0d0
        do j=1,gpr_len
            v = v + gpr_data(j) * gpr_K(j,i)
        end do
        gpr_kff1(i) = v
    end do

    if( fdebug ) then
        open(unit=7812,file='gpr-kernel.kff1-K',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, gpr_kff1(i)
        end do
        close(7812)
    end if

! -----------------------------
    do j=1,gpr_len
        gpr_data(j) = abf_init_gpr_kernel(gpr_kernel,real(gpr_mid,PMFDP)+gpr_delay,real(j,PMFDP),gpr_width)
    end do

    if( fdebug ) then
        open(unit=7812,file='gpr-kernel.kff2',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, gpr_data(i)
        end do
        close(7812)
    end if

! multiply by Kinv from right
    do i=1,gpr_len
        v = 0.0d0
        do j=1,gpr_len
            v = v + gpr_data(j) * gpr_K(j,i)
        end do
        gpr_kff2(i) = v
    end do

    if( fdebug ) then
        open(unit=7812,file='gpr-kernel.kff2-K',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, gpr_kff2(i)
        end do
        close(7812)
    end if

 10 format('>>> ABF GPR: (K+Sigma)^-1')

end subroutine abf_init_gpr

!===============================================================================
! Subroutine:  abf_init_gpr_invK
!===============================================================================

subroutine abf_init_gpr_invK(mat,logdet)

    use pmf_dat
    use abf_dat
    use pmf_utils

    implicit none
    real(PMFDP)                 :: mat(:,:)
    real(PMFDP)                 :: logdet
    integer                     :: i,alloc_failed,info,lwork,irank
    real(PMFDP)                 :: minv, maxv, a_rcond, r_sigma
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

    if( fdebug ) then
        open(unit=7812,file='gpr-kernel.isigma',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, 1.0d0/sig(i)
        end do
        close(7812)
    end if

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
    logdet = 0.0d0  ! this is logarithm of the determinant of original matrix
    r_sigma = 0.0d0
    if( gpr_rank .gt. 1 ) then
        do i=1,gpr_len
            if( i .le. gpr_rank ) then
               r_sigma = sig(i)
               sig_plus(i,i) = 1.0d0/sig(i)
               logdet        = logdet +  log(sig(i))
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    else if( gpr_rsigma .ne. 0.0d0 ) then
        do i=1,gpr_len
            if( sig(i) .gt. gpr_rsigma ) then
               r_sigma = sig(i)
               sig_plus(i,i) = 1.0d0/sig(i)
               logdet        = logdet +  log(sig(i))
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    else
        do i=1,gpr_len
            if( sig(i) .gt. gpr_rcond*maxv ) then
               r_sigma = sig(i)
               sig_plus(i,i) = 1.0d0/sig(i)
               logdet        = logdet +  log(sig(i))
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

    write(PMF_OUT,10) gpr_len,irank,r_sigma,a_rcond

 10 format('    Size = ',I5,'; Rank = ',I5,'; Rank sigma = ',E10.5,'; Real rcond = ',E10.5)

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
    integer         :: i,j,k,idx,alloc_failed
    real(PMFDP)     :: dx,u2,fac
    ! --------------------------------------------------------------------------

! use bigger distance buffer, fac is square of this buffer
    fac = 2**2 ! 2^2 = 4

! calculate the number of required pairs
    max_snb_size = 0
    do i=1,abfaccu%PMFAccuType%tot_nbins
        do j=1,abfaccu%PMFAccuType%tot_nbins
            u2 = 0.0d0
            do k=1,abfaccu%PMFAccuType%tot_cvs
                dx = abfaccu%PMFAccuType%sizes(k)%cv%get_deviation(abfaccu%binpos(k,i),abfaccu%binpos(k,j)) &
                   / (ABFCVList(k)%wfac * abfaccu%PMFAccuType%sizes(k)%bin_width)
                u2 = u2 + dx**2
            end do
            if( u2 .le. fac ) then
                max_snb_size = max_snb_size + 1
            end if
        end do
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
                dx = abfaccu%PMFAccuType%sizes(k)%cv%get_deviation(abfaccu%binpos(k,i),abfaccu%binpos(k,j)) &
                   / (ABFCVList(k)%wfac * abfaccu%PMFAccuType%sizes(k)%bin_width)
                u2 = u2 + dx**2
            end do
            if( u2 .le. fac ) then
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
