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

    fmode           = 0         ! 0 - disable ABF, 1 - enabled ABF
    fsample         = 5000      ! output sample pariod in steps
    frestart        = .false.   ! 1 - restart job with previous data, 0 - otherwise not
    frstupdate      = 5000      ! how often is restart file written
    feimode         = 3         ! extrapolation / interpolation mode
                                ! 1 - linear ramp I
                                ! 2 - linear ramp II
                                ! 3 - block averages
    ftrjsample      = 0         ! how often save accumulator to "accumulator evolution"
    fmask_mode      = 0         ! 0 - disable ABF mask, 1 - enable ABF mask
    fapply_abf      = .true.    ! on - apply ABF, off - do not apply ABF
    fprint_icf      = .false.   ! T - print instanteous collective forces, F - do not print
    fcache_icf      = .true.    ! T - cache icf into memory and dump them at the end,
                                ! F - write icf immediatelly at each time step
    frawicf         = .true.    ! T - use raw icf data (in internal units), F - transform them to user req. units

    faccupotene     = .false.   ! accumulate PotEne

    fhramp          = 100
    fhramp_min      = 200       ! definition of linear ramp
    fhramp_max      = 500       ! definition of linear ramp
    fblock_size     = 200

    NumOfABFCVs     = 0

    fserver_enabled = .false.
    fserverkey      = ''
    fserver         = ''  ! extracted from serverkey
    fserverupdate   = 500
    fconrepeats     = 0
    fabortonmwaerr  = .true.

    client_id       = -1
    failure_counter = 0

    insidesamples   = 0
    outsidesamples  = 0

    fdtx            = 0.0d0     ! time step in internal units

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
    write(PMF_OUT,120)  '      |-> Standard ABF algorithm'
    case(2)
    write(PMF_OUT,120)  '      |-> Numerical ABF algorithm'
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown fmode in abf_init_print_header!')
    end select
    write(PMF_OUT,125)  ' Coordinate definition file (fabfdef)    : ', trim(fabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfABFCVs
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' ABF mask mode (fmask_mode)              : ', fmask_mode
    write(PMF_OUT,125)  ' ABF mask file (fabfmask)                : ', trim(fabfmask)
    write(PMF_OUT,125)  ' Apply ABF force (fapply_abf)            : ', prmfile_onoff(fapply_abf)

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Interpolation/Extrapolation '
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Extra/interpolation mode (feimode)      : ', feimode
    select case(feimode)
    case(1)
    write(PMF_OUT,120)  '      |-> Simple linear ramp'
    write(PMF_OUT,130)  ' Linear ramp size (fhramp)               : ', fhramp
    case(2)
    write(PMF_OUT,120)  '      |-> Min/Max linear ramp'
    write(PMF_OUT,130)  ' Min of accu samples in bin (fhramp_min) : ', fhramp_min
    write(PMF_OUT,130)  ' Max of accu samples in bin (fhramp_max) : ', fhramp_max
    case(3)
    write(PMF_OUT,120)  '      |-> Block averages'
    write(PMF_OUT,130)  ' Block size (fblock_size)                : ', fblock_size
    case default
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode in abf_init_print_header!')
    end select

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fabfrst)                  : ', trim(fabfrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Final restart update (frstupdate)       : ', frstupdate
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fabfout)                   : ', trim(fabfout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,125)  ' Print instantaneous forces (fprint_icf) : ', prmfile_onoff(fprint_icf)
    write(PMF_OUT,125)  ' Cache ICF (fcache_icf)                  : ', prmfile_onoff(fcache_icf)
    write(PMF_OUT,125)  ' Print ICF in internal units (frawicf)   : ', prmfile_onoff(frawicf)
    write(PMF_OUT,125)  ' Accumulate PotEne (faccupotene)         : ', prmfile_onoff(faccupotene)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample
    write(PMF_OUT,125)  ' Trajectory file (fabftrj)               : ', trim(fabftrj)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF server options:'
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
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of ABF coordinates'
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

140 format(' == Collective variable #',I4.4)

end subroutine abf_init_print_header

!===============================================================================
! Subroutine:  abf_init_arrays
!===============================================================================

subroutine abf_init_arrays

    use pmf_utils
    use pmf_dat
    use abf_dat
    use abf_accumulator

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

    fdtx = fdt*PMF_DT2VDT

    ! general arrays --------------------------------
    allocate(a0(3,NumOfLAtoms),             &
          a1(3,NumOfLAtoms),                &
          pxi0(NumOfABFCVs),                &
          pxi1(NumOfABFCVs),                &
          pxip(NumOfABFCVs),                &
          pxim(NumOfABFCVs),                &
          avg_values(NumOfABFCVs),          &
          la(NumOfABFCVs),                  &
          fz(NumOfABFCVs,NumOfABFCVs),      &
          fzinv(NumOfABFCVs,NumOfABFCVs),   &
          zd0(3,NumOfLAtoms,NumOfABFCVs),   &
          cvaluehist0(NumOfABFCVs),         &
          cvaluehist1(NumOfABFCVs),         &
          cvaluehist2(NumOfABFCVs),         &
          cvaluehist3(NumOfABFCVs),         &
          stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory for arrays used in ABF calculation!')
    end if

    a0(:,:) = 0.0d0
    a1(:,:) = 0.0d0
    pxi0(:) = 0.0d0
    pxi1(:) = 0.0d0
    pxip(:) = 0.0d0
    pxim(:) = 0.0d0
    avg_values(:) = 0.0d0
    la(:) = 0.0d0
    fz(:,:) = 0.0d0
    fzinv(:,:) = 0.0d0
    zd0(:,:,:) = 0.0d0
    cvaluehist0(:) = 0.0d0
    cvaluehist1(:) = 0.0d0
    cvaluehist2(:) = 0.0d0
    cvaluehist3(:) = 0.0d0

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
    call abf_accumulator_init

end subroutine abf_init_arrays

!===============================================================================

end module abf_init
