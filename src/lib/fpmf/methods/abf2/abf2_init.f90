!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf2_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf2_init_method
!===============================================================================

subroutine abf2_init_method

    use abf2_output
    use abf2_restart
    use abf2_trajectory
    use abf2_client

    implicit none
    ! --------------------------------------------------------------------------

    call abf2_init_arrays
    call abf2_init_print_header
    call abf2_output_open
    call abf2_restart_read
    call abf2_trajectory_open
    call abf2_client_register
    call abf2_client_get_initial_data
    call abf2_output_write_header

end subroutine abf2_init_method

!===============================================================================
! Subroutine:  abf2_init_dat
!===============================================================================

subroutine abf2_init_dat

    use abf2_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0         ! 0 - disable ABF, 1 - enabled ABF
    fsample         = 500       ! output sample pariod in steps
    frestart        = .false.   ! 1 - restart job with previous data, 0 - otherwise not
    frstupdate      = 1000      ! how often is restart file written
    ftrjsample      = 0         ! how often save accumulator to "accumulator evolution"
    fmask_mode      = 0         ! 0 - disable ABF mask, 1 - enable ABF mask
    fblock_size     = 100       ! number of samples for blocking pre-averaging
    fintrpl         = 0         ! ABF force interpolation: 0 - no interpolation, 1 - linear interpolation
    fapply_abf      = .true.    ! on - apply ABF, off - do not apply ABF

    NumOfABFCVs     = 0

    fserver_enabled = .false.
    fserverkey      = ''
    fserver         = ''
    fpassword       = ''
    fserverupdate   = 500
    fconrepeats     = 0
    fabortonmwaerr  = .true.

    use_key         = .false.
    client_id       = -1
    failure_counter = 0

    insidesamples   = 0
    outsidesamples  = 0

    fdtx            = 0.0d0     ! time step in internal units

end subroutine abf2_init_dat

!===============================================================================
! Subroutine:  abf2_print_header
!===============================================================================

subroutine abf2_init_print_header

    use prmfile
    use pmf_constants
    use pmf_dat
    use abf2_dat
    use abf2_cvs
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
    write(PMF_OUT,125)  ' Coordinate definition file (fabfdef)    : ', trim(fabfdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfABFCVs
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABF Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' #samp. in block averages (fblock_size)  : ', fblock_size
    write(PMF_OUT,130)  ' ABF force interpolation mode (fintrpl)  : ', fintrpl
    write(PMF_OUT,130)  ' ABF mask mode (fmask_mode)              : ', fmask_mode
    write(PMF_OUT,125)  ' ABF mask file (fabfmask)                : ', trim(fabfmask)
    write(PMF_OUT,125)  ' Apply ABF force (fapply_abf)            : ', prmfile_onoff(fapply_abf)

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
    if( use_key ) then
    write(PMF_OUT,125)  ' Server key file name (fserverkey)            : ', trim(fserverkey)
    else
    write(PMF_OUT,125)  ' Server URL (fserver)                         : ', trim(fserver)
    write(PMF_OUT,125)  ' Server password (fpassword)                  : ', trim(fpassword)
    end if
    else
    write(PMF_OUT,125)  ' Server URL (fserver)                         : ', 'none'
    write(PMF_OUT,125)  ' Server password (fpassword)                  : ', 'none'
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
    call abf2_cvs_cv_info(ABFCVList(i))
    write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

 return

120 format(A)
125 format(A,A)
130 format(A,I6)

140 format(' == Collective variable #',I4.4)

end subroutine abf2_init_print_header

!===============================================================================
! Subroutine:  abf2_init_arrays
!===============================================================================

subroutine abf2_init_arrays

    use pmf_utils
    use pmf_dat
    use abf2_dat
    use abf2_accumulator

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
            '[ABF2] Unable to allocate memory for arrays used in ABF calculation!')
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
                 '[ABF2] Unable to allocate memory for arrays used in Z matrix inversion!')
        end if
    end if

    ! init accumulator ------------------------------
    call abf2_accumulator_init

end subroutine abf2_init_arrays

!===============================================================================

end module abf2_init
