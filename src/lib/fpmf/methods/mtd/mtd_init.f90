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

module mtd_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_init_method
!===============================================================================

subroutine mtd_init_method

    use mtd_output
    use mtd_restart
    use mtd_client
    use mtd_accu
    use mtd_trajectory

    implicit none
    ! -----------------------------------------------------------------------------

    call mtd_init_core
    call mtd_init_print_header
    call mtd_accu_init
    call mtd_output_open
    call mtd_trajectory_open
    call mtd_restart_read
    call mtd_client_register
    call mtd_client_get_initial_data
    call mtd_output_write_header

end subroutine mtd_init_method

!===============================================================================
! Subroutine:  mtd_init_dat
!===============================================================================

subroutine mtd_init_dat

    use mtd_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! control section ----------------------------------------------------------
    fmode           = 0
    fmetastep       = 1000
    fheight         = 0.10d0
    fmetatemp       = 0.0d0
    fsample         = 5000
    fwritehills     = .false.

    frestart        = .false.
    frstupdate      = 5000
    ftrjsample      = 0

    ! server part --------------------------------------------------------------
    fserver_enabled = .false.
    fserverkey      = ''
    fserver         = ''
    fserverupdate   = 500
    fconrepeats     = 0
    fabortonmwaerr  = .true.

    insidesamples   = 0
    outsidesamples  = 0

    numofhills      = 0

    NumOfMTDCVs     = 0

end subroutine mtd_init_dat

!===============================================================================
! Subroutine:  mtd_init_print_header
!===============================================================================

subroutine mtd_init_print_header

    use mtd_dat
    use pmf_dat
    use pmf_constants
    use mtd_cvs_mod
    use prmfile

    implicit none
    integer        :: i
    ! -----------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ******************************* METADYNAMICS ********************************* '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Metadynamics Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Metadynamics mode (fmode)               : ', fmode
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfMTDCVs
    write(PMF_OUT,125)  ' Coordinate definition file (fmtddef)    : ', trim(fmtddef)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Metadynamics Controls:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Metadynamics step size (fmetastep)      : ', fmetastep
    write(PMF_OUT,135)  ' Height (fheight)                        : ', pmf_unit_get_rvalue(EnergyUnit,fheight), &
                                                                       '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(PMF_OUT,135)  ' Meta temperature (fmetatemp)            : ', pmf_unit_get_rvalue(TemperatureUnit,fmetatemp), &
                                                                       '['//trim(pmf_unit_label(TemperatureUnit))//']'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fmtdrst)                  : ', trim(fmtdrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Restart file update (frstupdate)        : ', frstupdate
    write(PMF_OUT,120)

    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fmtdout)                   : ', trim(fmtdout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,125)  ' Output file (fmtdhills)                 : ', trim(fmtdhills)
    write(PMF_OUT,125)  ' Output MTD hills (fwritehills)          : ', prmfile_onoff(fwritehills)
    write(PMF_OUT,120)

    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Trajectory file (fmtdtrj)               : ', trim(fmtdtrj)
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample

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

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of MTD collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfMTDCVs
    write(PMF_OUT,140) i
    call mtd_cvs_cv_info(MTDCVList(i))
    write(PMF_OUT,120)
    enddo

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)
135 format(A,E12.5,1X,A)

140 format(' == Collective variable #',I4.4)

end subroutine mtd_init_print_header

!===============================================================================
! Subroutine:  mtd_init_core
!===============================================================================

subroutine mtd_init_core

    use pmf_dat
    use pmf_utils
    use mtd_dat
    use mtd_accu

    implicit none
    integer     :: alloc_failed
    ! -----------------------------------------------------------------------------

    meta_next_fstep = fmetastep ! initial value

! general arrays --------------------------------
    allocate(                               &
            CVValues(NumOfMTDCVs),          &
            MTDForce(NumOfMTDCVs),          &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[MTD] Unable to allocate memory for arrays used in MTD calculation!')
    end if

    CVValues(:) = 0.0d0
    MTDForce(:) = 0.0d0

    ! init accumulator ------------------------------
    call mtd_accu_init

    return

end subroutine mtd_init_core

!===============================================================================

end module mtd_init
