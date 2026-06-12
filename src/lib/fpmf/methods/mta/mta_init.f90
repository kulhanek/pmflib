!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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
module mtc_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtc_init_method
!===============================================================================

subroutine mtc_init_method

    use mtc_output
    use mtc_restart

    implicit none
    ! --------------------------------------------------------------------------

    call mtc_init_print_header
    call mtc_init_arrays
    call mtc_output_open
    call mtc_restart_read
    call mtc_output_write_header

end subroutine mtc_init_method

!===============================================================================
! Subroutine:  mtc_init_dat
!===============================================================================

subroutine mtc_init_dat

    use mtc_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode               = 0         ! 0 - disable MTC, 1 - enabled MTC
    fsample             = 500       ! output sample pariod in steps
    frestart            = .false.
    frstupdate          = 5000

    NumOfMTCCVs         = 0         ! number of CVs

end subroutine mtc_init_dat

!===============================================================================
! Subroutine:  mtc_init_print_header
!===============================================================================

subroutine mtc_init_print_header

    use mtc_dat
    use pmf_dat
    use pmf_utils
    use pmf_cvs
    use prmfile
    use mtc_cvs

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ************************** METRIC TENSOR CORRECTION ************************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' MTC Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' MTC mode (fmode)                        : ', fmode
    write(PMF_OUT,130)  ' Number of collective variables          : ', NumOfMTCCVs
    write(PMF_OUT,125)  ' CV definition file (fmtcdef)            : ', trim(fmtcdef)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fmtcout)                   : ', trim(fmtcout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fmtcrst)                  : ', trim(fmtcrst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Restart file update (frstupdate)        : ', frstupdate
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of MTC collective variables'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfMTCCVs
        write(PMF_OUT,140) i
        call mtc_cvs_cv_info(MTCCVList(i))
        write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)
140 format(' == Collective variable #',I2.2)

end subroutine mtc_init_print_header

!===============================================================================
! Subroutine:  mtc_init_arrays
!===============================================================================

subroutine mtc_init_arrays

    use pmf_utils
    use pmf_dat
    use mtc_dat
    use mtc_accu

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

! init accumulator
    call mtc_accu_init

! general arrays --------------------------------
    allocate(                                   &
            fz(NumOfMTCCVs,NumOfMTCCVs),        &
            fzinv(NumOfMTCCVs,NumOfMTCCVs),     &
            indx(NumOfMTCCVs),                  &
            vv(NumOfMTCCVs),                    &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[MTC] Unable to allocate memory for arrays used in MTC calculation!')
    end if

    fz(:,:)     = 0.0d0
    fzinv(:,:)  = 0.0d0
    indx(:)     = 1.0d0
    vv(:)       = 1.0d0

end subroutine mtc_init_arrays

!===============================================================================

end module mtc_init
