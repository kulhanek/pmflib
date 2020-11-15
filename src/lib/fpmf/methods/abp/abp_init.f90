!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_init_method
!===============================================================================

subroutine abp_init_method

 use abp_output
 use abp_restart
 use abp_trajectory

 implicit none
 ! -----------------------------------------------------------------------------

 call abp_init_arrays
 call abp_init_print_header
 call abp_output_open
 call abp_restart_read
 call abp_trajectory_open
 call abp_output_write_header

end subroutine abp_init_method

!===============================================================================
! Subroutine:  abp_init_dat
!===============================================================================

subroutine abp_init_dat

    use abp_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode              = 0          ! 0 - disable ABP, 1 - enabled ABP
    fsample            = 5000       ! output sample period in steps
    frestart           = .false.    ! 1 - restart job with previous data, 0 - otherwise not
    frstupdate         = 5000       ! how often is restart file written
    feimode            = 1          ! 1 - linear ramp I
    ftrjsample         = 0          ! how often save accumulator to "accumulator evolution"
    fhbias             = 3.0d0      ! bias height

    fhramp             = 10         ! definition of linear ramp

    NumOfABPCVs        = 0

end subroutine abp_init_dat

!===============================================================================
! Subroutine:  abp_print_header
!===============================================================================

subroutine abp_init_print_header

    use prmfile
    use pmf_constants
    use pmf_dat
    use abp_dat
    use abp_cvs
    use pmf_utils
    use pmf_unit

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ********************* ADAPTIVE BIASING POTENTIAL METHOD ********************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABP Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' ABP mode (fmode)                        : ', fmode
    write(PMF_OUT,125)  ' Coordinate definition file (fabpdef)    : ', trim(fabpdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfABPCVs
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' ABP Control'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,136)  ' Bias height (fhbias)                    : ', &
                          pmf_unit_get_rvalue(EnergyUnit,fhbias),trim(pmf_unit_label(EnergyUnit))
    write(PMF_OUT,130)  ' Extrapolation/interpolation mode        : ', feimode
    select case(feimode)
        case(1)
    write(PMF_OUT,130)  ' Linear ramp size (fhramp)               : ', fhramp
        case default
    call pmf_utils_exit(PMF_OUT,1,'[ABP] Unknown extrapolation/interpolation mode!')
    end select

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Restart options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Restart file (fabprst)                  : ', trim(fabprst)
    write(PMF_OUT,125)  ' Restart enabled (frestart)              : ', prmfile_onoff(frestart)
    write(PMF_OUT,130)  ' Final restart update (frstupdate)       : ', frstupdate
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fabpout)                   : ', trim(fabpout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Trajectory output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Trajectory sampling (ftrjsample)        : ', ftrjsample
    write(PMF_OUT,125)  ' Trajectory file (fabptrj)               : ', trim(fabptrj)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of ABP coordinates'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfABPCVs
        write(PMF_OUT,140) i
        call abp_cvs_cv_info(ABPCVList(i))
        write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

    ! setup data
    kt   = ftemp * PMF_Rgas
    cfac = fhbias/kt

    return

    120 format(A)
    125 format(A,A)
    130 format(A,I6)
    136 format(A,E12.5,' [',A,']')

    140 format(' == Collective variable #',I4.4)

end subroutine abp_init_print_header

!===============================================================================
! Subroutine:  abp_init_arrays
!===============================================================================

subroutine abp_init_arrays

    use pmf_utils
    use pmf_dat
    use abp_dat
    use abp_accumulator

    implicit none
    integer        :: alloc_failed
    ! --------------------------------------------------------------------------

    ! general arrays --------------------------------
    allocate(la(NumOfABPCVs),               &
            cvvalues(NumOfABPCVs),          &
            cvindx(NumOfABPCVs),            &
            gridindx(NumOfABPCVs),          &
            gridvalues(NumOfABPCVs),        &
            diffvalues(NumOfABPCVs),        &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
                '[ABP] Unable to allocate memory for arrays used in ABP calculation!')
    end if

    la(:) = 0.0d0
    cvvalues(:) = 0.0d0

    ! init accumulator ------------------------------
    call abp_accumulator_init

end subroutine abp_init_arrays

!===============================================================================

end module abp_init
