!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module stm_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  stm_init_method
!===============================================================================

subroutine stm_init_method

    use stm_output
    use stm_client

    implicit none
    ! --------------------------------------------------------------------------

    call stm_init_arrays
    call stm_init_print_header
    call stm_output_open
    call stm_client_register
    !call stm_client_get_initial_data
    call stm_output_write_header

end subroutine stm_init_method

!===============================================================================
! Subroutine:  stm_init_dat
!===============================================================================

subroutine stm_init_dat

    use stm_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode           = 0         ! 0 - disable STM, 1 - enabled STM
    fsample         = 500       ! output sample pariod in steps
    fbeadidfile     = 'beadid'  ! bead id file definition
    ftensor         = 2         ! 0 - unity, 1 - normal, 2 - massweighted

    NumOfSTMCVs     = 0

    fserverkey      = ''
    fserver         = ''
    fpassword       = ''

    use_key         = .false.
    client_id       = -1

    stmmode         = BMO_UNKNOWN   ! current STM mode, 0 - undefined (U)
    bead_id         = -1            ! bead ID
    stmsteps        = 0             ! number of steps for current mode
    curstep         = 0             ! current step

end subroutine stm_init_dat

!===============================================================================
! Subroutine:  stm_print_header
!===============================================================================

subroutine stm_init_print_header

    use prmfile
    use pmf_constants
    use pmf_dat
    use stm_dat
    use stm_cvs

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' ******************************* String Method ******************************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' STM Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' STM mode (fmode)                        : ', fmode
    write(PMF_OUT,125)  ' Bead ID definition file (fbeadidfile)   : ', trim(fbeadidfile)
    write(PMF_OUT,125)  ' Coordinate definition file (fstmdef)    : ', trim(fstmdef)
    write(PMF_OUT,130)  ' Number of coordinates                   : ', NumOfSTMCVs
    write(PMF_OUT,130)  ' Bead ID                                 : ', bead_id
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fstmout)                   : ', trim(fstmout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' STM server options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    if( use_key ) then
    write(PMF_OUT,125)  ' Server key file name (fserverkey)            : ', trim(fserverkey)
    else
    write(PMF_OUT,125)  ' Server URL (fserver)                         : ', trim(fserver)
    write(PMF_OUT,125)  ' Server password (fpassword)                  : ', trim(fpassword)
    end if

    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of STM coordinates'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfSTMCVs
    write(PMF_OUT,140) i
    call stm_cvs_cv_info(STMCVList(i))
    write(PMF_OUT,120)
    end do

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)

140 format(' == Collective variable #',I2.2)

end subroutine stm_init_print_header

!===============================================================================
! Subroutine:  stm_init_arrays
!===============================================================================

subroutine stm_init_arrays

    use pmf_utils
    use pmf_dat
    use stm_dat

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

    ! general arrays --------------------------------
    allocate(beadpos(NumOfSTMCVs),          &
          PMF(NumOfSTMCVs),                 &
          MTZ(NumOfSTMCVs,NumOfSTMCVs),     &
          stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
             '[STM] Unable to allocate memory for arrays!')
    end if

    beadpos(:) = 0.0d0
    PMF(:) = 0.0d0
    MTZ(:,:) = 0.0d0

    return

end subroutine stm_init_arrays

!===============================================================================

end module stm_init
