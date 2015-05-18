!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pdrv_init

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  pdrv_init_method
!===============================================================================

subroutine pdrv_init_method

    use pdrv_output

    implicit none
    ! --------------------------------------------------------------------------

    call pdrv_init_print_header
    call pdrv_output_open
    call pdrv_output_write_header

end subroutine pdrv_init_method

!===============================================================================
! Subroutine:  pdrv_init_dat
!===============================================================================

subroutine pdrv_init_dat

    use pdrv_dat

    implicit none
    ! --------------------------------------------------------------------------

    fmode              = 0         ! 0 - disable PDRV, 1 - enabled PDRV
    fsample            = 500       ! output sample pariod in steps

    NumOfPDRVItems      = 0         ! number of monitored CVs

end subroutine pdrv_init_dat

!===============================================================================
! Subroutine:  pdrv_init_print_header
!===============================================================================

subroutine pdrv_init_print_header

    use pdrv_dat
    use pmf_dat
    use pmf_utils
    use pmf_paths
    use pdrv_paths

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    write(PMF_OUT,120)
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)  ' **************************** PATH DRIVING METHOD ***************************** '
    write(PMF_OUT,120)  '================================================================================'
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Path Driving Mode'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,130)  ' Path driving mode (fmode)               : ', fmode
    write(PMF_OUT,130)  ' Number of paths                         : ', NumOfPDRVItems
    write(PMF_OUT,125)  ' Path driving definition file (fpdrvdef) : ', trim(fpdrvdef)
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' Output options:'
    write(PMF_OUT,120)  ' ------------------------------------------------------'
    write(PMF_OUT,125)  ' Output file (fpdrvout)                  : ', trim(fpdrvout)
    write(PMF_OUT,130)  ' Output sampling (fsample)               : ', fsample
    write(PMF_OUT,120)
    write(PMF_OUT,120)  ' List of driven paths'
    write(PMF_OUT,120)  ' -------------------------------------------------------'
    write(PMF_OUT,120)

    do i=1,NumOfPDRVItems
        write(PMF_OUT,140) i
        call pdrv_paths_pdrv_info(PDRVCVList(i))
        write(PMF_OUT,120)
    enddo

    write(PMF_OUT,120)  '================================================================================'

    return

120 format(A)
125 format(A,A)
130 format(A,I6)

140 format(' == Driven path #',I2.2)

end subroutine pdrv_init_print_header

!===============================================================================

end module pdrv_init
