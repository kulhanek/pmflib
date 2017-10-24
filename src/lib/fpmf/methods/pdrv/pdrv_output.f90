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

module pdrv_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  pdrv_output_open
!===============================================================================

subroutine pdrv_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(PDRV_OUT,fpdrvout,'R')

    write(PDRV_OUT,10)
    write(PDRV_OUT,20)
    write(PDRV_OUT,30)
    write(PDRV_OUT,40)

    return

10 format('#===============================================================================')
20 format('# Path Driving Method                                                           ')
30 format('#===============================================================================')
40 format('# NOTE: The requested values are reported from the previous time step!')

end subroutine pdrv_output_open

!===============================================================================
! Subroutine:  pdrv_output_write_header
!===============================================================================

subroutine pdrv_output_write_header

    use pmf_constants
    use pmf_dat
    use pdrv_dat
    use pmf_cvs

    implicit none
    integer        :: i,j,off
    ! --------------------------------------------------------------------------

    write(PDRV_OUT,15,advance='NO')
    do i=1,NumOfPDRVItems
        write(PDRV_OUT,55,advance='NO') trim(PDRVCVList(i)%path%name)
        do j=1,PDRVCVList(i)%path%ncvs
            write(PDRV_OUT,50,advance='NO') trim(PDRVCVList(i)%path%name) // ' / ' &
                                            // trim(PDRVCVList(i)%path%cvs(j)%cv%name) &
                                            // ' ['//trim(PDRVCVList(i)%path%cvs(j)%cv%get_ulabel())//']'
        end do
    end do
    write(PDRV_OUT,*)

    write(PDRV_OUT,15,advance='NO')
    do i=1,NumOfPDRVItems
        write(PDRV_OUT,11,advance='NO')
        do j=1,PDRVCVList(i)%path%ncvs
            write(PDRV_OUT,20,advance='NO')
        end do
    end do
    write(PDRV_OUT,*)

    write(PDRV_OUT,5,advance='NO')
    do i=1,NumOfPDRVItems
        write(PDRV_OUT,6,advance='NO')
        do j=1,PDRVCVList(i)%path%ncvs
            write(PDRV_OUT,30,advance='NO')
        end do
    end do
    write(PDRV_OUT,*)

    write(PDRV_OUT,10,advance='NO')
    do i=1,NumOfPDRVItems
        write(PDRV_OUT,12,advance='NO')
        do j=1,PDRVCVList(i)%path%ncvs
            write(PDRV_OUT,35,advance='NO')
        end do
    end do
    write(PDRV_OUT,*)

    write(PDRV_OUT,2,advance='NO')
    off = 2
    do i=1,NumOfPDRVItems
        write(PDRV_OUT,45,advance='NO') off
        off = off + 1
        write(PDRV_OUT,45,advance='NO') off
        off = off + 1
        write(PDRV_OUT,46,advance='NO') off
        off = off + 1
        do j=1,PDRVCVList(i)%path%ncvs*3
            write(PDRV_OUT,40,advance='NO') off
            off = off + 1
        end do
    end do
    write(PDRV_OUT,*)

    write(PDRV_OUT,10,advance='NO')
    do i=1,NumOfPDRVItems
        write(PDRV_OUT,12,advance='NO')
        do j=1,PDRVCVList(i)%path%ncvs
            write(PDRV_OUT,35,advance='NO')
        end do
    end do
    write(PDRV_OUT,*)

    return

 2 format('#       1')
 5 format('#  NSTEP ')
10 format('# -------')
11 format(' ------------------------------')
 6 format(' ReqAlpha CurrAlph  DevToPath  ')
12 format(' -------- -------- ------------')
15 format('#        ')
20 format(' -----------------------------------------------')
30 format(' requested value  current value    difference   ')
35 format(' --------------- --------------- ---------------')
40 format(1X,I15)
45 format(1X,I8)
46 format(1X,I12)
50 format(1X,A47)
55 format(1X,A17)

end subroutine pdrv_output_write_header

!===============================================================================
! Subroutine:  pdrv_output_write_output
!===============================================================================

subroutine pdrv_output_write_output

    use pmf_constants
    use pmf_dat
    use pdrv_dat
    use pmf_cvs
    use pmf_paths
    use pmf_timers

    implicit none
    integer     :: i, j
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(PDRV_OUT,10,advance='NO') fstep

    if( .not. fmonitor_paths ) then
        ! if the path monitoring is not enabled calculate the aplha only for driven paths
        call pmf_timers_start_timer(PMFLIB_PATH_TIMER)
            ! update PATHs
            do i=1,NumOfPDRVItems
                call pmf_paths_get_path_current_alpha(PDRVCVList(i)%path,CVContext)
            end do
        call pmf_timers_stop_timer(PMFLIB_PATH_TIMER)
    end if

    do i=1,NumOfPDRVItems
        write(PDRV_OUT,20,advance='NO') PDRVCVList(i)%req_alpha
        write(PDRV_OUT,20,advance='NO') PDRVCVList(i)%path%current_alpha
        write(PDRV_OUT,25,advance='NO') pmf_paths_get_dist(PDRVCVList(i)%path,CVContext,PDRVCVList(i)%path%current_alpha)
        do j=1,PDRVCVList(i)%path%ncvs
            write(PDRV_OUT,30,advance='NO') &
                PDRVCVList(i)%path%cvs(j)%cv%get_rvalue(PDRVCVList(i)%path%rpos(j))
            write(PDRV_OUT,30,advance='NO') &
                PDRVCVList(i)%path%cvs(j)%cv%get_rvalue(PDRVCVList(i)%path%cpos(j))
            write(PDRV_OUT,30,advance='NO') &
                PDRVCVList(i)%path%cvs(j)%cv%get_rvalue(PDRVCVList(i)%path%rpos(j)-PDRVCVList(i)%path%cpos(j))
        end do
    end do
    write(PDRV_OUT,*)

    return

10 format(I9)
20 format(1X,F8.4)
25 format(1X,E12.5)
30 format(1X,F15.8)


end subroutine pdrv_output_write_output

!===============================================================================
! Subroutine:  pdrv_output_close
!===============================================================================

subroutine pdrv_output_close

    use pmf_constants
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(PDRV_OUT)

    return

end subroutine pdrv_output_close

!===============================================================================

end module pdrv_output
