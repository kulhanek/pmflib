!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_output_open
!===============================================================================

subroutine abp_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(ABP_OUT,fabpout,'R')

    write(ABP_OUT,10)
    write(ABP_OUT,20)
    write(ABP_OUT,30)
    write(ABP_OUT,*)

    return

10 format('#===============================================================================')
20 format('# Adaptive Biasing Potential Method')
30 format('#===============================================================================')

end subroutine abp_output_open

!===============================================================================
! Subroutine:  abp_output_write_header
!===============================================================================

subroutine abp_output_write_header

    use pmf_constants
    use pmf_dat
    use abp_dat
    use pmf_cvs

    implicit none
    integer        :: i, off, ci
    ! --------------------------------------------------------------------------

    write(ABP_OUT,*)

    write(ABP_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfABPCVs
        write(ABP_OUT,15,advance='NO') i
    end do
    write(ABP_OUT,*)

    write(ABP_OUT,10,advance='NO') '#  NSTEP '
    do i=1,NumOfABPCVs
        ci = ABPCVList(i)%cvindx
        write(ABP_OUT,20,advance='NO') CVList(ci)%cv%name
    end do
    write(ABP_OUT,*)

    write(ABP_OUT,10,advance='NO') '#        '
    do i=1,NumOfABPCVs
        ci = ABPCVList(i)%cvindx
        write(ABP_OUT,25,advance='NO') trim(CVList(ci)%cv%get_ulabel())
    end do
    write(ABP_OUT,*)

    write(ABP_OUT,10,advance='NO') '#--------'
    do i=1,NumOfABPCVs
        write(ABP_OUT,30,advance='NO') '---------------'
    end do
    write(ABP_OUT,*)

    return

10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
25 format(1X,'[',A13,']')
30 format(1X,A15)

end subroutine abp_output_write_header

!===============================================================================
! Subroutine:  abp_output_write
!===============================================================================

subroutine abp_output_write

    use pmf_constants
    use pmf_dat
    use abp_dat
    use pmf_cvs

    implicit none
    integer        :: i,ci
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0

    if( mod(fstep,fsample) .ne. 0 ) return

    write(ABP_OUT,10,advance='NO') fstep

    do i=1,NumOfABPCVs
        ci = ABPCVList(i)%cvindx
        write(ABP_OUT,20,advance='NO') ABPCVList(i)%cv%get_rvalue(CVContext%CVsValues(ci))
    end do

    write(ABP_OUT,*)

    return

    10 format(I9)
    20 format(1X,F15.8)

end subroutine abp_output_write

!===============================================================================
! Subroutine:  abp_output_close
!===============================================================================

subroutine abp_output_close

    use pmf_constants
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(ABP_OUT)

    return

end subroutine abp_output_close

!===============================================================================

end module abp_output
