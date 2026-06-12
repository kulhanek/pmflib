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

module mtc_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtc_output_open
!===============================================================================

subroutine mtc_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(MTC_OUT,fmtcout,'R')

    write(MTC_OUT,10)
    write(MTC_OUT,20)
    write(MTC_OUT,30)

    return

10 format('#===============================================================================')
20 format('# Metric Tensor Correction                                                      ')
30 format('#===============================================================================')

end subroutine mtc_output_open

!===============================================================================
! Subroutine:  mtc_output_write_header
!===============================================================================

subroutine mtc_output_write_header

    use pmf_constants
    use pmf_dat
    use mtc_dat
    use pmf_cvs

    implicit none
    integer        :: i,off
    ! --------------------------------------------------------------------------

    write(MTC_OUT,1) '#'
    write(MTC_OUT,10,advance='NO') '#  NSTEP '
    do i=1,NumOfMTCCVs
        write(MTC_OUT,20,advance='NO') trim(MTCCVList(i)%cv%name)
    end do
    write(MTC_OUT,20,advance='NO') 'MTC'
    write(MTC_OUT,*)

    write(MTC_OUT,10,advance='NO') '#        '
    do i=1,NumOfMTCCVs
        write(MTC_OUT,30,advance='NO') '['//trim(MTCCVList(i)%cv%get_ulabel())//']'
    end do
    write(MTC_OUT,20,advance='NO') '[i.u.]'
    write(MTC_OUT,*)

    write(MTC_OUT,10,advance='NO') '#--------'
    do i=1,NumOfMTCCVs
        write(MTC_OUT,40,advance='NO') '---------------'
    end do
    write(MTC_OUT,40,advance='NO') '---------------'
    write(MTC_OUT,*)

    write(MTC_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfMTCCVs
        write(MTC_OUT,15,advance='NO') i
    end do
    write(MTC_OUT,*)

    write(MTC_OUT,10,advance='NO') '#--------'
    do i=1,NumOfMTCCVs+1
        write(MTC_OUT,40,advance='NO') '---------------'
    end do
    write(MTC_OUT,*)

    flush(MTC_OUT)

    return

 1 format(A)
10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
30 format(1X,A15)
40 format(1X,A15)

end subroutine mtc_output_write_header

!===============================================================================
! Subroutine:  mtc_output_write_output
!===============================================================================

subroutine mtc_output_write_output

    use pmf_constants
    use pmf_dat
    use mtc_dat
    use pmf_cvs

    implicit none
    integer         :: i
    real(PMFDP)     :: mtc
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(MTC_OUT,10,advance='NO') fstep

    mtc = sqrt(fzdet)

    do i=1,NumOfMTCCVs
         write(MTC_OUT,20,advance='NO') &
            MTCCVList(i)%cv%get_rvalue(CVContext%CVsValues(MTCCVList(i)%cvindx))
    end do
    write(MTC_OUT,20,advance='NO') mtc
    write(MTC_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)

end subroutine mtc_output_write_output

!===============================================================================
! Subroutine:  mtc_output_close
!===============================================================================

subroutine mtc_output_close

    use pmf_constants
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(MTC_OUT)

    return

end subroutine mtc_output_close

!===============================================================================

end module mtc_output
