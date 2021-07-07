!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module usabf_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  usabf_output_open
!===============================================================================

subroutine usabf_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use usabf_dat

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(USABF_OUT,fusabfout,'R')

    write(USABF_OUT,10)
    write(USABF_OUT,20)
    write(USABF_OUT,30)

    return

10 format('#===============================================================================')
20 format('# Umbrella Sampling / Adaptive Biasing Force (US-ABF) Method')
30 format('#===============================================================================')

end subroutine usabf_output_open

!===============================================================================
! Subroutine:  usabf_output_write_header
!===============================================================================

subroutine usabf_output_write_header

    use pmf_constants
    use pmf_dat
    use usabf_dat
    use pmf_cvs

    implicit none
    integer        :: i, off
    ! --------------------------------------------------------------------------

    write(USABF_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfUSABFCVs
        write(USABF_OUT,15,advance='NO') i
    end do
    write(USABF_OUT,*)

    write(USABF_OUT,10,advance='NO') '#   NSTEP'
    do i=1,NumOfUSABFCVs
        write(USABF_OUT,20,advance='NO') trim(USABFCVList(i)%cv%name)
    end do
    write(USABF_OUT,*)

    write(USABF_OUT,10,advance='NO') '#        '
    do i=1,NumOfUSABFCVs
        write(USABF_OUT,25,advance='NO') '[' // trim(USABFCVList(i)%cv%get_ulabel()) // ']'
    end do
    write(USABF_OUT,*)

    write(USABF_OUT,10,advance='NO') '#--------'
    do i=1,NumOfUSABFCVs
        write(USABF_OUT,30,advance='NO') '---------------'
    end do
    write(USABF_OUT,*)

    return

10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
25 format(1X,A15)
30 format(1X,A15)

end subroutine usabf_output_write_header

!===============================================================================
! Subroutine:  usabf_output_write
!===============================================================================

subroutine usabf_output_write

    use pmf_constants
    use pmf_dat
    use usabf_dat
    use pmf_cvs

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(USABF_OUT,10,advance='NO') fstep

    do i=1,NumOfUSABFCVs
        write(USABF_OUT,20,advance='NO') &
                USABFCVList(i)%cv%get_rvalue(CVContext%CVsValues(USABFCVList(i)%cvindx))
    end do

    write(USABF_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)

end subroutine usabf_output_write

!===============================================================================
! Subroutine:  usabf_output_close
!===============================================================================

subroutine usabf_output_close

    use pmf_constants
    use pmf_dat
    use usabf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(USABF_OUT)

    return

end subroutine usabf_output_close

!===============================================================================

end module usabf_output
