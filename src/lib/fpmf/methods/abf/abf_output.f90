!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_output_open
!===============================================================================

subroutine abf_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use abf_dat

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(ABF_OUT,fabfout,'R')

    write(ABF_OUT,10)
    write(ABF_OUT,20)
    write(ABF_OUT,30)

    return

10 format('#===============================================================================')
20 format('# Adaptive Biasing Force Method')
30 format('#===============================================================================')

end subroutine abf_output_open

!===============================================================================
! Subroutine:  abf_output_write_header
!===============================================================================

subroutine abf_output_write_header

    use pmf_constants
    use pmf_dat
    use abf_dat
    use pmf_cvs

    implicit none
    integer        :: i, off
    ! --------------------------------------------------------------------------

    write(ABF_OUT,1) '#'
    write(ABF_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfABFCVs
        write(ABF_OUT,15,advance='NO') i
    end do
    write(ABF_OUT,*)

    write(ABF_OUT,10,advance='NO') '#   NSTEP'
    do i=1,NumOfABFCVs
        write(ABF_OUT,20,advance='NO') trim(ABFCVList(i)%cv%name)
    end do
    write(ABF_OUT,*)

    write(ABF_OUT,10,advance='NO') '#        '
    do i=1,NumOfABFCVs
        write(ABF_OUT,25,advance='NO') '[' // trim(ABFCVList(i)%cv%get_ulabel()) // ']'
    end do
    write(ABF_OUT,*)

    write(ABF_OUT,10,advance='NO') '#--------'
    do i=1,NumOfABFCVs
        write(ABF_OUT,30,advance='NO') '---------------'
    end do
    write(ABF_OUT,*)

    return

 1 format(A)
10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
25 format(1X,A15)
30 format(1X,A15)

end subroutine abf_output_write_header

!===============================================================================
! Subroutine:  abf_output_write
!===============================================================================

subroutine abf_output_write

    use pmf_constants
    use pmf_dat
    use abf_dat
    use pmf_cvs

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(ABF_OUT,10,advance='NO') fstep

    do i=1,NumOfABFCVs
        write(ABF_OUT,20,advance='NO') &
                ABFCVList(i)%cv%get_rvalue(CVContext%CVsValues(ABFCVList(i)%cvindx))
    end do

    write(ABF_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)

end subroutine abf_output_write

!===============================================================================
! Subroutine:  abf_output_close
!===============================================================================

subroutine abf_output_close

    use pmf_constants
    use pmf_dat
    use abf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(ABF_OUT)
    return

end subroutine abf_output_close

!===============================================================================

end module abf_output
