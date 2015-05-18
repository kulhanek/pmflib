!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module mon_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mon_output_open
!===============================================================================

subroutine mon_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(MON_OUT,fmonout,'R')

    write(MON_OUT,10)
    write(MON_OUT,20)
    write(MON_OUT,30)

    return

10 format('#===============================================================================')
20 format('# Monitoring Method                                                             ')
30 format('#===============================================================================')

end subroutine mon_output_open

!===============================================================================
! Subroutine:  mon_output_write_header
!===============================================================================

subroutine mon_output_write_header

    use pmf_constants
    use pmf_dat
    use mon_dat
    use pmf_cvs

    implicit none
    integer        :: i,off
    ! --------------------------------------------------------------------------

    write(MON_OUT,10,advance='NO') '#  NSTEP '
    do i=1,NumOfMONItems
        write(MON_OUT,20,advance='NO') trim(MONCVList(i)%cv%name)
    end do
    write(MON_OUT,*)

    write(MON_OUT,10,advance='NO') '#        '
    do i=1,NumOfMONItems
        write(MON_OUT,30,advance='NO') '['//trim(MONCVList(i)%cv%get_ulabel())//']'
    end do
    write(MON_OUT,*)

    write(MON_OUT,10,advance='NO') '#--------'
    do i=1,NumOfMONItems
        write(MON_OUT,40,advance='NO') '---------------'
    end do
    write(MON_OUT,*)

    write(MON_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfMONItems
        write(MON_OUT,15,advance='NO') i
    end do
    write(MON_OUT,*)

    write(MON_OUT,10,advance='NO') '#--------'
    do i=1,NumOfMONItems
        write(MON_OUT,40,advance='NO') '---------------'
    end do
    write(MON_OUT,*)

    return

10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
30 format(1X,A15)
40 format(1X,A15)

end subroutine mon_output_write_header

!===============================================================================
! Subroutine:  mon_output_write_output
!===============================================================================

subroutine mon_output_write_output

    use pmf_constants
    use pmf_dat
    use mon_dat
    use pmf_cvs

    implicit none
    integer        :: i
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(MON_OUT,10,advance='NO') fstep

    do i=1,NumOfMONItems
     write(MON_OUT,20,advance='NO') &
            MONCVList(i)%cv%get_rvalue(CVContext%CVsValues(MONCVList(i)%cvindx))
    end do
    write(MON_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)

end subroutine mon_output_write_output

!===============================================================================
! Subroutine:  mon_output_close
!===============================================================================

subroutine mon_output_close

    use pmf_constants
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(MON_OUT)

    return

end subroutine mon_output_close

!===============================================================================

end module mon_output
