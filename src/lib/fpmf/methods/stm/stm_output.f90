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

module stm_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  stm_output_open
!===============================================================================

subroutine stm_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(STM_OUT,fstmout,'R')

    write(STM_OUT,10)
    write(STM_OUT,20)
    write(STM_OUT,30)

    return

10 format('#===============================================================================')
20 format('# String Method (STM)')
30 format('#===============================================================================')

end subroutine stm_output_open

!===============================================================================
! Subroutine:  stm_output_write_header
!===============================================================================

subroutine stm_output_write_header

    use pmf_constants
    use pmf_dat
    use stm_dat
    use pmf_cvs

    implicit none
    integer        :: i, off
    ! --------------------------------------------------------------------------

    write(STM_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfSTMCVs
        write(STM_OUT,15,advance='NO') i
    end do
    write(STM_OUT,*)

    write(STM_OUT,10,advance='NO') '#   NSTEP'
    do i=1,NumOfSTMCVs
        write(STM_OUT,20,advance='NO') trim(STMCVList(i)%cv%name)
    end do
    write(STM_OUT,*)

    write(STM_OUT,10,advance='NO') '#        '
    do i=1,NumOfSTMCVs
        write(STM_OUT,20,advance='NO') trim('['//trim(STMCVList(i)%cv%get_ulabel())//']')
    end do
    write(STM_OUT,*)

    write(STM_OUT,10,advance='NO') '#--------'
    do i=1,NumOfSTMCVs
        write(STM_OUT,30,advance='NO') '---------------'
    end do
    write(STM_OUT,*)

    return

10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
30 format(1X,A15)

end subroutine stm_output_write_header

!===============================================================================
! Subroutine:  stm_output_write
!===============================================================================

subroutine stm_output_write

    use pmf_constants
    use pmf_dat
    use stm_dat
    use pmf_cvs

    implicit none
    integer        :: i,ci
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0

    if( mod(fstep,fsample) .ne. 0 ) return

    write(STM_OUT,10,advance='NO') fstep

    do i=1,NumOfSTMCVs
        ci = STMCVList(i)%cvindx
        write(STM_OUT,20,advance='NO') STMCVList(i)%cv%get_rvalue(CVContext%CVsValues(STMCVList(i)%cvindx))
    end do

    write(STM_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)

end subroutine stm_output_write

!===============================================================================
! Subroutine:  stm_output_close
!===============================================================================

subroutine stm_output_close

    use pmf_constants
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(STM_OUT)

    return

end subroutine stm_output_close

!===============================================================================

end module stm_output
