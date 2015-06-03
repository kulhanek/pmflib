!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module gap_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  gap_output_open
!===============================================================================

subroutine gap_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_utils_open(GAP_OUT,fgapout,'R')

    write(GAP_OUT,10)
    write(GAP_OUT,20)
    write(GAP_OUT,30)

    return

10 format('#===============================================================================')
20 format('# GAP Method                                                                    ')
30 format('#===============================================================================')

end subroutine gap_output_open

!===============================================================================
! Subroutine:  gap_output_write_header
!===============================================================================

subroutine gap_output_write_header

    use pmf_constants
    use pmf_dat
    use pmf_cvs
    use gap_dat

    implicit none
    integer        :: off, k, i
    ! --------------------------------------------------------------------------

    write(GAP_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfGAPCVs+NumOfGAPGroups
        write(GAP_OUT,15,advance='NO') i
    end do
    write(GAP_OUT,*)

    write(GAP_OUT,10,advance='NO') '#  NSTEP '
    do k=1,NumOfGAPGroups
       do i=1,GAPGroupList(k)%nindexes
          write(GAP_OUT,20,advance='NO') trim(GAPCVList(GAPGroupList(k)%gapcvindx(i))%cv%name)
       end do
       write(GAP_OUT,20,advance='NO') 'E('//trim(GAPGroupList(k)%groupname)//')'
    end do
    write(GAP_OUT,*)

    write(GAP_OUT,10,advance='NO') '#        '
    do k=1,NumOfGAPGroups
       do i=1,GAPGroupList(k)%nindexes
          write(GAP_OUT,30,advance='NO') '['//trim(GAPCVList(GAPGroupList(k)%gapcvindx(i))%cv%get_ulabel())//']'
       end do
    end do
    write(GAP_OUT,30,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
    write(GAP_OUT,*)

    write(GAP_OUT,10,advance='NO') '#--------'

    do k=1,NumOfGAPGroups
       do i=1,GAPGroupList(k)%nindexes
          write(GAP_OUT,40,advance='NO') '---------------'
       end do
    end do

    write(GAP_OUT,40,advance='NO') '---------------'

    write(GAP_OUT,*)

    return

10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
30 format(1X,A15)
40 format(1X,A15)

end subroutine gap_output_write_header

!===============================================================================
! Subroutine:  gap_output_write
!===============================================================================

subroutine gap_output_write

    use pmf_constants
    use pmf_dat
    use gap_dat
    use pmf_cvs

    implicit none
    integer        :: k, i, ci
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(GAP_OUT,10,advance='NO') fstep

    do k=1,NumOfGAPGroups
       do i=1,GAPGroupList(k)%nindexes
          ci = GAPGroupList(k)%gapcvindx(i)
          write(GAP_OUT,20,advance='NO') GAPCVList(ci)%cv%get_rvalue(CVContext%CVsValues(GAPCVList(ci)%cvindx))
       end do
       write(GAP_OUT,30,advance='NO') GAPGroupList(k)%energy
    end do

    write(GAP_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)
30 format(1X,F15.8)

end subroutine gap_output_write

!===============================================================================
! Subroutine:  gap_output_close
!===============================================================================

subroutine gap_output_close

    use pmf_constants
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(GAP_OUT)

    return

end subroutine gap_output_close

!===============================================================================

end module gap_output
