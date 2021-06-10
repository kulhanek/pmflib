!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz 
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtd_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_output_open
!===============================================================================

subroutine mtd_output_open

    use mtd_dat
    use pmf_dat
    use pmf_utils
    use pmf_constants

    implicit none
    integer        :: hills_starts,fstat
    logical        :: eread
    ! --------------------------------------------------------------------------

    ! open output file
    call pmf_utils_open(MTD_OUT,fmtdout,'R')

    ! method header --------------------------------
    write(MTD_OUT,10)
    select case(fmode)
        case(1,2,3)
            write(MTD_OUT,20)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmode!')
    end select
    write(MTD_OUT,30)

    ! extended output -------------------------------
    select case(fextout)
        case(0)
            ! do nothing
        case(1)
            ! open both files - as new one
            call pmf_utils_open(MTD_CVS,fmtdcvs,'R')
            call pmf_utils_open(MTD_HILLS,fmtdhills,'R')
            cvs_starts = 0
        case(2)
            ! open both files - and append new contents
            call pmf_utils_open(MTD_CVS,fmtdcvs,'U')
            ! rewind to the end and determine the number of hills
            eread       = .true.
            cvs_starts  = 0
            do while(eread)
                read(MTD_CVS,*,iostat = fstat) cvs_starts
                eread = fstat .eq. 0
            end do
            call pmf_utils_open(MTD_HILLS,fmtdhills,'U')
            eread        = .true.
            hills_starts = 0
            do while(eread)
                read(MTD_HILLS,*,iostat = fstat) hills_starts
                eread = fstat .eq. 0
            end do
            if( hills_starts .ne. cvs_starts ) then
                write(PMF_OUT,100) cvs_starts, hills_starts
                call pmf_utils_exit(PMF_OUT,1)
            end if
    end select

 return

 10 format('#===============================================================================')
 20 format('# Direct Metadynamics')
 30 format('#===============================================================================')

100 format('[MTD] CVS file contains ',I9,' records but HILLS file contains only ',I9,' records!')

end subroutine mtd_output_open

!===============================================================================
! Subroutine:  mtd_output_write_header
!===============================================================================

subroutine mtd_output_write_header

    use pmf_constants
    use pmf_dat
    use mtd_dat
    use pmf_cvs
    use pmf_utils

    implicit none
    integer        :: i,off
    ! -----------------------------------------------------------------------------

    write(MTD_OUT,*)

    ! data header - related to mtddynamics type
    select case(fmode)
        case(1,2,3)
            write(MTD_OUT,100,advance='NO') '#       1'
            off = 1
            do i=off+1,off+NumOfMTDCVs
                write(MTD_OUT,205,advance='NO') i
            end do
            write(MTD_OUT,*)

            write(MTD_OUT,100,advance='NO') '#   NSTEP'
            do i=1,NumOfMTDCVs
                write(MTD_OUT,200,advance='NO') trim(MTDCVList(i)%cv%name)

            end do
            write(MTD_OUT,*)

            write(MTD_OUT,100,advance='NO') '#        '
            do i=1,NumOfMTDCVs
                write(MTD_OUT,201,advance='NO') trim(MTDCVList(i)%cv%get_ulabel())

            end do
            write(MTD_OUT,*)

            write(MTD_OUT,100,advance='NO') '#--------'
            do i=1,NumOfMTDCVs
                write(MTD_OUT,300,advance='NO') '------------------'
            end do
            write(MTD_OUT,*)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmode!')
    end select

    return

100 format(A)
200 format(1X,A15)
201 format(1X,A15)
205 format(13X,I3)
300 format(1X,A15)

end subroutine mtd_output_write_header

!===============================================================================
! Subroutine:  mtd_output_write
!===============================================================================

subroutine mtd_output_write

 use pmf_constants
 use pmf_dat
 use mtd_dat
 use pmf_utils

 implicit none
 ! -----------------------------------------------------------------------------

 if( fpsample .le. 0 ) return ! output is written only if fpsample > 0
 if( mod(fstep,fpsample) .ne. 0 ) return

 select case(fmode)
    case(0)
        return
    case(1,2,3)
        call mtd_output_write_direct
    case default
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmode!')
 end select

 return

end subroutine mtd_output_write

!===============================================================================
! Subroutine:  mtd_output_write_direct
!===============================================================================

subroutine mtd_output_write_direct

 use pmf_constants
 use pmf_dat
 use mtd_dat
 use pmf_cvs

 implicit none
 integer        :: i
 ! -----------------------------------------------------------------------------

 write(MTD_OUT,10,advance='NO') fstep

 do i=1,NumOfMTDCVs
     write(MTD_OUT,20,advance='NO') &
           MTDCVList(i)%cv%get_rvalue(CVContext%CVsValues(MTDCVList(i)%cvindx))
 end do

 write(MTD_OUT,*)

 return

10 format(I9)
20 format(1X,F15.8)

end subroutine mtd_output_write_direct

!===============================================================================
! Subroutine:  mtd_output_close
!===============================================================================

subroutine mtd_output_close

 use mtd_dat
 use pmf_dat
 use pmf_utils
 use pmf_constants

 implicit none
 ! -----------------------------------------------------------------------------

 if( .not. mtd_enabled ) return ! method is not enabled

 select case(fmode)
    case(1,2,3)
        ! nothing to do
    case default
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported fmode!')
 end select


 select case(fextout)
    case(0)
        ! do nothing
    case(1,2)
        close(MTD_CVS)
        close(MTD_HILLS)
 end select

 close(MTD_OUT)

 return

end subroutine mtd_output_close

!===============================================================================

end module mtd_output
