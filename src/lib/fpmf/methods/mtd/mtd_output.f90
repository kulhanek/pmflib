!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    ! --------------------------------------------------------------------------

    ! open output file
    call pmf_utils_open(MTD_OUT,fmtdout,'R')

    ! method header --------------------------------
    write(MTD_OUT,10)
    write(MTD_OUT,20)
    write(MTD_OUT,30)

 return

 10 format('#===============================================================================')
 20 format('# Direct Metadynamics')
 30 format('#===============================================================================')

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


    write(MTD_OUT,100,advance='NO') '#   NSTEP'
    write(MTD_OUT,300,advance='NO') 'HEIGHT'
    do i=1,NumOfMTDCVs
        write(MTD_OUT,300,advance='NO') 'x('//trim(MTDCVList(i)%cv%name)//')'
        write(MTD_OUT,300,advance='NO') 'w('//trim(MTDCVList(i)%cv%name)//')'
    end do
    write(MTD_OUT,*)

    write(MTD_OUT,100,advance='NO') '#--------'
    write(MTD_OUT,300,advance='NO') '---------------'
    do i=1,NumOfMTDCVs
        write(MTD_OUT,300,advance='NO') '---------------'
        write(MTD_OUT,300,advance='NO') '---------------'
    end do
    write(MTD_OUT,*)

    write(MTD_OUT,100,advance='NO') '#        '
    write(MTD_OUT,300,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
    do i=1,NumOfMTDCVs
        write(MTD_OUT,305,advance='NO') '['//trim(MTDCVList(i)%cv%get_ulabel())//']'
    end do
    write(MTD_OUT,*)

    write(MTD_OUT,100,advance='NO') '#--------'
    write(MTD_OUT,300,advance='NO') '---------------'
    do i=1,NumOfMTDCVs
        write(MTD_OUT,100,advance='NO') ' -------------------------------'
    end do
    write(MTD_OUT,*)


    write(MTD_OUT,100,advance='NO') '#       1'
    write(MTD_OUT,205,advance='NO') 2
    off = 3
    do i=off+1,off+NumOfMTDCVs
        write(MTD_OUT,205,advance='NO') off
        off = off + 1
        write(MTD_OUT,205,advance='NO') off
        off = off + 1
    end do
    write(MTD_OUT,*)

    write(MTD_OUT,100,advance='NO') '#--------'
    write(MTD_OUT,300,advance='NO') '---------------'
    do i=1,NumOfMTDCVs
        write(MTD_OUT,300,advance='NO') '---------------'
        write(MTD_OUT,300,advance='NO') '---------------'
    end do
    write(MTD_OUT,*)

    return

100 format(A)
205 format(13X,I3)
300 format(1X,A15)
305 format(1X,A31)

end subroutine mtd_output_write_header

!===============================================================================
! Subroutine:  mtd_output_write
!===============================================================================

subroutine mtd_output_write(cvs,height,widths)

    use pmf_constants
    use pmf_dat
    use mtd_dat
    use pmf_utils

    implicit none
    real(PMFDP)     :: cvs(:)
    real(PMFDP)     :: height
    real(PMFDP)     :: widths(:)
    ! --------------------------------------------
    integer         :: i
    ! -----------------------------------------------------------------------------

    write(MTD_OUT,10,advance='NO') fstep
    write(MTD_OUT,15,advance='NO') height

    do i=1,NumOfMTDCVs
        write(MTD_OUT,20,advance='NO') &
            MTDCVList(i)%cv%get_rvalue(cvs(MTDCVList(i)%cvindx)), &
            MTDCVList(i)%cv%get_rvalue(widths(MTDCVList(i)%cvindx))
    end do

    write(MTD_OUT,*)

    return

10 format(I9)
15 format(1X,F15.8)
20 format(1X,F15.8,1X,F15.8)

end subroutine mtd_output_write

!===============================================================================
! Subroutine:  mtd_output_close
!===============================================================================

subroutine mtd_output_close

    use mtd_dat
    use pmf_dat
    use pmf_constants

    implicit none
    ! -----------------------------------------------------------------------------

    if( .not. mtd_enabled ) return ! method is not enabled
    close(MTD_OUT)

    return

end subroutine mtd_output_close

!===============================================================================

end module mtd_output
