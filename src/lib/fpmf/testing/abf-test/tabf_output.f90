!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module tabf_output

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  tabf_output_open
!===============================================================================

subroutine tabf_output_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use tabf_dat

    implicit none
    integer     :: alloc_failed
    ! --------------------------------------------------------------------------

    call pmf_utils_open(TABF_OUT,ftabfout,'R')

    write(TABF_OUT,10)
    write(TABF_OUT,20)
    write(TABF_OUT,30)

    if( fprint_icf ) then
        call pmf_utils_open(TABF_ICF,ftabficf,'R')

        write(TABF_ICF,10)
        write(TABF_ICF,20)
        write(TABF_ICF,30)

        call tabf_output_write_header_icf

        allocate(icf_cache(2*NumOfTABFCVs,fnstlim), stat= alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                '[TABF] Unable to allocate memory for icf_cache!')
        end if
        ! clear cache
        icf_cache(:,:) = 0.0d0
    end if

    return

10 format('#===============================================================================')
20 format('# Adaptive Biasing Force Method')
30 format('#===============================================================================')

end subroutine tabf_output_open

!===============================================================================
! Subroutine:  tabf_output_write_header
!===============================================================================

subroutine tabf_output_write_header

    use pmf_constants
    use pmf_dat
    use tabf_dat
    use pmf_cvs

    implicit none
    integer        :: i, off
    ! --------------------------------------------------------------------------

    write(TABF_OUT,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+NumOfTABFCVs
        write(TABF_OUT,15,advance='NO') i
    end do
    write(TABF_OUT,*)

    write(TABF_OUT,10,advance='NO') '#   NSTEP'
    do i=1,NumOfTABFCVs
        write(TABF_OUT,20,advance='NO') trim(TABFCVList(i)%cv%name)
    end do
    write(TABF_OUT,*)

    write(TABF_OUT,10,advance='NO') '#        '
    do i=1,NumOfTABFCVs
        write(TABF_OUT,25,advance='NO') '[' // trim(TABFCVList(i)%cv%get_ulabel()) // ']'
    end do
    write(TABF_OUT,*)

    write(TABF_OUT,10,advance='NO') '#--------'
    do i=1,NumOfTABFCVs
        write(TABF_OUT,30,advance='NO') '---------------'
    end do
    write(TABF_OUT,*)

    ! write header for IFC file if requested
    call tabf_output_write_header_icf

    return

10 format(A9)
15 format(1X,I15)
20 format(1X,A15)
25 format(1X,A15)
30 format(1X,A15)

end subroutine tabf_output_write_header

!===============================================================================
! Subroutine:  tabf_output_write
!===============================================================================

subroutine tabf_output_write

    use pmf_constants
    use pmf_dat
    use tabf_dat
    use pmf_cvs

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    if( fsample .le. 0 ) return ! output is written only of fsample > 0
    if( mod(fstep,fsample) .ne. 0 ) return

    write(TABF_OUT,10,advance='NO') fstep

    do i=1,NumOfTABFCVs
        write(TABF_OUT,20,advance='NO') &
                TABFCVList(i)%cv%get_rvalue(CVContext%CVsValues(TABFCVList(i)%cvindx))
    end do

    write(TABF_OUT,*)

    return

10 format(I9)
20 format(1X,F15.8)

end subroutine tabf_output_write

!===============================================================================
! Subroutine:  tabf_output_write_header_icf
!===============================================================================

subroutine tabf_output_write_header_icf

    use pmf_constants
    use pmf_dat
    use tabf_dat
    use pmf_cvs

    implicit none
    integer         :: i, off
    ! --------------------------------------------------------------------------

    if( .not. fprint_icf ) return

    write(TABF_ICF,10,advance='NO') '#       1'
    off = 1
    do i=off+1,off+2*NumOfTABFCVs,2
        write(TABF_ICF,15,advance='NO') i
        write(TABF_ICF,16,advance='NO') i+1
    end do
    write(TABF_ICF,*)

    write(TABF_ICF,10,advance='NO') '#   NSTEP'
    do i=1,NumOfTABFCVs
        write(TABF_ICF,20,advance='NO') trim(TABFCVList(i)%cv%name)
        write(TABF_ICF,25,advance='NO') 'ICF (' // trim(TABFCVList(i)%cv%name) // ')'
    end do
    write(TABF_ICF,*)

    write(TABF_ICF,10,advance='NO') '#        '
    do i=1,NumOfTABFCVs
        write(TABF_ICF,20,advance='NO') '[i.u.]'
        write(TABF_ICF,25,advance='NO') '[i.u.]'
    end do
    write(TABF_ICF,*)

    write(TABF_ICF,10,advance='NO') '#--------'
    do i=1,NumOfTABFCVs
        write(TABF_ICF,30,advance='NO') '---------------'
        write(TABF_ICF,35,advance='NO') '--------------------'
    end do
    write(TABF_ICF,*)

    write(TABF_ICF,40)

    return

10 format(A9)
15 format(1X,I15)
16 format(1X,I20)
20 format(1X,A15)
25 format(1X,A20)
30 format(1X,A15)
35 format(1X,A20)
40 format('# ICF accumulation is cached, data will be printed at the end at once')

end subroutine tabf_output_write_header_icf

!===============================================================================
! Subroutine:  tabf_output_write_icf
!===============================================================================

subroutine tabf_output_write_icf(cvs,gfx)

    use pmf_constants
    use pmf_dat
    use tabf_dat
    use pmf_cvs

    implicit none
    real(PMFDP)     :: cvs(:)
    real(PMFDP)     :: gfx(:)
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    if( .not. fprint_icf ) return

    ! write into memory cache
    if( (fstep .le. 0) .or. (fstep .gt. fnstlim) ) return ! out of cache size
    do i=1,NumOfTABFCVs
            icf_cache((i-1)*2+1,fstep) = cvs(i)
            icf_cache((i-1)*2+2,fstep) = gfx(i)
    end do

    return

end subroutine tabf_output_write_icf

!===============================================================================
! Subroutine:  tabf_output_dump_icf
!===============================================================================

subroutine tabf_output_dump_icf()

    use pmf_constants
    use pmf_dat
    use tabf_dat
    use pmf_cvs

    implicit none
    integer         :: i,j
    ! --------------------------------------------------------------------------

    if( .not. fprint_icf ) return

    ! 4 - do to numerical integration of ABF forces
    do j=4,fnstlim
        ! write into TABF_ICF
        write(TABF_ICF,10,advance='NO') j
        do i=1,NumOfTABFCVs
            write(TABF_ICF,20,advance='NO') icf_cache((i-1)*2+1,j)
            write(TABF_ICF,25,advance='NO') icf_cache((i-1)*2+2,j)
        end do
        write(TABF_ICF,*)
    end do

    return

10 format(I9)
20 format(1X,F15.8)
25 format(1X,F20.8)

end subroutine tabf_output_dump_icf

!===============================================================================
! Subroutine:  tabf_output_close
!===============================================================================

subroutine tabf_output_close

    use pmf_constants
    use pmf_dat
    use tabf_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(TABF_OUT)

    if( fprint_icf ) then
        call tabf_output_dump_icf
        close(TABF_ICF)
    end if

    return

end subroutine tabf_output_close

!===============================================================================

end module tabf_output
