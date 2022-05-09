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
    integer         :: i
    real(PMFDP)     :: mlogml, slogml, mlogpl, slogpl
    ! --------------------------------------------------------------------------

    if( gpr_calc_logxx ) then
        write(ABF_OUT,*)
        write(ABF_OUT,5)
        write(ABF_OUT,10)
        write(ABF_OUT,20)
        if( gpr_icf_enabled ) then
            do i=1,NumOfABFCVs
                mlogml = gpr_icf_mlogml(i)
                slogml = sqrt(gpr_icf_m2logml(i)/gpr_icf_nlogxx)
                mlogpl = gpr_icf_mlogpl(i)
                slogpl = sqrt(gpr_icf_m2logpl(i)/gpr_icf_nlogxx)
                write(ABF_OUT,30) i, trim(CVList(ABFCVList(i)%cvindx)%cv%name), mlogml, slogml, mlogpl, slogpl

            end do
        end if
        if( gpr_ene_enabled ) then
            mlogml = gpr_ene_mlogml
            slogml = sqrt(gpr_ene_m2logml/gpr_ene_nlogxx)
            mlogpl = gpr_ene_mlogpl
            slogpl = sqrt(gpr_ene_m2logpl/gpr_ene_nlogxx)
            write(ABF_OUT,30) i, 'etot', mlogml, slogml, mlogpl, slogpl
        end if
    end if

    close(ABF_OUT)

    return

  5 format('# GPR statistics')
 10 format('# N  CV name     <logML>      s(logML)     <logPL>      s(logPL)  ')
 20 format('#-- ---------- ------------ ------------ ------------ ------------')
 30 format(I3,1X,A10,1X,E12.5,1X,E12.5,1X,E12.5,1X,E12.5)

end subroutine abf_output_close

!===============================================================================

end module abf_output
