!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module rst_output

implicit none
contains

!===============================================================================
! Subroutine:  rst_output_open
!===============================================================================

subroutine rst_output_open

    use pmf_constants
    use pmf_dat
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    ! open output file
    call pmf_utils_open(RST_OUT,frstout,'R')

    write(RST_OUT,10)
    write(RST_OUT,20)
    write(RST_OUT,30)

    return

    10 format('#===============================================================================')
    20 format('# Restrained Dynamics')
    30 format('#===============================================================================')

end subroutine rst_output_open

!===============================================================================
! Subroutine:  rst_output_write_header
!===============================================================================

subroutine rst_output_write_header

    use pmf_constants
    use pmf_dat
    use rst_dat
    use pmf_unit
    use pmf_cvs

    implicit none
    integer        :: i, off, iend
    ! --------------------------------------------------------------------------

    write(RST_OUT,10,advance='NO') '#          '
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,20,advance='NO') trim(RSTCVList(i)%cv%name)
            case(1)
                write(RST_OUT,25,advance='NO') trim(RSTCVList(i)%cv%name)
            case(2)
                write(RST_OUT,26,advance='NO') trim(RSTCVList(i)%cv%name)
            case(3)
                write(RST_OUT,27,advance='NO') trim(RSTCVList(i)%cv%name)
        end select
    end do
    write(RST_OUT,*)

    write(RST_OUT,10,advance='NO') '#          '
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,21,advance='NO') '---------------'
            case(1)
                write(RST_OUT,20,advance='NO') '---------------'
                write(RST_OUT,21,advance='NO') '----------------'
            case(2)
                write(RST_OUT,20,advance='NO') '---------------'
                write(RST_OUT,21,advance='NO') '----------------'
                write(RST_OUT,21,advance='NO') '----------------'
            case(3)
                write(RST_OUT,20,advance='NO') '---------------'
                write(RST_OUT,21,advance='NO') '----------------'
                write(RST_OUT,21,advance='NO') '----------------'
                write(RST_OUT,21,advance='NO') '----------------'
        end select
    end do
    write(RST_OUT,*)

    write(RST_OUT,10,advance='NO') '#          '
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,20,advance='NO') 'Value'
            case(1)
                write(RST_OUT,20,advance='NO') 'Value'
                write(RST_OUT,20,advance='NO') 'Energy'
            case(2)
                write(RST_OUT,20,advance='NO') 'Value'
                write(RST_OUT,20,advance='NO') 'Energy'
                write(RST_OUT,20,advance='NO') 'Deviation'
            case(3)
                write(RST_OUT,20,advance='NO') 'Value'
                write(RST_OUT,20,advance='NO') 'Energy'
                write(RST_OUT,20,advance='NO') 'Deviation'
                write(RST_OUT,20,advance='NO') 'Target value'
        end select
    end do
    write(RST_OUT,*)

    write(RST_OUT,10,advance='NO') '#  NSTEP  F'
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
            case(1)
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
                write(RST_OUT,40,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
            case(2)
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
                write(RST_OUT,40,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
            case(3)
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
                write(RST_OUT,40,advance='NO') '['//trim(pmf_unit_label(EnergyUnit))//']'
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
                write(RST_OUT,40,advance='NO') '['//trim(RSTCVList(i)%cv%get_ulabel())//']'
        end select
    end do
    write(RST_OUT,*)

    write(RST_OUT,10,advance='NO') '#-------- -'
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,30,advance='NO') '---------------'
            case(1)
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
            case(2)
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
            case(3)
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
        end select
    end do
    write(RST_OUT,*)

    write(RST_OUT,10,advance='NO') '#       1 2'
    off = 2
    select case(fplevel)
        case(0)
            iend = NumOfRSTCVs
        case(1)
            iend = 2*NumOfRSTCVs
        case(2)
            iend = 3*NumOfRSTCVs
        case(3)
            iend = 4*NumOfRSTCVs
    end select
    do i=off+1,off+iend
        write(RST_OUT,15,advance='NO') i
    end do
    write(RST_OUT,*)

    write(RST_OUT,10,advance='NO') '#--------'
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,30,advance='NO') '---------------'
            case(1)
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
            case(2)
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
            case(3)
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
                write(RST_OUT,30,advance='NO') '---------------'
        end select
    end do
    write(RST_OUT,*)

    return

10 format(A11)
15 format(1X,I15)
20 format(1X,A15)
21 format(A16)
25 format(1x,A31)
26 format(1x,A47)
27 format(1x,A63)
30 format(1X,A15)
40 format(1X,A15)

end subroutine rst_output_write_header

!===============================================================================
! Subroutine:  rst_output_write_data
!===============================================================================

subroutine rst_output_write_data

    use pmf_cvs
    use rst_dat
    use pmf_unit
    use pmf_dat

    implicit none
    integer             :: i
    character(len=1)    :: flag
    ! --------------------------------------------------------------------------

    flag = 'N'

    if( fwarnlevel .ge. 0.0d0 ) then
        do i=1,NumOfRSTCVs
            if( RSTCVList(i)%energy .gt. fwarnlevel ) flag = 'W'
        end do
    end if

    if( flag .eq. 'N' ) then
        if( fsample .le. 0 ) return ! output is written only of fsample > 0
        if( mod(fstep,fsample) .ne. 0 ) return
    end if

    ! print data to output --------------------------
    write(RST_OUT,200,ADVANCE='NO') fstep, flag
    do i=1,NumOfRSTCVs
        select case(fplevel)
            case(0)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(CVContext%CVsValues(RSTCVList(i)%cvindx))
            case(1)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(CVContext%CVsValues(RSTCVList(i)%cvindx))
                write(RST_OUT,210,ADVANCE='NO') pmf_unit_get_rvalue(EnergyUnit,RSTCVList(i)%energy)
            case(2)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(CVContext%CVsValues(RSTCVList(i)%cvindx))
                write(RST_OUT,210,ADVANCE='NO') pmf_unit_get_rvalue(EnergyUnit,RSTCVList(i)%energy)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(RSTCVList(i)%deviation)
            case(3)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(CVContext%CVsValues(RSTCVList(i)%cvindx))
                write(RST_OUT,210,ADVANCE='NO') pmf_unit_get_rvalue(EnergyUnit,RSTCVList(i)%energy)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(RSTCVList(i)%deviation)
                write(RST_OUT,210,ADVANCE='NO') RSTCVList(i)%cv%get_rvalue(RSTCVList(i)%target_value)
        end select
    end do
    write(RST_OUT,*)

    return

200 format(I9,1X,A1)
210 format(1X,E15.8)

end subroutine rst_output_write_data

!===============================================================================
! Subroutine:  rst_output_close
!===============================================================================

subroutine rst_output_close

    use pmf_constants

    implicit none
    ! --------------------------------------------------------------------------

    close(RST_OUT)

    return

end subroutine rst_output_close

!===============================================================================

end module rst_output
