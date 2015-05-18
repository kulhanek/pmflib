!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module mtd_restart

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_restart_open_and_read
!===============================================================================

subroutine mtd_restart_open_and_read

 use pmf_utils
 use pmf_dat
 use mtd_dat
 use mtd_history

 implicit none
 character(5)               :: mtd_trac
 character(PMF_MAX_TYPE)    :: cvs_trac
 character(PMF_MAX_CV_NAME) :: name_trac
 character(2)               :: ver_trac
 integer                    :: oldnitem,i,id,nbins,total_hills,nstep
 integer                    :: fiostat
 double precision           :: min_value,max_value
 type(MTDHistType),pointer  :: current_history
 ! -----------------------------------------------------------------------------

 ! test if restart file exists
 if( frestart .and. .not. pmf_utils_fexist(fmtdrst) ) then
    frestart = .false.
    write(MTD_OUT,30) trim(fmtdrst)
 end if

 if( frestart ) then
    ! open restart file -----------------------------------------------------------
    call pmf_utils_open(MTD_RST,fmtdrst,'O')

    ! read header -----------------------------------------------------------------
    read(MTD_RST,10) mtd_trac, ver_trac, oldnitem

    if( mtd_trac .ne. 'MTD' ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Attempt to read non-MTD restart file (MTD trac is missing)!')
    end if

    if( ver_trac .ne. 'V3' ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Attempt to read non-MTD restart file (not V3 version)!')
    end if

    if(oldnitem .ne. NumOfMTDCVs) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Previous and present number of CVs is different!')
    end if

    ! read rest of header ---------------------------------------------------------
    do i=1,NumOfMTDCVs
        read(MTD_RST,20) id, cvs_trac, min_value, max_value, nbins
        if( i .ne. id ) then
            call pmf_utils_exit(PMF_OUT, 1,'[MTD] Illegal CV record!')
        end if
        if( trim(adjustl(cvs_trac)) .ne. trim(MTDCVList(i)%cv%ctype) ) then
            write(PMF_OUT,*) '[MTD] CV type = [',trim(adjustl(cvs_trac)), &
                             '] should be [',trim(MTDCVList(i)%cv%ctype),']'
            call pmf_utils_exit(PMF_OUT,1,'[MTD] CV type was redefined in MTD restart file!')
        end if
        read(MTD_RST,25) id, name_trac
        if( i .ne. id ) then
            call pmf_utils_exit(PMF_OUT, 1,'[MTD] Illegal CV record!')
        end if
        if( trim(adjustl(name_trac)) .ne. trim(MTDCVList(i)%cv%name) ) then
            write(PMF_OUT,*) '[MTD] CV name = [',trim(adjustl(cvs_trac)),'] should be [',trim(MTDCVList(i)%cv%name),']'
            call pmf_utils_exit(PMF_OUT,1,'[MTD] CV name was redefined in MTD restart file!')
        end if
    end do

    ! now process hills ----------------------------------------------------------
    if( .not. associated(hill_history) ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] hill_history is not allocated!')
    end if

    total_hills = 0
    current_history => hill_history

    do while(.true.)
        if( current_history%nrst_of_values .eq. current_history%length_of_buffer ) then
            call mtd_history_allocate_buffer(current_history%next_history_buffer,fbuffersize)
            current_history => current_history%next_history_buffer
            cycle
        end if
        ! this is just fusion
        if( current_history%nrst_of_values .gt. current_history%length_of_buffer ) then
            call pmf_utils_exit(PMF_OUT, 1,'[MTD] Buffer length is inconsistent (should not happen)!')
        end if

        current_history%nrst_of_values = current_history%nrst_of_values + 1
        read(MTD_RST,'(2X,I5,2X,I10,1X,F10.3,1X)',ADVANCE='NO',iostat=fiostat) &
                    frstlevel, nstep, current_history%heights(current_history%nrst_of_values)
        if( fiostat .lt. 0 ) then
            current_history%nrst_of_values = current_history%nrst_of_values - 1
            exit       ! end of file
        end if
        if( fiostat .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to read line from history list!')
        end if
        do i=1,NumOfMTDCVs
            read(MTD_RST,'(F10.4,1X,F10.4,1X)',ADVANCE='NO',iostat=fiostat) &
                    current_history%values(current_history%nrst_of_values,i), &
                    current_history%widths(current_history%nrst_of_values,i)
        end do
        if( fiostat .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to read line from history list!')
        end if
        read(MTD_RST,*,iostat=fiostat)
        total_hills = total_hills + 1
    end do

    write(MTD_OUT,40) total_hills

 else
    ! open output file ------------------------------------------------------------
    call pmf_utils_open(MTD_RST,fmtdrst,'R')

    ! write header ----------------------------------------------------------------
    write(MTD_RST,10) 'MTD', 'V3', NumOfMTDCVs

    do i=1,NumOfMTDCVs
        write(MTD_RST,20) i,trim(MTDCVList(i)%cv%ctype), &
                      MTDCVList(i)%min_value, &
                      MTDCVList(i)%max_value, &
                      MTDCVList(i)%nbins
        write(MTD_RST,25) i,trim(MTDCVList(i)%cv%name)
    end do

    write(MTD_OUT,50)
 end if

 frstlevel = frstlevel + 1  ! increase restart level

 return

 10  format(A3,1X,A2,1X,I2)
 20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
 25  format(I2,1X,A55)

 30 format('# WARNING: frestart = on, but file (',A,') does not exist! => frestart = off')
 40 format('# RST: Number of hills read from restart file: ',I9)
 50 format('# RST: frestart = off')

 return

end subroutine mtd_restart_open_and_read

!===============================================================================
! Subroutine:  mtd_restart_close
!===============================================================================

subroutine mtd_restart_close

    use pmf_constants
    use mtd_dat

    implicit none
    ! --------------------------------------------------------------------------

    close(MTD_RST)

    return

end subroutine mtd_restart_close

!===============================================================================

end module mtd_restart
