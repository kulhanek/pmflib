!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module mtd_history

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine: mtd_history_print_buffer
!===============================================================================

subroutine mtd_history_print_buffer(printed_buffer,level,start)

 use pmf_dat
 use mtd_dat

 implicit none
 type(MTDHistType),pointer          :: printed_buffer
 integer                            :: level
 integer                            :: start
 ! -----------------------------------------------
 integer                            :: i,j
 ! -----------------------------------------------------------------------------

 if( .not. associated(printed_buffer) ) return ! nothing to print

 ! write info to restart file -------------------------------------------------
 do i=1,printed_buffer%nrst_of_values
    write(MTD_RST,10,ADVANCE='NO') level, start+i-1, printed_buffer%heights(i)
    do j=1, NumOfMTDCVs
        write(MTD_RST,20,ADVANCE='NO') printed_buffer%values(i,j),printed_buffer%widths(i,j)
    end do
    write(MTD_RST,*)
 end do

 return

 10 format(2X,I5,2X,I10,1X,F10.3,1X)
 20 format(F10.4,1X,F10.4,1X)

end subroutine mtd_history_print_buffer

!===============================================================================
! Subroutine: mtd_history_add_new_hill
!===============================================================================

subroutine mtd_history_add_new_hill(mtd_energy)

 use pmf_dat
 use mtd_dat
 use mtd_grid

 implicit none
 real(PMFDP), intent(in)            :: mtd_energy
 integer                            :: i,ci
 type(MTDHistType),pointer          :: last_buffer
 real(PMFDP)                        :: hill_height
 real(PMFDP)                        :: values(NumOfMTDCVs), widths(NumOfMTDCVs)
 ! -----------------------------------------------------------------------------

 ! conventional MTD
 hill_height = fheight

 ! well-tempered MTD
 if( fmetatemp .ne. 0 ) then
     ! what to change?
     select case(fmetavary)
        case(0)
        	! no metavary
        case(1)
        	! vary the height
            hill_height = fheight * exp( -mtd_energy/PMF_Rgas/fmetatemp )
        case(2)
        	! vary the step
            meta_next_fstep = fstep + nint( fmetastep * exp( mtd_energy/PMF_Rgas/fmetatemp ) )
     end select
 end if

 if (fmode .ne. 2) then
    ! find last active buffer ----------------------------------------------------
    last_buffer => hill_history
    do while( associated(last_buffer) )
       if( .not. associated(last_buffer%next_history_buffer) ) exit
       last_buffer => last_buffer%next_history_buffer
    end do

    !do we need new buffer --------------------------------------------------------
    if( last_buffer%nrst_of_values .ge. last_buffer%length_of_buffer ) then
       call mtd_history_allocate_buffer(last_buffer%next_history_buffer,fbuffersize)
       last_buffer => last_buffer%next_history_buffer
    end if

    ! add hill --------------------------------------------------------------------
    last_buffer%nrst_of_values = last_buffer%nrst_of_values + 1

    if( last_buffer%nrst_of_values .eq. last_buffer%length_of_buffer ) then
       fserverupdate = .true.  ! buffer will be ready - we should exchange data
    end if

    do i = 1,NumOfMTDCVs
       ci = MTDCVList(i)%cvindx
       last_buffer%values(last_buffer%nrst_of_values,i) = CVContext%CVsValues(ci)
       last_buffer%widths(last_buffer%nrst_of_values,i) = MTDCVList(i)%width
    end do
    last_buffer%heights(last_buffer%nrst_of_values) = hill_height
    ! record fstep for delaying and hill scaling
    last_buffer%steps(last_buffer%nrst_of_values) = fstep

    ! write info to restart file -------------------------------------------------
    if( fserver_enabled ) then
       ! as "time" we use meta_step
       write(MTD_RST,10,ADVANCE='NO') fclient_id, meta_step, &
                                      last_buffer%heights(last_buffer%nrst_of_values)
    else
       ! as time we use real time step
       write(MTD_RST,10,ADVANCE='NO') frstlevel, fstep, &
                                      last_buffer%heights(last_buffer%nrst_of_values)
    end if

    do i=1, NumOfMTDCVs
       write(MTD_RST,20,ADVANCE='NO') last_buffer%values(last_buffer%nrst_of_values,i), &
                                      last_buffer%widths(last_buffer%nrst_of_values,i)
    end do
    write(MTD_RST,*)

    select case(fextout)
       case(1,2)
           ! write CVs
           write(MTD_CVS,30,ADVANCE='NO') fstep+cvs_starts
           do i=1, NumOfMTDCVs
               write(MTD_CVS,40,ADVANCE='NO') last_buffer%values(last_buffer%nrst_of_values,i)
           end do
           write(MTD_CVS,*)

           ! write hills
           write(MTD_HILLS,30,ADVANCE='NO') fstep+cvs_starts
           do i=1, NumOfMTDCVs
               write(MTD_HILLS,40,ADVANCE='NO') last_buffer%values(last_buffer%nrst_of_values,i)
           end do
           do i=1, NumOfMTDCVs
               write(MTD_HILLS,40,ADVANCE='NO') last_buffer%widths(last_buffer%nrst_of_values,i)
           end do
           write(MTD_HILLS,50,ADVANCE='NO') last_buffer%heights(last_buffer%nrst_of_values)
           write(MTD_HILLS,*)
    end select
 else
    do i = 1,NumOfMTDCVs
       ci = MTDCVList(i)%cvindx
       values(i) = CVContext%CVsValues(ci)
       widths(i) = MTDCVList(i)%width
    end do

    call mtd_grid_add_data(values, hill_height, widths) 

    ! write info to restart file -------------------------------------------------
    if( fserver_enabled ) then
       ! as "time" we use meta_step
       write(MTD_RST,10,ADVANCE='NO') fclient_id, meta_step, &
                                      hill_height
    else
       ! as time we use real time step
       write(MTD_RST,10,ADVANCE='NO') frstlevel, fstep, &
                                      hill_height
    end if

    do i=1, NumOfMTDCVs
       write(MTD_RST,20,ADVANCE='NO') values(i), &
                                      widths(i)
    end do
    write(MTD_RST,*)

    select case(fextout)
       case(1,2)
           ! write CVs
           write(MTD_CVS,30,ADVANCE='NO') fstep+cvs_starts
           do i=1, NumOfMTDCVs
               write(MTD_CVS,40,ADVANCE='NO') values(i)
           end do
           write(MTD_CVS,*)

           ! write hills
           write(MTD_HILLS,30,ADVANCE='NO') fstep+cvs_starts
           do i=1, NumOfMTDCVs
               write(MTD_HILLS,40,ADVANCE='NO') values(i)
           end do
           do i=1, NumOfMTDCVs
               write(MTD_HILLS,40,ADVANCE='NO') widths(i)
           end do
           write(MTD_HILLS,50,ADVANCE='NO') hill_height
           write(MTD_HILLS,*)
    end select
 end if

 return

 10 format(2X,I5,2X,I10,1X,F10.3,1X)
 20 format(F10.4,1X,F10.4,1X)
 30 format(2X,I10,1X)
 40 format(F18.11,1X)
 50 format(F18.11)

end subroutine mtd_history_add_new_hill

!===============================================================================
! Subroutine: mtd_history_allocate_buffer
!===============================================================================

subroutine mtd_history_allocate_buffer(current_history,buffer_size)

 use pmf_utils
 use mtd_dat

 implicit none
 type(MTDHistType),pointer          :: current_history
 integer                            :: buffer_size
 ! -----------------------------------------------
 integer                            :: alloc_failed
 ! -----------------------------------------------------------------------------

 ! create first buffer ---------------------------------------------------------
 ! history node
 allocate(current_history,stat=alloc_failed)

 if( alloc_failed .ne. 0 ) then
   call pmf_utils_exit(PMF_OUT, 1,'>>> ERROR:Unable to allocate memory for history node!')
 endif

 ! node items
 allocate(current_history%values(buffer_size,NumOfMTDCVs), &
          current_history%widths(buffer_size,NumOfMTDCVs), &
          current_history%heights(buffer_size), &
          current_history%steps(buffer_size), &
          stat=alloc_failed)

 if( alloc_failed .ne. 0 ) then
   call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for history items!')
 endif

 nullify(current_history%next_history_buffer)
 current_history%length_of_buffer = fbuffersize
 current_history%nrst_of_values = 0

 return

end subroutine mtd_history_allocate_buffer

!===============================================================================

end module mtd_history
