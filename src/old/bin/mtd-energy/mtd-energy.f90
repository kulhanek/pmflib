! ===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! -------------------------------------------------------------------------------
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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
! ===============================================================================

program metadyn_energy

 use pmf_sizes
 use pmf_constants
 use pmf_core

 implicit none
 character(PMF_MAX_PATH)             :: frstfile        ! name of restart file
 integer                             :: fnitem          ! number of CVs
 double precision, allocatable       :: min_values(:)   ! minimal values of CVs
 double precision, allocatable       :: max_values(:)   ! maximal values of CVs
 integer, allocatable                :: grid_sizes(:)   ! number of grid points

 integer                             :: meta_time       ! use only specified meta time

 double precision, allocatable       :: point(:)        ! point on surface

 type history_buffer_rec
    type(history_buffer_rec),pointer    :: next_history_buffer
    integer                             :: length_of_buffer
    integer                             :: num_of_values
    double precision,pointer            :: values(:,:)
    double precision,pointer            :: widths(:,:)
    double precision,pointer            :: heights(:)
  end type history_buffer_rec

 integer,parameter                  :: MAX_BUFFER_SIZE = 1000
 type(history_buffer_rec),pointer   :: hill_history
 integer                            :: total_hills

 character(30)                      :: tmp
 integer                            :: num_args, liostat
 ! -----------------------------------------------------------------------------

 ! read input parameters
 num_args = iargc()
 if( num_args .ne. 2 ) then
    call libpmf_exit(PMF_OUT,1,'>>> ERROR: Two arguments are required!')
 end if

 call getarg(1, frstfile)
 call getarg(2, tmp)

 read(tmp,*, iostat=liostat) meta_time
 if( liostat .ne. 0 ) then
    call libpmf_exit(PMF_OUT,1,'>>> ERROR: Second argument has to be an integer number!')
 end if

 ! load dat ----------------------------
 call read_metadyn_restart

 if( meta_time .eq. 0 ) meta_time = total_hills

 ! write header
 write(*,10) fnitem
 write(*,20) total_hills
 write(*,30) meta_time

 if( meta_time .gt. total_hills ) then
    call libpmf_exit(PMF_OUT,1,'>>> ERROR: Meta-time is in future!')
 end if


 ! calculate profile--------------------
 call calculate_fes(1)

 stop

10 format('# Number of CVs         : ',I9)
20 format('# Total number of hills : ',I9)
30 format('# Used meta time        : ',I9)

contains

!===============================================================================
! ------------------------------------------------------------------------------
!===============================================================================

subroutine read_metadyn_restart

 use pmf_core

 implicit none
 character                          :: tmp
 integer                            :: i, j, alloc_failed,rstlevel,nstep
 integer                            :: fiostat
 type(history_buffer_rec),pointer   :: current_history
 ! -----------------------------------------------------------------------------

 ! open restart file -----------------------------------------------------------
 call libpmf_open(MED_RST,frstfile,'O')

 ! read header -----------------------------------------------------------------
 read(MED_RST,*) tmp, fnitem

 if( fnitem .le. 0 ) then
    write(PMF_OUT,'(/,A)') '>>> ERROR: Number of CVs is zero or negative number!'
    call libpmf_exit(PMF_OUT,1)
 end if

 ! allocate data ---------------------------------------------------------------
 allocate(  min_values(fnitem), &
            max_values(fnitem), &
            grid_sizes(fnitem), &
            point(fnitem), &
            stat=alloc_failed)

 if( alloc_failed .ne. 0 ) then
    write(PMF_OUT,*) '>>> ERROR: Unable to allocate memory for description of CVs!'
    call libpmf_exit(PMF_OUT, 1)
 endif

 min_values(:) = 0.0
 max_values(:) = 0.0
 grid_sizes(:) = 0.0
 point(:) = 0.0

 ! read rest of header ---------------------------------------------------------
 do i=1,fnitem
    read(MED_RST,*) tmp, j, min_values(i), max_values(i), grid_sizes(i)
    if( j .ne. i ) then
        write(PMF_OUT,'(/,A,I3,A,I3)') '>>> ERROR: Illegal CV record! Found ',j,' but should be ',i
        call libpmf_exit(PMF_OUT, 1)
    end if
 end do

 ! now process hills ----------------------------------------------------------

 call allocate_buffer(hill_history)

 total_hills = 0
 current_history => hill_history

 do while(.true.)
    if( current_history%num_of_values .eq. current_history%length_of_buffer ) then
        call allocate_buffer(current_history%next_history_buffer)
        current_history => current_history%next_history_buffer
        cycle
    end if
    ! this is just fusion
    if( current_history%num_of_values .gt. current_history%length_of_buffer ) then
        call libpmf_exit(PMF_OUT, 1,'>>> ERROR: Buffer length is inconsistent (should not happen)!')
    end if

    current_history%num_of_values = current_history%num_of_values + 1
    read(MED_RST,'(2X,I5,2X,I10,1X,F10.3,1X)',ADVANCE='NO',iostat=fiostat) &
                rstlevel, nstep, current_history%heights(current_history%num_of_values)
    if( fiostat .lt. 0 ) then
        current_history%num_of_values = current_history%num_of_values - 1
        exit       ! end of file
    end if
    if( fiostat .ne. 0 ) then
        call libpmf_exit(PMF_OUT, 1,'>>> ERROR: Unable to read line from history list!')
    end if
    do i=1,fnitem
        read(MED_RST,'(F10.4,1X,F10.4,1X)',ADVANCE='NO',iostat=fiostat) &
                current_history%values(current_history%num_of_values,i), &
                current_history%widths(current_history%num_of_values,i)
    end do
    if( fiostat .ne. 0 ) then
        call libpmf_exit(PMF_OUT, 1,'>>> ERROR: Unable to read line from history list!')
    end if
    read(MED_RST,*,iostat=fiostat)
    total_hills = total_hills + 1
 end do

 ! close restart file ---------------------------------------------------------
 close(MED_RST)

 return

end subroutine read_metadyn_restart

!===============================================================================
! ------------------------------------------------------------------------------
!===============================================================================

subroutine allocate_buffer(current_history)

 use pmf_core

 implicit none
 type(history_buffer_rec),pointer   :: current_history
 integer                            :: alloc_failed
 ! -----------------------------------------------------------------------------

 ! create first buffer --------------------------------------------------------
 ! history node
 allocate(current_history,stat=alloc_failed)

 if( alloc_failed .ne. 0 ) then
   call libpmf_exit(PMF_OUT, 1,'>>> ERROR: Unable to allocate memory for history node!')
 endif

 ! node items
 allocate(current_history%values(MAX_BUFFER_SIZE,fnitem), &
          current_history%widths(MAX_BUFFER_SIZE,fnitem), &
          current_history%heights(MAX_BUFFER_SIZE), &
          stat=alloc_failed)

 if( alloc_failed .ne. 0 ) then
   call libpmf_exit(PMF_OUT, 1,'>>> ERROR: Unable to allocate memory for history items!')
 endif

 nullify(current_history%next_history_buffer)
 current_history%length_of_buffer = MAX_BUFFER_SIZE
 current_history%num_of_values = 0

 return

end subroutine allocate_buffer

!===============================================================================
! ------------------------------------------------------------------------------
!===============================================================================

recursive subroutine calculate_fes(coordinate)

 implicit none
 integer            :: coordinate
 integer            :: i
 double precision   :: sp, value
 ! -----------------------------------------------------------------------------

 if( coordinate .gt. fnitem ) then
    ! calculate value
    value = calculate_value()
    ! print point position and value
    write(*,*) (point(i),i=1,fnitem), -value
    return
 end if

 ! cycle through variable
 point(coordinate) = min_values(coordinate)
 sp = (max_values(coordinate) - min_values(coordinate)) / grid_sizes(coordinate)

 do while( point(coordinate) .le. max_values(coordinate) )
    call calculate_fes(coordinate+1)
    point(coordinate) = point(coordinate) + sp
 end do

 ! write block delimiter
 write(*,*)

 return

end subroutine calculate_fes

!===============================================================================
! ------------------------------------------------------------------------------
!===============================================================================

function calculate_value()

 implicit none
 double precision                   :: calculate_value
 double precision                   :: fexparg
 integer                            :: j,k,processed
 type(history_buffer_rec),pointer   :: current_history
 ! -----------------------------------------------------------------------------

 calculate_value = 0.0d0

 ! go through history
 processed = 0
 current_history => hill_history

 do while( associated(current_history) )
    !for every item in buffer
    do j = 1,current_history%num_of_values
        fexparg = 0.0d0
        do k = 1,fnitem
            fexparg = fexparg + (point(k) - current_history%values(j,k))**2 / &
                                (2.0d0 * current_history%widths(j,k)**2)
        end do
        calculate_value = calculate_value + current_history%heights(j)*exp(-fexparg)
        processed = processed + 1
        if( processed .ge. meta_time ) return
    end do
    current_history => current_history%next_history_buffer
 end do

 return

end function calculate_value

!===============================================================================
! ------------------------------------------------------------------------------
!===============================================================================

end program metadyn_energy