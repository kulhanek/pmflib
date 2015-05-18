!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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

module mtd_client

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interfaces
!===============================================================================

interface
    ! set number of coordinates
    subroutine cpmf_mtd_client_set_header(ret_st,nitems)
        integer         :: ret_st
        integer         :: nitems
    end subroutine cpmf_mtd_client_set_header

    ! set coordinate data
    subroutine cpmf_mtd_client_set_coord(ret_st,id,ctype,name,min_value,max_value,nbins)
        integer         :: ret_st
        integer         :: id
        character(*)    :: ctype
        character(*)    :: name
        real(8)         :: min_value
        real(8)         :: max_value
        integer         :: nbins
    end subroutine cpmf_mtd_client_set_coord

    ! register client on server
    subroutine cpmf_mtd_client_reg_by_name(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_mtd_client_reg_by_name

    ! register client on server
    subroutine cpmf_mtd_client_reg_by_key(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_mtd_client_reg_by_key

    ! transfer initial data from server to client and store them in cpmd buffer
    subroutine cpmf_mtd_client_initial_data(ret_st)
        integer         :: ret_st
    end subroutine cpmf_mtd_client_initial_data

    ! exchange data with server
    subroutine cpmf_mtd_client_exchange_data(ret_st)
        integer         :: ret_st
    end subroutine cpmf_mtd_client_exchange_data

    ! get info about cpmd buffer with index id
    subroutine cpmf_mtd_client_get_buffer_info(ret_st,id,level,start,nrst_of_hills)
        integer         :: ret_st
        integer         :: id
        integer         :: level
        integer         :: start
        integer         :: nrst_of_hills
    end subroutine cpmf_mtd_client_get_buffer_info

    ! add new buffer to cpmf buffer list
    subroutine cpmf_mtd_client_add_buffer_data(ret_st,level,start,nrst_of_hills,history_buffer)
        integer         :: ret_st
        integer         :: level
        integer         :: start
        integer         :: nrst_of_hills
        real(8)         :: history_buffer(*)
    end subroutine cpmf_mtd_client_add_buffer_data

    ! get cpmf buffer data with index id to history_buffer
    subroutine cpmf_mtd_client_get_buffer_data(ret_st,id,history_buffer)
        integer         :: ret_st
        integer         :: id
        real(8)         :: history_buffer(*)
    end subroutine cpmf_mtd_client_get_buffer_data

    ! clear cpmd buffer list
    subroutine cpmf_mtd_client_clear_buf_list()
    end subroutine cpmf_mtd_client_clear_buf_list

    ! unregister client on server
    subroutine cpmf_mtd_client_unregister()
    end subroutine cpmf_mtd_client_unregister
end interface

contains

!===============================================================================
! Subroutine:  mtd_client_register
!===============================================================================

subroutine mtd_client_register

 use pmf_dat
 use mtd_dat
 use pmf_utils
 use pmf_timers
 use pmf_unit

 implicit none
#ifdef PMFLIB_NETWORK
 integer        :: ret_st, i
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'Multiple-walkers MTD method','=')

 if( fserver_key_enabled ) then
    write(PMF_OUT,15) trim(fserverkey)
 else
    write(PMF_OUT,10) trim(fserver)
 end if
 write(PMF_OUT,*)
 write(PMF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! register coordinates
    call cpmf_mtd_client_set_header(ret_st,NumOfMTDCVs)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1)
    end if

    do i=1,NumOfMTDCVs
        call cpmf_mtd_client_set_coord(ret_st, &
                        i, &
                        MTDCVList(i)%cv%ctype, &
                        MTDCVList(i)%cv%name, &
                        MTDCVList(i)%min_value, &
                        MTDCVList(i)%max_value, &
                        MTDCVList(i)%nbins)

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1)
        end if
    end do

    ! register client
    if( fserver_key_enabled ) then
        call cpmf_mtd_client_reg_by_key(fserverkey,fserver,fclient_id)
    else
        call cpmf_mtd_client_reg_by_name(fserver,fpassword,fclient_id)
    end if

    write(MTD_OUT,*)

    if( fclient_id .le. 0 ) then
        fserver_enabled = .false.
        write(PMF_OUT,30)
        write(MTD_OUT,50) trim(fserver)

        call pmf_utils_exit(PMF_OUT,1)
    else
        write(PMF_OUT,40) fclient_id
        write(MTD_OUT,60) trim(fserver)
        write(MTD_OUT,70) fclient_id
    end if

#endif

 call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)

 return

 10 format(' MTD Server          : ', A)
 15 format(' MTD Server Key file : ', A)
 20 format(' Registering client on server, please wait .... ')

 30 format(' Registration FAILED!')
 40 format(' Registration SUCCESSFULL! (Client ID: ',I6,')')

 50 format('# [MTD-CLIENT] Registration to server ',A,' failed!')
 60 format('# [MTD-CLIENT] Registration to server ',A,' successful.')
 70 format('# [MTD-CLIENT] Client ID: ',I6)

end subroutine mtd_client_register

!===============================================================================
! Subroutine:  mtd_client_get_initial_data
!===============================================================================

subroutine mtd_client_get_initial_data

 use pmf_utils
 use mtd_dat
 use mtd_history
 use pmf_timers

 implicit none
#ifdef PMFLIB_NETWORK
 integer                    :: ret_st = 0
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

#ifdef PMFLIB_NETWORK
    call cpmf_mtd_client_initial_data(ret_st)

    if( ret_st .ne. 0 ) then
        fserver_enabled = .false.   ! disable server
        write(MTD_OUT,20)
        write(PMF_OUT,*)
        call pmf_utils_exit(PMF_OUT,1)
    end if

    ! add initial data
    call mtd_client_get_cpmf_buffers(hill_history)

    ! did we receive any data?
    if( .not. associated(hill_history) )then
        ! if not - then allocate basicbuffer
        call mtd_history_allocate_buffer(hill_history,fbuffersize)
    end if

    ! clear buffer list
    call cpmf_mtd_client_clear_buf_list

    write(MTD_OUT,10)
    write(PMF_OUT,*)
#endif

 call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)

 return

 10 format('# [MTD-CLIENT] Initial data from server were uploaded')
 20 format('# [MTD-CLIENT] Unable to get initial data from server - disabling method')

end subroutine mtd_client_get_initial_data

!===============================================================================
! Subroutine:  mtd_client_exchange_data
!===============================================================================

subroutine mtd_client_exchange_data

 use pmf_utils
 use pmf_dat
 use mtd_dat
 use pmf_timers
 use pmf_exit

 implicit none
#ifdef PMFLIB_NETWORK
 integer                    :: ret_st = 0
 type(MTDHistType),pointer  :: last_buffer
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 ! do we need to exchange data?
 if( .not. fserverupdate ) return
 fserverupdate = .false.

 call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

 write(MTD_OUT,10) fstep

#ifdef PMFLIB_NETWORK
    ! find last buffer
    last_buffer => hill_history
    do while( associated(last_buffer) )
        if( .not. associated(last_buffer%next_history_buffer) ) exit
        last_buffer => last_buffer%next_history_buffer
    end do

    ! add buffer to cpmd list
    call mtd_client_add_cpmf_buffer(last_buffer)

    ! now exchange data
    call cpmf_mtd_client_exchange_data(ret_st)

    if( ret_st .ne. 0 ) then
        fserver_enabled = .false.   ! disable server in the case of error
        write(MTD_OUT,20)
        write(PMF_OUT,*)
        write(PMF_OUT,20)
        call pmf_exit_mdloop(PMF_OUT,1)
        return
    end if

    ! add receive buffers
    call mtd_client_get_cpmf_buffers(last_buffer)

    ! clear buffer list
    call cpmf_mtd_client_clear_buf_list
#endif

 call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)

 return

 10 format('# [MTD-CLIENT] ',I9,' Exchanging data with server')
 20 format('# [MTD-CLIENT] Unable to exchange data with server - disabling method')

end subroutine mtd_client_exchange_data

!===============================================================================
! Subroutine:  mtd_client_add_cpmf_buffer
!===============================================================================

subroutine mtd_client_add_cpmf_buffer(last_buffer)

 use pmf_utils
 use mtd_dat

 implicit none
 type(MTDHistType),pointer          :: last_buffer
 ! -----------------------------------------------
#ifdef PMFLIB_NETWORK
 integer                            :: ret_st
 integer                            :: nrst_of_hills,nrst_of_records
 integer                            :: alloc_failed,i,j,k
 real(PMFDP),allocatable            :: loc_hbuffer(:)
#endif
 !-----------------------------------------------------------------------------

 if( .not. associated(last_buffer) ) return ! nothing to set

#ifdef PMFLIB_NETWORK

    nrst_of_hills = last_buffer%nrst_of_values
    nrst_of_records = nrst_of_hills*(1 + 2*NumOfMTDCVs)
 
    ! allocate buffer for transmittion
    allocate(loc_hbuffer(nrst_of_records),stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for transmitted history buffer!')
    end if

    ! copy data
    k = 1
    do i=1,nrst_of_hills
        loc_hbuffer(k) = last_buffer%heights(i)
        k = k + 1
        do j=1,NumOfMTDCVs
            loc_hbuffer(k) = last_buffer%values(i,j)
            k = k + 1
            loc_hbuffer(k) = last_buffer%widths(i,j)
            k = k + 1
        end do
    end do

    call cpmf_mtd_client_add_buffer_data(ret_st, &
                                        fclient_id,meta_step-nrst_of_hills+1,nrst_of_hills,loc_hbuffer)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to add history buffer to cpmd library!')
    end if

    deallocate(loc_hbuffer)

#endif

 return

end subroutine mtd_client_add_cpmf_buffer

!===============================================================================
! Subroutine:  mtd_client_get_cpmf_buffers
!===============================================================================

subroutine mtd_client_get_cpmf_buffers(last_buffer)

 use pmf_utils
 use mtd_dat
 use mtd_history

 implicit none
 type(MTDHistType),pointer   :: last_buffer
 ! -----------------------------------------------
#ifdef PMFLIB_NETWORK
 integer                            :: ret_st, buffer_id
 integer                            :: level,start,nrst_of_hills,nrst_of_records
 integer                            :: alloc_failed,i,j,k
 real(PMFDP),allocatable            :: loc_hbuffer(:)
#endif
 !-----------------------------------------------------------------------------

 if( .not. associated(last_buffer) ) return ! nothing to set

#ifdef PMFLIB_NETWORK
    buffer_id = 1
    do while(.true.)
        ! get buffer info - e.g. size
        call cpmf_mtd_client_get_buffer_info(ret_st,buffer_id,level,start,nrst_of_hills)

        if( ret_st .ne. 0 ) exit    ! no more buffers

        nrst_of_records = nrst_of_hills*(1 + 2*NumOfMTDCVs);

        ! allocate local buffer
        allocate(loc_hbuffer(nrst_of_records),stat = alloc_failed)

        if ( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for history buffer!')
        end if

        ! get data
        call cpmf_mtd_client_get_buffer_data(ret_st,buffer_id,loc_hbuffer)

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'Unable to get cpmd history buffer!')
        end if

        ! allocate new MTD buffer
        call mtd_history_allocate_buffer(last_buffer%next_history_buffer,nrst_of_hills)

        ! move to newly added buffer
        last_buffer => last_buffer%next_history_buffer

        ! copy data
        k = 1
        do i=1,nrst_of_hills
            last_buffer%heights(i) = loc_hbuffer(k)
            k = k + 1
            do j=1,NumOfMTDCVs
                last_buffer%values(i,j) = loc_hbuffer(k)
                k = k + 1
                last_buffer%widths(i,j) = loc_hbuffer(k)
                k = k + 1
            end do
        end do
        last_buffer%nrst_of_values = nrst_of_hills

        deallocate(loc_hbuffer)

        ! print data to restart file
        call mtd_history_print_buffer(last_buffer,level,start)

        buffer_id = buffer_id + 1
    end do

#endif

 return

end subroutine mtd_client_get_cpmf_buffers

!===============================================================================
! Subroutine:  mtd_client_unregister
!===============================================================================

subroutine mtd_client_unregister

 use pmf_utils
 use mtd_dat
 use pmf_timers

 implicit none
 type(MTDHistType),pointer   :: last_buffer
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

 write(PMF_OUT,*)
 write(PMF_OUT,10) trim(fserver)

 write(MTD_OUT,*)
 write(MTD_OUT,20)

#ifdef PMFLIB_NETWORK
    ! find last active buffer ----------------------------------------------------
    last_buffer => hill_history
    do while( associated(last_buffer) )
        if( .not. associated(last_buffer%next_history_buffer) ) exit
        last_buffer => last_buffer%next_history_buffer
    end do

    ! is the last buffer unprocessed?
    ! note: received buffers (from other walkers) have always nrst_of_values .eq. length_of_buffer
    if( last_buffer%nrst_of_values .ne. last_buffer%length_of_buffer ) then
        ! buffer is not completed - exchande data with server
        fserverupdate = .true.
        write(MTD_OUT,30)
        call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)
        call mtd_client_exchange_data()
        call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)
    end if

    call cpmf_mtd_client_unregister()
    write(MTD_OUT,40)
#endif

 call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)

 return

 10 format('>>> INFO: Removing registration from MTD server ',A)
 20 format('# [MTD-CLIENT] Removing registration from server')
 30 format('# [MTD-CLIENT]    Sending last unprocessed data')
 40 format('# [MTD-CLIENT]    Client was unregistered')

end subroutine mtd_client_unregister

!===============================================================================

end module mtd_client
