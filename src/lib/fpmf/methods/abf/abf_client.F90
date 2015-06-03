!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module abf_client

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interfaces
!===============================================================================

interface
    ! set number of coordinates
    subroutine cpmf_abf_client_set_header(ret_st,nitems,totnbins)
        integer         :: ret_st
        integer         :: nitems
        integer         :: totnbins
    end subroutine cpmf_abf_client_set_header

    ! set coordinate data
    subroutine cpmf_abf_client_set_coord(ret_st,id,name,ctype,min_value,max_value,nbins)
        integer         :: ret_st
        integer         :: id
        character(*)    :: ctype
        character(*)    :: name
        real(8)         :: min_value
        real(8)         :: max_value
        integer         :: nbins
    end subroutine cpmf_abf_client_set_coord

    ! register client on server
    subroutine cpmf_abf_client_reg_by_name(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_abf_client_reg_by_name

    ! register client on server
    subroutine cpmf_abf_client_reg_by_key(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_abf_client_reg_by_key

    ! get initial data from server
    subroutine cpmf_abf_client_initial_data(ret_st,nisamples,iabfforce,iabfforce2)
        integer         :: ret_st
        integer         :: nisamples(*)
        real(8)         :: iabfforce(*)
        real(8)         :: iabfforce2(*)
    end subroutine cpmf_abf_client_initial_data

    ! exchange data with server
    subroutine cpmf_abf_client_exchange_data(ret_st,nisamples,iabfforce,iabfforce2)
        integer         :: ret_st
        integer         :: nisamples(*)
        real(8)         :: iabfforce(*)
        real(8)         :: iabfforce2(*)
    end subroutine cpmf_abf_client_exchange_data

    ! unregister client on server
    subroutine cpmf_abf_client_unregister()
    end subroutine cpmf_abf_client_unregister
end interface

contains

!===============================================================================
! Subroutine:  abf_register_client
!===============================================================================

subroutine abf_client_register

 use abf_dat
 use pmf_utils
 use pmf_timers
 use pmf_unit

 implicit none
#ifdef PMFLIB_NETWORK
 integer        :: i
 integer        :: ret_st = 0
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'Multiple-walkers ABF method','=')

 if( use_key ) then
    write(PMF_OUT,15) trim(fserverkey)
 else
    write(PMF_OUT,10) trim(fserver)
 end if
 write(PMF_OUT,*)
 write(PMF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! register coordinates
    call cpmf_abf_client_set_header(ret_st,NumOfABFCVs,accumulator%tot_nbins)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1)
    end if

    do i=1,NumOfABFCVs
        call cpmf_abf_client_set_coord(ret_st, &
                        i, &
                        ABFCVList(i)%cv%name, &
                        ABFCVList(i)%cv%ctype, &
                        ABFCVList(i)%min_value, &
                        ABFCVList(i)%max_value, &
                        ABFCVList(i)%nbins)

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1)
        end if
    end do

    ! register client
    if( use_key ) then
        call cpmf_abf_client_reg_by_key(fserverkey,fserver,client_id)
    else
        call cpmf_abf_client_reg_by_name(fserver,fpassword,client_id) 
    end if

    write(ABF_OUT,*)

    if( client_id .le. 0 ) then
        fserver_enabled = .false.
        write(PMF_OUT,30)
        write(ABF_OUT,50) trim(fserver)
        call pmf_utils_exit(PMF_OUT,1)
    else
        write(PMF_OUT,40) client_id
        write(ABF_OUT,60) trim(fserver)
        write(ABF_OUT,70) client_id
        write(ABF_OUT,*)
    end if

#endif

 call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)

 return

 10 format(' ABF Server          : ', A)
 15 format(' ABF Server Key file : ', A)
 20 format(' Registering client on server, please wait .... ')

 30 format(' Registration FAILED!')
 40 format(' Registration SUCCESSFULL! (Client ID: ',I6,')')

 50 format('# [ABF-CLIENT] Registration to server ',A,' failed!')
 60 format('# [ABF-CLIENT] Registration to server ',A,' successful.')
 70 format('# [ABF-CLIENT] Client ID: ',I6)

end subroutine abf_client_register

!===============================================================================
! Subroutine:  abf_client_get_initial_data
!===============================================================================

subroutine abf_client_get_initial_data

 use abf_dat
 use pmf_utils
 use pmf_timers
 use abf_accumulator

 implicit none
#ifdef PMFLIB_NETWORK
 integer        :: ret_st = 0
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

#ifdef PMFLIB_NETWORK
    call cpmf_abf_client_initial_data(ret_st, &
                                        accumulator%nsamples, &
                                        accumulator%abfforce, &
                                        accumulator%abfforce2)

    if( ret_st .ne. 0 ) then
        write(ABF_OUT,20)
        write(PMF_OUT,*)
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmask_mode .eq. 1 ) then
        write(ABF_OUT,40)
        call abf_accumulator_apply_mask()
    end if

    write(ABF_OUT,10)
    write(PMF_OUT,*)
#endif

 call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)

 return

 10 format('# [ABF-CLIENT] Initial data from server were uploaded')
 20 format('# [ABF-CLIENT] Unable to get initial data from server - disabling method')
 40 format('# MASK: applying mask to current ABF accumulator')

end subroutine abf_client_get_initial_data

!===============================================================================
! Subroutine:  abf_client_exchange_data
!===============================================================================

subroutine abf_client_exchange_data(force_exchange)

 use abf_dat
 use pmf_utils
 use pmf_timers
 use pmf_exit

 implicit none
 logical        :: force_exchange       ! do exchange even if mod(fstep,fserverupdate) .ne. 0
 ! -----------------------------------------------
#ifdef PMFLIB_NETWORK
 integer        :: ret_st = 0
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 if( .not. force_exchange ) then
    ! do we need to exchange data?
    if( fserverupdate .le. 0 ) return
    if( mod(fstep,fserverupdate) .ne. 0 ) return
 end if

 call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

 write(ABF_OUT,10) fstep

#ifdef PMFLIB_NETWORK
    call cpmf_abf_client_exchange_data(ret_st, &
                                       accumulator%nisamples, &
                                       accumulator%iabfforce, &
                                       accumulator%iabfforce2)

    if( ret_st .ne. 0 ) then
        failure_counter = failure_counter + 1
        write(ABF_OUT,20) failure_counter,fconrepeats
        if( failure_counter .gt. fconrepeats ) then
            call pmf_exit_mdloop(PMF_OUT,1,'MWA connection failures reach the treshold!')
        end if
        call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)
        return
    end if
#endif

 ! move received data to main accumulator
 accumulator%nsamples(:) = accumulator%nisamples(:)
 accumulator%abfforce(:,:) = accumulator%iabfforce(:,:)
 accumulator%abfforce2(:,:) = accumulator%iabfforce2(:,:)

 ! and reset incremental data
 accumulator%nisamples(:) = 0
 accumulator%iabfforce(:,:) = 0.0d0
 accumulator%iabfforce2(:,:) = 0.0d0

 failure_counter = 0

 call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)

 return

 10 format('# [ABF-CLIENT] ',I9,' Exchanging data with server')
 20 format('# [ABF-CLIENT] Unable to exchange data with server - number of connection failures: ',I3,' (Treshold:',I3,')')

end subroutine abf_client_exchange_data

!===============================================================================
! Subroutine:  abf_client_unregister
!===============================================================================

subroutine abf_client_unregister

 use abf_dat
 use pmf_utils
 use pmf_timers

 implicit none
#ifdef PMFLIB_NETWORK
 integer        :: i
 logical        :: new_data
#endif
 ! -----------------------------------------------------------------------------

 if( .not. fserver_enabled ) return

 call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

 write(PMF_OUT,*)
 write(PMF_OUT,10) trim(fserver)

 write(ABF_OUT,*)
 write(ABF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! test if there are new accumulated data
    new_data = .false.
    do i=1,accumulator%tot_nbins
        if( accumulator%nisamples(i) .gt. 0 ) new_data = .true.
    end do
    if( new_data ) then
        write(ABF_OUT,30)
        call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)
        call abf_client_exchange_data(.true.)
        call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)
    end if

    call cpmf_abf_client_unregister()
    write(ABF_OUT,40)
#endif

 call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)

 return

 10 format('>>> INFO: Removing registration from ABF server ',A)
 20 format('# [ABF-CLIENT] Removing registration from server')
 30 format('# [ABF-CLIENT]    Sending last unprocessed data')
 40 format('# [ABF-CLIENT]    Client was unregistered')

end subroutine abf_client_unregister

!===============================================================================

end module abf_client
