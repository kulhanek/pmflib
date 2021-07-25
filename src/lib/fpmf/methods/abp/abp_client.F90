!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_client

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interfaces
!===============================================================================

interface
    ! set number of coordinates
    subroutine cpmf_abp_client_set_header(ret_st,ncvs,nbins,version,driver, &
                    temp,temp_unit,temp_fconv, &
                    ene_unit,ene_fconv,widths)
        integer         :: ret_st
        integer         :: ncvs
        integer         :: nbins
        character(*)    :: version
        character(*)    :: driver
        real(8)         :: temp
        character(*)    :: temp_unit
        real(8)         :: temp_fconv
        character(*)    :: ene_unit
        real(8)         :: ene_fconv
        real(8)         :: widths(*)
    end subroutine cpmf_abp_client_set_header

    ! set coordinate data
    subroutine cpmf_abp_client_set_coord(ret_st,id,name,ctype,min_value,max_value,nbins,fconv,unit)
        integer         :: ret_st
        integer         :: id
        character(*)    :: ctype
        character(*)    :: name
        real(8)         :: min_value
        real(8)         :: max_value
        integer         :: nbins
        real(8)         :: fconv
        character(*)    :: unit
    end subroutine cpmf_abp_client_set_coord

    ! register client on server
    subroutine cpmf_abp_client_reg_by_key(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_abp_client_reg_by_key

    ! get initial data from server
    subroutine cpmf_abp_client_initial_data(ret_st,nisamples,idpop,ipop)
        integer         :: ret_st
        integer         :: nisamples(*)
        real(8)         :: idpop(*)
        real(8)         :: ipop(*)
    end subroutine cpmf_abp_client_initial_data

    ! exchange data with server
    subroutine cpmf_abp_client_exchange_data(ret_st,nisamples,idpop,ipop)
        integer         :: ret_st
        integer         :: nisamples(*)
        real(8)         :: idpop(*)
        real(8)         :: ipop(*)
    end subroutine cpmf_abp_client_exchange_data

    ! unregister client on server
    subroutine cpmf_abp_client_unregister()
    end subroutine cpmf_abp_client_unregister
end interface

contains

!===============================================================================
! Subroutine:  abp_register_client
!===============================================================================

subroutine abp_client_register

    use abp_dat
    use pmf_utils
    use pmf_timers
    use pmf_unit
    use pmf_ver

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: i
    integer        :: ret_st = 0
#endif
    ! --------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_ABP_MWA_TIMER)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Multiple-walkers ABP method','=')
    write(PMF_OUT,15) trim(fserverkey)

    write(PMF_OUT,*)
    write(PMF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! register coordinates
    call cpmf_abp_client_set_header(ret_st,abpaccu%tot_cvs,abpaccu%tot_nbins,PMFLIBVER,DriverName, &
                                    ftemp,trim(pmf_unit_label(TemperatureUnit)),pmf_unit_get_rvalue(TemperatureUnit,1.0d0), &
                                    trim(pmf_unit_label(EnergyUnit)),pmf_unit_get_rvalue(EnergyUnit,1.0d0), &
                                    abpaccu%widths)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1)
    end if

    do i=1,abpaccu%tot_cvs
        call cpmf_abp_client_set_coord(ret_st,      &
                        i,                          &
                        abpaccu%sizes(i)%cv%name,   &
                        abpaccu%sizes(i)%cv%ctype,  &
                        abpaccu%sizes(i)%min_value, &
                        abpaccu%sizes(i)%max_value, &
                        abpaccu%sizes(i)%nbins,     &
                        pmf_unit_get_rvalue(abpaccu%sizes(i)%cv%unit,1.0d0), &
                        pmf_unit_label(abpaccu%sizes(i)%cv%unit)             &
                        )

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1)
        end if
    end do

    ! register client
    call cpmf_abp_client_reg_by_key(fserverkey,fserver,client_id)

    write(ABP_OUT,*)

    if( client_id .le. 0 ) then
        fserver_enabled = .false.
        write(PMF_OUT,30)
        write(ABP_OUT,50) trim(fserver)
        call pmf_utils_exit(PMF_OUT,1)
    else
        write(PMF_OUT,40) client_id
        write(ABP_OUT,60) trim(fserver)
        write(ABP_OUT,70) client_id
        write(ABP_OUT,*)
    end if

#endif

    call pmf_timers_stop_timer(PMFLIB_ABP_MWA_TIMER)

    return

 15 format(' ABP Server Key file : ', A)
 20 format(' Registering client on server, please wait .... ')

 30 format(' Registration FAILED!')
 40 format(' Registration SUCCESSFULL! (Client ID: ',I6,')')

 50 format('# [ABP-CLIENT] Registration to server ',A,' failed!')
 60 format('# [ABP-CLIENT] Registration to server ',A,' successful.')
 70 format('# [ABP-CLIENT] Client ID: ',I6)

end subroutine abp_client_register

!===============================================================================
! Subroutine:  abp_client_get_initial_data
!===============================================================================

subroutine abp_client_get_initial_data

    use abp_dat
    use pmf_utils
    use pmf_timers
    use abp_accu

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: ret_st = 0
#endif
    integer        :: i
    ! -----------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_ABP_MWA_TIMER)

#ifdef PMFLIB_NETWORK
    call cpmf_abp_client_initial_data(ret_st,               &
                                        abpaccu%nsamples,   &
                                        abpaccu%dpop,       &
                                        abpaccu%pop         &
                                        )

    if( ret_st .ne. 0 ) then
        write(ABP_OUT,20)
        write(PMF_OUT,*)
        call pmf_utils_exit(PMF_OUT,1)
    end if

    write(ABP_OUT,10)
    write(PMF_OUT,*)
#endif

    ! update M
    do i=1,abpaccu%tot_nbins
        if( (abpaccu%pop(i)+1.0d0) .gt. abpaccu%M ) abpaccu%M = abpaccu%pop(i) + 1.0d0
    end do

    call pmf_timers_stop_timer(PMFLIB_ABP_MWA_TIMER)

    return

 10 format('# [ABP-CLIENT] Initial data from server were uploaded')
 20 format('# [ABP-CLIENT] Unable to get initial data from server - disabling method')

end subroutine abp_client_get_initial_data

!===============================================================================
! Subroutine:  abp_client_exchange_data
!===============================================================================

subroutine abp_client_exchange_data(force_exchange)

    use abp_dat
    use pmf_utils
    use pmf_timers
    use pmf_exit

    implicit none
    logical        :: force_exchange       ! do exchange even if mod(fstep,fserverupdate) .ne. 0
    ! -----------------------------------------------
#ifdef PMFLIB_NETWORK
    integer        :: ret_st = 0
#endif
    integer        :: i
    ! --------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    if( .not. force_exchange ) then
        ! do we need to exchange data?
        if( fserverupdate .le. 0 ) return
        if( mod(fstep,fserverupdate) .ne. 0 ) return
    end if

    call pmf_timers_start_timer(PMFLIB_ABP_MWA_TIMER)

    write(ABP_OUT,10) fstep

#ifdef PMFLIB_NETWORK
    call cpmf_abp_client_exchange_data(ret_st,                 &
                                       abpaccu%inc_nsamples,   &
                                       abpaccu%inc_dpop,       &
                                       abpaccu%inc_pop         &
                                       )

    if( ret_st .ne. 0 ) then
        failure_counter = failure_counter + 1
        write(ABP_OUT,20) failure_counter,fconrepeats
        if( failure_counter .gt. fconrepeats ) then
            if( fabortonmwaerr ) then
                call pmf_utils_exit(PMF_OUT,1,'MWA connection failures reach the treshold!')
            else
                call pmf_exit_mdloop(PMF_OUT,1,'MWA connection failures reach the treshold!')
            end if
        end if
        call pmf_timers_stop_timer(PMFLIB_ABP_MWA_TIMER)
        return
    end if

    ! move received data to main abpaccu
    abpaccu%nsamples(:)    = abpaccu%inc_nsamples(:)
    abpaccu%dpop(:,:)      = abpaccu%inc_dpop(:,:)
    abpaccu%pop(:)         = abpaccu%inc_pop(:)
#endif

    ! update M
    do i=1,abpaccu%tot_nbins
        if( (abpaccu%pop(i)+1.0d0) .gt. abpaccu%M ) abpaccu%M = abpaccu%pop(i) + 1.0d0
    end do

    ! and reset incremental data
    abpaccu%inc_nsamples(:)       = 0
    abpaccu%inc_dpop(:,:)         = 0.0d0
    abpaccu%inc_pop(:)            = 0.0d0

    failure_counter = 0

    call pmf_timers_stop_timer(PMFLIB_ABP_MWA_TIMER)

    return

 10 format('# [ABP-CLIENT] ',I9,' Exchanging data with server')
 20 format('# [ABP-CLIENT] Unable to exchange data with server - number of connection failures: ',I3,' (Treshold:',I3,')')

end subroutine abp_client_exchange_data

!===============================================================================
! Subroutine:  abp_client_unregister
!===============================================================================

subroutine abp_client_unregister

    use abp_dat
    use pmf_utils
    use pmf_timers

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: i
    logical        :: new_data
#endif
    ! --------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_ABP_MWA_TIMER)

    write(PMF_OUT,*)
    write(PMF_OUT,10) trim(fserver)

    write(ABP_OUT,*)
    write(ABP_OUT,20)

#ifdef PMFLIB_NETWORK
    ! test if there are new accumulated data
    new_data = .false.
    do i=1,abpaccu%tot_nbins
        if( abpaccu%inc_nsamples(i) .gt. 0 ) new_data = .true.
    end do
    if( new_data ) then
        write(ABP_OUT,30)
        call pmf_timers_stop_timer(PMFLIB_ABP_MWA_TIMER)
        call abp_client_exchange_data(.true.)
        call pmf_timers_start_timer(PMFLIB_ABP_MWA_TIMER)
    end if

    call cpmf_abp_client_unregister()
    write(ABP_OUT,40)
#endif

    call pmf_timers_stop_timer(PMFLIB_ABP_MWA_TIMER)

    return

 10 format('>>> INFO: Removing registration from ABP server ',A)
 20 format('# [ABP-CLIENT] Removing registration from server')
 30 format('# [ABP-CLIENT]    Sending last unprocessed data')
 40 format('# [ABP-CLIENT]    Client was unregistered')

end subroutine abp_client_unregister

!===============================================================================

end module abp_client
