!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    subroutine cpmf_abf_client_set_header(ret_st,ncvs,nbins,version,driver, &
                    temp,temp_unit,temp_fconv, &
                    ene_unit,ene_fconv,enthalpy_enabled,mwa_mode)
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
        integer         :: enthalpy_enabled
        integer         :: mwa_mode
    end subroutine cpmf_abf_client_set_header

    ! set coordinate data
    subroutine cpmf_abf_client_set_coord(ret_st,id,name,ctype,min_value,max_value,nbins,fconv,unit)
        integer         :: ret_st
        integer         :: id
        character(*)    :: ctype
        character(*)    :: name
        real(8)         :: min_value
        real(8)         :: max_value
        integer         :: nbins
        real(8)         :: fconv
        character(*)    :: unit
    end subroutine cpmf_abf_client_set_coord

    ! register client on server
    subroutine cpmf_abf_client_reg_by_key(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_abf_client_reg_by_key

    ! get initial data from server
    subroutine cpmf_abf_client_initial_data(ret_st,inc_nsamples,inc_micf,inc_m2icf, &
                                            inc_mepot,inc_m2epot)
        integer         :: ret_st
        real(8)         :: inc_nsamples(*)
        real(8)         :: inc_micf(*)
        real(8)         :: inc_m2icf(*)
        real(8)         :: inc_mepot(*)
        real(8)         :: inc_m2epot(*)
    end subroutine cpmf_abf_client_initial_data

    ! exchange data with server
    subroutine cpmf_abf_client_exchange_data(ret_st,inc_nsamples,inc_micf,inc_m2icf, &
                                            inc_mepot,inc_m2epot)
        integer         :: ret_st
        real(8)         :: inc_nsamples(*)
        real(8)         :: inc_micf(*)
        real(8)         :: inc_m2icf(*)
        real(8)         :: inc_mepot(*)
        real(8)         :: inc_m2epot(*)
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
    use pmf_ver

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: i
    integer        :: ret_st = 0
#endif
    integer        :: enthalpy_enabled
    ! --------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Multiple-walkers ABF method','=')
    write(PMF_OUT,15) trim(fserverkey)

    write(PMF_OUT,*)
    write(PMF_OUT,20)

    enthalpy_enabled = 0
    if( fenthalpy ) enthalpy_enabled = 1

#ifdef PMFLIB_NETWORK

        ! register coordinates
        call cpmf_abf_client_set_header(ret_st,abfaccu%tot_cvs,abfaccu%tot_nbins,PMFLIBVER,DriverName, &
                                        ftemp,trim(pmf_unit_label(TemperatureUnit)),pmf_unit_get_rvalue(TemperatureUnit,1.0d0), &
                                        trim(pmf_unit_label(EnergyUnit)),pmf_unit_get_rvalue(EnergyUnit,1.0d0), &
                                        enthalpy_enabled,fmwamode)

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1)
        end if

        do i=1,abfaccu%tot_cvs
            call cpmf_abf_client_set_coord(ret_st,      &
                            i,                          &
                            abfaccu%sizes(i)%cv%name,   &
                            abfaccu%sizes(i)%cv%ctype,  &
                            abfaccu%sizes(i)%min_value, &
                            abfaccu%sizes(i)%max_value, &
                            abfaccu%sizes(i)%nbins,     &
                            pmf_unit_get_rvalue(abfaccu%sizes(i)%cv%unit,1.0d0), &
                            pmf_unit_label(abfaccu%sizes(i)%cv%unit)             &
                            )

            if( ret_st .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1)
            end if
        end do

        ! register client
        call cpmf_abf_client_reg_by_key(fserverkey,fserver,client_id)

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
    use abf_accu

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: ret_st = 0
#endif
    ! -----------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

#ifdef PMFLIB_NETWORK
    call cpmf_abf_client_initial_data(ret_st,               &
                                    abfaccu%inc_nsamples,   &
                                    abfaccu%inc_micf,       &
                                    abfaccu%inc_m2icf,      &
                                    abfaccu%inc_mepot,      &
                                    abfaccu%inc_m2epot)
    if( ret_st .ne. 0 ) then
        write(ABF_OUT,20)
        write(PMF_OUT,*)
        call pmf_utils_exit(PMF_OUT,1)
    end if

    select case(fmwamode)
        case(0)
            ! move received data to main abfaccu
            abfaccu%nsamples(:)         = abfaccu%inc_nsamples(:)
            abfaccu%micf(:,:)           = abfaccu%inc_micf(:,:)
            abfaccu%m2icf(:,:)          = abfaccu%inc_m2icf(:,:)

            abfaccu%bnsamples(:)        = abfaccu%inc_nsamples(:)
            abfaccu%bmicf(:,:)          = abfaccu%inc_micf(:,:)

            if( fenthalpy ) then
                abfaccu%mepot(:)        = abfaccu%inc_mepot(:)
                abfaccu%m2epot(:)       = abfaccu%inc_m2epot(:)
            end if

        case(1)
            abfaccu%bnsamples(:)        = abfaccu%inc_nsamples(:)
            abfaccu%bmicf(:,:)          = abfaccu%inc_micf(:,:)
        case default
            call pmf_utils_exit(PMF_OUT,1,'fmwamode not implemented in abf_client_get_initial_data')
    end select

    abfaccu%inc_nsamples(:) = 0
    abfaccu%inc_micf(:,:)   = 0.0d0
    abfaccu%inc_m2icf(:,:)  = 0.0d0
    abfaccu%inc_mepot(:)    = 0.0d0
    abfaccu%inc_m2epot(:)   = 0.0d0

    write(ABF_OUT,10)
    write(PMF_OUT,*)
#endif

    call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)

    return

 10 format('# [ABF-CLIENT] Initial data from server were uploaded')
 20 format('# [ABF-CLIENT] Unable to get initial data from server - disabling method')

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
    ! --------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    if( .not. force_exchange ) then
        ! do we need to exchange data?
        if( fserverupdate .le. 0 ) return
        if( mod(fstep,fserverupdate) .ne. 0 ) return
    end if

    call pmf_timers_start_timer(PMFLIB_ABF_MWA_TIMER)

    write(ABF_OUT,10) fstep

#ifdef PMFLIB_NETWORK
    call cpmf_abf_client_exchange_data(ret_st,                  &
                                        abfaccu%inc_nsamples,   &
                                        abfaccu%inc_micf,       &
                                        abfaccu%inc_m2icf,      &
                                        abfaccu%inc_mepot,      &
                                        abfaccu%inc_m2epot )

    if( ret_st .ne. 0 ) then
        failure_counter = failure_counter + 1
        write(ABF_OUT,20) failure_counter,fconrepeats
        if( failure_counter .gt. fconrepeats ) then
            if( fabortonmwaerr ) then
                call pmf_utils_exit(PMF_OUT,1,'MWA connection failures reach the treshold!')
            else
                call pmf_exit_mdloop(PMF_OUT,1,'MWA connection failures reach the treshold!')
            end if
        end if
        call pmf_timers_stop_timer(PMFLIB_ABF_MWA_TIMER)
        return
    end if

! move received data to main abfaccu
    select case(fmwamode)
        case(0)
            abfaccu%nsamples(:)         = abfaccu%inc_nsamples(:)
            abfaccu%micf(:,:)           = abfaccu%inc_micf(:,:)
            abfaccu%m2icf(:,:)          = abfaccu%inc_m2icf(:,:)

            abfaccu%bnsamples(:)        = abfaccu%inc_nsamples(:)
            abfaccu%bmicf(:,:)          = abfaccu%inc_micf(:,:)

            if( fenthalpy ) then
                abfaccu%mepot(:)        = abfaccu%inc_mepot(:)
                abfaccu%m2epot(:)       = abfaccu%inc_m2epot(:)
            end if
        case(1)
            abfaccu%bnsamples(:)        = abfaccu%inc_nsamples(:)
            abfaccu%bmicf(:,:)          = abfaccu%inc_micf(:,:)
        case default
            call pmf_utils_exit(PMF_OUT,1,'fmwamode not implemented in abf_client_get_initial_data')
    end select

    abfaccu%inc_nsamples(:) = 0
    abfaccu%inc_micf(:,:)   = 0.0d0
    abfaccu%inc_m2icf(:,:)  = 0.0d0
    abfaccu%inc_mepot(:)    = 0.0d0
    abfaccu%inc_m2epot(:)   = 0.0d0

#endif

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
    do i=1,abfaccu%tot_nbins
        if( abfaccu%inc_nsamples(i) .gt. 0 ) new_data = .true.
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
