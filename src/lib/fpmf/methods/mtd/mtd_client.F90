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

module mtd_client

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interfaces
!===============================================================================

interface
    ! set number of coordinates
    subroutine cpmf_mtd_client_set_header(ret_st,ncvs,nbins,version,driver, &
                    temp,temp_unit,temp_fconv, &
                    ene_unit,ene_fconv,widths,wt_temp)
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
        real(8)         :: wt_temp
    end subroutine cpmf_mtd_client_set_header

    ! set coordinate data
    subroutine cpmf_mtd_client_set_coord(ret_st,id,name,ctype,min_value,max_value,nbins,fconv,unit)
        integer         :: ret_st
        integer         :: id
        character(*)    :: ctype
        character(*)    :: name
        real(8)         :: min_value
        real(8)         :: max_value
        integer         :: nbins
        real(8)         :: fconv
        character(*)    :: unit
    end subroutine cpmf_mtd_client_set_coord

    ! register client on server
    subroutine cpmf_mtd_client_reg_by_key(str1,str2,id)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_mtd_client_reg_by_key

    ! get initial data from server
    subroutine cpmf_mtd_client_initial_data(ret_st,insamples,imtdpot,imtdforces)
        integer         :: ret_st
        integer         :: insamples(*)
        real(8)         :: imtdforces(*)
        real(8)         :: imtdpot(*)
    end subroutine cpmf_mtd_client_initial_data

    ! exchange data with server
    subroutine cpmf_mtd_client_exchange_data(ret_st,insamples,imtdpot,imtdforces)
        integer         :: ret_st
        integer         :: insamples(*)
        real(8)         :: imtdforces(*)
        real(8)         :: imtdpot(*)
    end subroutine cpmf_mtd_client_exchange_data

    ! unregister client on server
    subroutine cpmf_mtd_client_unregister()
    end subroutine cpmf_mtd_client_unregister
end interface

contains

!===============================================================================
! Subroutine:  mtd_register_client
!===============================================================================

subroutine mtd_client_register

    use mtd_dat
    use pmf_utils
    use pmf_timers
    use pmf_unit
    use pmf_ver
    use pmf_dat

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: i
    integer        :: ret_st = 0
#endif
    ! -----------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Multiple-walkers MTD method','=')
    write(PMF_OUT,15) trim(fserverkey)

    write(PMF_OUT,*)
    write(PMF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! register coordinates
    call cpmf_mtd_client_set_header(ret_st,mtdaccu%tot_cvs,mtdaccu%tot_nbins,PMFLIBVER,DriverName, &
                                    ftemp,trim(pmf_unit_label(TemperatureUnit)),pmf_unit_get_rvalue(TemperatureUnit,1.0d0), &
                                    trim(pmf_unit_label(EnergyUnit)),pmf_unit_get_rvalue(EnergyUnit,1.0d0), &
                                    mtdaccu%widths,fmetatemp)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1)
    end if

    do i=1,NumOfMTDCVs
        call cpmf_mtd_client_set_coord(ret_st,      &
                        i,                          &
                        mtdaccu%sizes(i)%cv%name,   &
                        mtdaccu%sizes(i)%cv%ctype,  &
                        mtdaccu%sizes(i)%min_value, &
                        mtdaccu%sizes(i)%max_value, &
                        mtdaccu%sizes(i)%nbins,     &
                        pmf_unit_get_rvalue(mtdaccu%sizes(i)%cv%unit,1.0d0), &
                        pmf_unit_label(mtdaccu%sizes(i)%cv%unit)             &
                        )

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1)
        end if
    end do

    ! register client
    call cpmf_mtd_client_reg_by_key(fserverkey,fserver,client_id)

    write(MTD_OUT,*)

    if( client_id .le. 0 ) then
        fserver_enabled = .false.
        write(PMF_OUT,30)
        write(MTD_OUT,50) trim(fserver)
        call pmf_utils_exit(PMF_OUT,1)
    else
        write(PMF_OUT,40) client_id
        write(MTD_OUT,60) trim(fserver)
        write(MTD_OUT,70) client_id
        write(MTD_OUT,*)
    end if

#endif

    call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)

    return

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

    use mtd_dat
    use pmf_utils
    use pmf_timers
    use mtd_accu

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: ret_st = 0
#endif
    ! -----------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

#ifdef PMFLIB_NETWORK
    call cpmf_mtd_client_initial_data(ret_st,               &
                                        mtdaccu%nsamples,   &
                                        mtdaccu%mtdforce,   &
                                        mtdaccu%mtdpot      &
                                        )

    if( ret_st .ne. 0 ) then
        write(MTD_OUT,20)
        write(PMF_OUT,*)
        call pmf_utils_exit(PMF_OUT,1)
    end if

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

subroutine mtd_client_exchange_data(force_exchange)

    use mtd_dat
    use pmf_dat
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

    call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

    write(MTD_OUT,10) fstep

#ifdef PMFLIB_NETWORK
    call cpmf_mtd_client_exchange_data(ret_st,                  &
                                        mtdaccu%inc_nsamples,   &
                                        mtdaccu%inc_mtdforce,   &
                                        mtdaccu%inc_mtdpot      &
                                       )

    if( ret_st .ne. 0 ) then
        failure_counter = failure_counter + 1
        write(MTD_OUT,20) failure_counter,fconrepeats
        if( failure_counter .gt. fconrepeats ) then
            if( fabortonmwaerr ) then
                call pmf_utils_exit(PMF_OUT,1,'MWA connection failures reach the treshold!')
            else
                call pmf_exit_mdloop(PMF_OUT,1,'MWA connection failures reach the treshold!')
            end if
        end if
        call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)
        return
    end if

    ! move received data to main accumulator
    mtdaccu%nsamples(:)     = mtdaccu%inc_nsamples(:)
    mtdaccu%mtdpot(:)       = mtdaccu%inc_mtdpot(:)
    mtdaccu%mtdforce(:,:)   = mtdaccu%inc_mtdforce(:,:)
#endif

    ! and reset incremental data
    mtdaccu%inc_nsamples(:)     = 0
    mtdaccu%inc_mtdpot(:)       = 0.0d0
    mtdaccu%inc_mtdforce(:,:)   = 0.0d0

    failure_counter = 0

    call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)

    return

10 format('# [MTD-CLIENT] ',I9,' Exchanging data with server')
20 format('# [MTD-CLIENT] Unable to exchange data with server - number of connection failures: ',I3,' (Treshold:',I3,')')

end subroutine mtd_client_exchange_data

!===============================================================================
! Subroutine:  mtd_client_unregister
!===============================================================================

subroutine mtd_client_unregister

    use mtd_dat
    use pmf_utils
    use pmf_timers

    implicit none
    ! -----------------------------------------------------------------------------

    if( .not. fserver_enabled ) return

    call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

    write(PMF_OUT,*)
    write(PMF_OUT,10) trim(fserver)

    write(MTD_OUT,*)
    write(MTD_OUT,20)

#ifdef PMFLIB_NETWORK
    ! test if there are new accumulated data
    write(MTD_OUT,30)
    call pmf_timers_stop_timer(PMFLIB_MTD_MWA_TIMER)
    call mtd_client_exchange_data(.true.)
    call pmf_timers_start_timer(PMFLIB_MTD_MWA_TIMER)

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
