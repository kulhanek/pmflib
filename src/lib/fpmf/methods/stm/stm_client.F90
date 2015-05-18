!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module stm_client

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interfaces
!===============================================================================

interface
    ! set number of coordinates
    subroutine cpmf_stm_client_set_header(ret_st,nitems)
        integer         :: ret_st
        integer         :: nitems
    end subroutine cpmf_stm_client_set_header

    ! set coordinate data
    subroutine cpmf_stm_client_set_coord(ret_st,id,ctype,name)
        integer         :: ret_st
        integer         :: id
        character(*)    :: ctype
        character(*)    :: name
    end subroutine cpmf_stm_client_set_coord

    ! register client on server
    subroutine cpmf_stm_client_reg_by_name(str1,str2,cid,bid)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: cid
        integer         :: bid
    end subroutine cpmf_stm_client_reg_by_name

    ! register client on server
    subroutine cpmf_stm_client_reg_by_key(str1,str2,cid,bid)
        character(*)    :: str1
        character(*)    :: str2
        integer         :: cid
        integer         :: bid
    end subroutine cpmf_stm_client_reg_by_key

    ! exchange data with server
    subroutine cpmf_stm_client_exchange_data(ret_st,mode,isteps,ipos,rpmf,rfz)
        integer         :: ret_st
        integer         :: mode
        integer         :: isteps
        real(8)         :: ipos(*)
        real(8)         :: rpmf(*)
        real(8)         :: rfz(*)
    end subroutine cpmf_stm_client_exchange_data

    ! unregister client on server
    subroutine cpmf_stm_client_unregister()
    end subroutine cpmf_stm_client_unregister
end interface

contains

!===============================================================================
! Subroutine:  stm_register_client
!===============================================================================

subroutine stm_client_register

    use stm_dat
    use pmf_utils
    use pmf_timers
    use pmf_unit

    implicit none
#ifdef PMFLIB_NETWORK
    integer        :: i
    integer        :: ret_st = 0
#endif
    ! --------------------------------------------------------------------------

    call pmf_timers_start_timer(PMFLIB_STM_NET_TIMER)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'String method','=')

    if( use_key ) then
        write(PMF_OUT,15) trim(fserverkey)
    else
        write(PMF_OUT,10) trim(fserver)
    end if

    write(PMF_OUT,*)
    write(PMF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! register coordinates
    call cpmf_stm_client_set_header(ret_st,NumOfSTMCVs)

    if( ret_st .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1)
    end if

    do i=1,NumOfSTMCVs
        call cpmf_stm_client_set_coord(ret_st, i, &
                        STMCVList(i)%cv%ctype, &
                        STMCVList(i)%cv%name)

        if( ret_st .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1)
        end if
    end do

    ! register client
    if( use_key ) then
        call cpmf_stm_client_reg_by_key(fserverkey,fserver,client_id,bead_id)
    else
        call cpmf_stm_client_reg_by_name(fserver,fpassword,client_id,bead_id)
    end if

    if( client_id .le. 0 ) then
        write(PMF_OUT,30)
        write(STM_OUT,45) bead_id
        write(STM_OUT,50) trim(fserver)
        call pmf_utils_exit(PMF_OUT,1)
    else
        write(PMF_OUT,40) client_id
        write(STM_OUT,45) bead_id
        write(STM_OUT,60) trim(fserver)
        write(STM_OUT,70) client_id
    end if

    ! init start values - this is important for clients that join the server several times
    ! in the next round they directly start by Equilibration
    do i=1,NumOfSTMCVs
        STMCVList(i)%target_value = CVContext%CVsValues(STMCVList(i)%cvindx)
    end do

#endif

    call pmf_timers_stop_timer(PMFLIB_STM_NET_TIMER)

    return

 10 format(' STM Server          : ', A)
 15 format(' STM Server Key file : ', A)
 20 format(' Registering client on server, please wait .... ')

 30 format(' Registration FAILED!')
 40 format(' Registration SUCCESSFULL! (Client ID: ',I6,')')

 45 format('# [STM-CLIENT] Bead ID  : ',I6)
 50 format('# [STM-CLIENT] Registration to server ',A,' failed!')
 60 format('# [STM-CLIENT] Registration to server ',A,' successful.')
 70 format('# [STM-CLIENT] Client ID: ',I6)

end subroutine stm_client_register

!===============================================================================
! Subroutine:  stm_client_exchange_data
!===============================================================================

subroutine stm_client_exchange_data()

    use stm_dat
    use pmf_utils
    use pmf_timers
    use pmf_cvs
    use pmf_exit

    implicit none
    ! -----------------------------------------------
#ifdef PMFLIB_NETWORK
    integer         :: ret_st = 0
#endif
    integer         :: i
    ! -----------------------------------------------------------------------------

    if( (curstep .ne. stmsteps) .and. (fstep .ne. 1) ) return

    call pmf_timers_start_timer(PMFLIB_STM_NET_TIMER)

    write(STM_OUT,10) fstep,stmsteps

    ! scale down PMF and MTZ
    if( stmsteps .gt. 0 ) then
        MTZ = MTZ/stmsteps
        PMF = PMF/stmsteps
    end if

#ifdef PMFLIB_NETWORK
    call cpmf_stm_client_exchange_data(ret_st, stmmode, stmsteps, &
                                       beadpos, pmf, MTZ)

    if( ret_st .ne. 0 ) then
        call pmf_exit_mdloop(PMF_OUT,1,'Unable to exchange data with the server')
    end if
#endif

    ! print info about new mode
    select case(stmmode)
        case(BMO_INITIALIZATION)
            write(STM_OUT,100)
            do i=1,NumOfSTMCVs
                STMCVList(i)%startvalue = CVContext%CVsValues(STMCVList(i)%cvindx)
                STMCVList(i)%stopvalue  = beadpos(i)
            end do
            write(STM_OUT,400,advance='NO') '# Start  '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%startvalue)
            end do
            write(STM_OUT,*)
            write(STM_OUT,400,advance='NO') '# Stop   '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%stopvalue)
            end do
            write(STM_OUT,*)
            write(STM_OUT,400,advance='NO') '# Diff   '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%stopvalue-STMCVList(i)%startvalue)
            end do
            write(STM_OUT,*)
            write(STM_OUT,350) stmsteps
        case(BMO_ACCUMULATION)
            write(STM_OUT,200)
            do i=1,NumOfSTMCVs
                STMCVList(i)%target_value = beadpos(i)
            end do
            write(STM_OUT,400,advance='NO') '# Target '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%target_value)
            end do
            write(STM_OUT,*)
            write(STM_OUT,350) stmsteps
        case(BMO_EQUILIBRATION)
            write(STM_OUT,300)
            do i=1,NumOfSTMCVs
                STMCVList(i)%startvalue = STMCVList(i)%target_value
                STMCVList(i)%stopvalue  = beadpos(i)
            end do
            write(STM_OUT,400,advance='NO') '# Start  '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%startvalue)
            end do
            write(STM_OUT,*)
            write(STM_OUT,400,advance='NO') '# Stop   '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%stopvalue)
            end do
            write(STM_OUT,*)
            write(STM_OUT,400,advance='NO') '# Diff   '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%stopvalue-STMCVList(i)%startvalue)
            end do
            write(STM_OUT,*)
            write(STM_OUT,350) stmsteps
        case(BMO_PRODUCTION)
            write(STM_OUT,310)
            do i=1,NumOfSTMCVs
                STMCVList(i)%target_value = beadpos(i)
            end do
            write(STM_OUT,400,advance='NO') '# Target '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%target_value)
            end do
            write(STM_OUT,*)
            write(STM_OUT,350) stmsteps
        case(BMO_TERMINATE)
            write(STM_OUT,320)
            do i=1,NumOfSTMCVs
                STMCVList(i)%target_value = beadpos(i)
            end do
            write(STM_OUT,400,advance='NO') '# Target '
            do i=1,NumOfSTMCVs
                write(STM_OUT,410,advance='NO') STMCVList(i)%cv%get_rvalue(STMCVList(i)%target_value)
            end do
            write(STM_OUT,*)
            call pmf_exit_mdloop(PMF_OUT,1,'STM server requested client termination!') ! force MD loop termination
        case default
            call pmf_utils_exit(PMF_OUT,1,'Illegal STM mode after data exchange!')
    end select

    ! reset accumulators
    curstep = 0
    pmf(:) = 0.0d0
    MTZ(:,:) = 0.0d0

    call pmf_timers_stop_timer(PMFLIB_STM_NET_TIMER)

    return

 10 format('# [STM-CLIENT] Exchanging data with the server in step: ',I9,1X,I9)

100 format('# [STM-CLIENT] Initialization mode (I)')
200 format('# [STM-CLIENT] PMF accumulation (A)')
300 format('# [STM-CLIENT] Equilibration for the next run (E)')
310 format('# [STM-CLIENT] Final production run (P)')
320 format('# [STM-CLIENT] Client termination requested (T)')
350 format('# [STM-CLIENT] Number of steps : ',I9)

400 format(A9)
410 format(1X,F15.8)

end subroutine stm_client_exchange_data

!===============================================================================
! Subroutine:  stm_client_unregister
!===============================================================================

subroutine stm_client_unregister

    use stm_dat
    use pmf_utils
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_timers_start_timer(PMFLIB_STM_NET_TIMER)

    write(PMF_OUT,*)
    write(PMF_OUT,10) trim(fserver)

    write(STM_OUT,20)

#ifdef PMFLIB_NETWORK
    ! any intermediate data are discarded
    ! so directly unregister client
    call cpmf_stm_client_unregister()
    write(STM_OUT,40)
#endif

    call pmf_timers_stop_timer(PMFLIB_STM_NET_TIMER)

    return

10 format('>>> INFO: Removing registration from STM server ',A)
20 format('# [STM-CLIENT] Removing registration from server')
40 format('# [STM-CLIENT]    Client was unregistered')

end subroutine stm_client_unregister

!===============================================================================

end module stm_client
