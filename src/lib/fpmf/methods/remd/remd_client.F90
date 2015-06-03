!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module remd_client

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================
! CPMF/FPMF interfaces
!===============================================================================

interface
    ! register client on server
    subroutine cpmf_remd_client_reg_by_name(natoms,str1,str2,id)
        integer         :: natoms
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_remd_client_reg_by_name

    ! register client on server
    subroutine cpmf_remd_client_reg_by_key(natoms,str1,str2,id)
        integer         :: natoms
        character(*)    :: str1
        character(*)    :: str2
        integer         :: id
    end subroutine cpmf_remd_client_reg_by_key

    ! exchange data with server
    subroutine cpmf_remd_client_initial_data(ret_st,mode,period,bath_id,temp)
        integer         :: ret_st
        integer         :: mode
        integer         :: period
        integer         :: bath_id
        real(8)         :: temp
    end subroutine cpmf_remd_client_initial_data

    ! exchange data with server
    subroutine cpmf_remd_client_exchange_data(ret_st,epot,bath_id,ctemp,otemp)
        integer         :: ret_st
        real(8)         :: epot
        integer         :: bath_id
        real(8)         :: ctemp
        real(8)         :: otemp
    end subroutine cpmf_remd_client_exchange_data

    ! unregister client on server
    subroutine cpmf_remd_client_unregister()
    end subroutine cpmf_remd_client_unregister
end interface

contains

!===============================================================================
! Subroutine:  remd_register_client
!===============================================================================

subroutine remd_client_register

 use pmf_dat
 use remd_dat
 use pmf_utils

 implicit none
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'REMD walker','=')

 if( use_key ) then
    write(PMF_OUT,15) trim(fserverkey)
 else
    write(PMF_OUT,10) trim(fserver)
 end if
 write(PMF_OUT,*)
 write(PMF_OUT,20)

#ifdef PMFLIB_NETWORK
    ! register client
    if( use_key ) then
        call cpmf_remd_client_reg_by_key(fnatoms,fserverkey,fserver,ReplicaId)
    else
        call cpmf_remd_client_reg_by_name(fnatoms,fserver,fpassword,ReplicaId)
    end if

    write(REMD_OUT,*)

    if( ReplicaId .lt. 0 ) then
        write(PMF_OUT,30)
        write(REMD_OUT,50) trim(fserver)
        call pmf_utils_exit(PMF_OUT,1,'Client/server communication error!')
    end if

    write(PMF_OUT,40) ReplicaId
    write(PMF_OUT,*)

    ! write info to REMDOutput
    write(REMD_OUT,60) trim(fserver)
    write(REMD_OUT,70) ReplicaId

#endif

 return

 10 format(' REMD Server          : ', A)
 15 format(' REMD Server Key file : ', A)
 20 format(' Registering client on server, please wait .... ')

 30 format(' Registration FAILED!')
 40 format(' Registration SUCCESSFULL! (Replica ID: ',I6,')')

 50 format('# Registration to server ',A,' failed!')
 60 format('# Registration to server ',A,' successful.')
 70 format('# Replica ID: ',I6)

end subroutine remd_client_register

!===============================================================================
! Subroutine:  remd_client_initial_data
!===============================================================================

subroutine remd_client_get_initial_data()

 use pmf_dat
 use remd_dat
 use pmf_utils

 implicit none
#ifdef PMFLIB_NETWORK
 integer        :: ret_st = 0
 integer        :: remd_mode
#endif
 ! -----------------------------------------------------------------------------

#ifdef PMFLIB_NETWORK
    write(REMD_OUT,*)

    remd_mode = 0
    call cpmf_remd_client_initial_data(ret_st,remd_mode,REMDPeriod,BathId,CurBathTemp)
    OldBathTemp = CurBathTemp
    REMDFullSwap = remd_mode .eq. 1
    if( ret_st .ne. 0 ) then
        write(REMD_OUT,20)
        call pmf_utils_exit(PMF_OUT,1,'Client/server communication error!')
    end if

    if( REMDFullSwap ) then
        write(REMD_OUT,35)
    else
        write(REMD_OUT,30)
    end if

    write(REMD_OUT,40) REMDPeriod
    write(REMD_OUT,45) BathId
    write(REMD_OUT,50) CurBathTemp
#endif

 return

 20 format('# ERROR: Unable to get initial data from server!')

 30 format('# Initial data: Exchange mode   = temperature swap')
 35 format('# Initial data: Exchange mode   = phase point swap')
 40 format('# Initial data: Exchange period = ',I6)
 45 format('# Initial data: Bath ID         = ',I6)
 50 format('# Initial data: Temperature     = ',F10.3)

end subroutine remd_client_get_initial_data

!===============================================================================
! Subroutine:  remd_client_exchange_data
!===============================================================================

subroutine remd_client_exchange_data()

 use pmf_dat
 use remd_dat
 use pmf_utils
 use pmf_pbc

 implicit none
#ifdef PMFLIB_NETWORK
 integer        :: ret_st = 0
 real(PMFDP)    :: box_abc(3),box_angs(3)
#endif
 ! -----------------------------------------------------------------------------

 write(REMD_OUT,10) fstep

#ifdef PMFLIB_NETWORK
    if( REMDFullSwap ) then
        ! get box dimmensions
        call pmf_pbc_get_box(box_abc(1),box_abc(2),box_abc(3),box_angs(1),box_angs(2),box_angs(3))
        ! set snapshot
        call cpmf_remd_client_set_data(ret_st,Crd,Vel,box_abc,box_angs)
    end if

    OldBathTemp = CurBathTemp
    call cpmf_remd_client_exchange_data(ret_st,PotEne,BathId,CurBathTemp,OldBathTemp)

    if( ret_st .ne. 0 ) then
        write(REMD_OUT,20)
        call pmf_utils_exit(PMF_OUT,1,'Client/server communication error!')
    end if

    ! get snapshot data if necessary
    if( REMDFullSwap ) then
        ! get snapshot
        call cpmf_remd_client_get_data(ret_st,Crd,Vel,box_abc,box_angs)
        ! set box dimmensions
        call pmf_pbc_set_box(box_abc(1),box_abc(2),box_abc(3),box_angs(1),box_angs(2),box_angs(3))
    end if

    ! rescale velocities according to exchange mode
    if( REMDFullSwap ) then
        write(REMD_OUT,35) fstep,OldBathTemp
        ! data are always exchanged
        ! all coordinates and velocities were exchanged by server
    else
        write(REMD_OUT,30) fstep,CurBathTemp,OldBathTemp
        ! rescale velocities
        Vel = sqrt(CurBathTemp/OldBathTemp)*Vel
    end if
#endif

 return

 10 format('# ',I10,' Exchanging data with server')
 20 format('# ERROR: Unable to exchange data with server!')
 30 format('# ',I10,' New temperature  = ',F10.3,' (Previous temperature = ',F10.3,')')
 35 format('# ',I10,' Previous temperature = ',F10.3)

end subroutine remd_client_exchange_data

!===============================================================================
! Subroutine:  remd_client_unregister
!===============================================================================

subroutine remd_client_unregister

 use remd_dat
 use pmf_utils

 implicit none
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,*)
 write(PMF_OUT,10) trim(fserver)

 write(REMD_OUT,*)
 write(REMD_OUT,20)

#ifdef PMFLIB_NETWORK
    call cpmf_remd_client_unregister()
    write(REMD_OUT,40)
#endif

 return

 10 format('>>> INFO: Removing registration from REMD server ',A)
 20 format('# Removing registration from server')
 40 format('#    Client was unregistered')

end subroutine remd_client_unregister

!===============================================================================

end module remd_client
