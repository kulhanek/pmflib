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

module remd_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  remd_control_read_remd
! load only [remd] section
!===============================================================================

subroutine remd_control_read_remd(prm_fin)

 use prmfile
 use pmf_dat
 use pmf_utils
 use remd_dat
 use remd_init

 implicit none
 type(PRMFILE_TYPE),intent(inout)   :: prm_fin
 ! -----------------------------------------------------------------------------

 call remd_init_dat

 write(PMF_OUT,'(/,a)') '--- [remd] ---------------------------------------------------------------------'

 ! try open group
 if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
    write(PMF_OUT,100)
    return
 end if

 ! try open section
 if( .not. prmfile_open_section(prm_fin,'remd') ) then
    write(PMF_OUT,100)
    return
 end if

 ! process options from section
 if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
    call pmf_utils_exit(PMF_OUT,1,'[REMD] fmode item is mandatory in this section')
 else
    write(PMF_OUT,10) fmode
 end if

 if (fmode .ne. 0 .and. fmode .ne. 1) then
    write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0 or 1'
    call pmf_utils_exit(PMF_OUT,1)
 end if

 if( fmode .eq. 0 ) then
    write(PMF_OUT,100)
    ! no remd - rest of section is skipped
    call prmfile_set_sec_as_processed(prm_fin)
    return
 end if

! network setup ----------------------------------------------------------------
#ifdef PMFLIB_NETWORK
    if( prmfile_get_string_by_key(prm_fin,'fserverkey',fserverkey)) then
        write(PMF_OUT,110) trim(fserverkey)
        use_key = .true.
    end if

    if( .not. use_key ) then
        if( prmfile_get_string_by_key(prm_fin,'fserver', fserver) ) then
            write(PMF_OUT,120) fserver
        else
            call pmf_utils_exit(PMF_OUT,1,'fserver/fserverkey is required when [remd] is requested')
        end if
        if( prmfile_get_string_by_key(prm_fin,'fpassword', fpassword) ) then
            write(PMF_OUT,130) fpassword
        else
            call pmf_utils_exit(PMF_OUT,1,'fpassword/fserverkey is required when [remd] is requested')
        end if
    end if
#else
    use_key = .false.
    write(PMF_OUT,105)
    fmode = 0
    return
#endif

 if( prmfile_get_integer_by_key(prm_fin,'fsample',fsample) ) then
    write(PMF_OUT,20) fsample
 else
    write(PMF_OUT,25) fsample
 end if

 remd_enabled = fmode .gt. 0

 return

 10      format ('fmode                                  = ',i12)
 15      format ('fmode                                  = ',i12,'                  (default)')
 20      format ('fsample                                = ',i12)
 25      format ('fsample                                = ',i12,'                  (default)')

100      format (' >> Replica exchange molecular dynamics is disabled!')
105      format (' >> Network subsystem is not built in PMFLib - disabling REMD!')

110      format ('fserverkey                             = ',a)
120      format ('fserver                                = ',a)
130      format ('fpassword                              = ',a16)

end subroutine remd_control_read_remd

!===============================================================================

end module remd_control
