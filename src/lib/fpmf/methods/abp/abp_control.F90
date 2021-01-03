!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_control_read_abp
!===============================================================================

subroutine abp_control_read_abp(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use abp_dat
    use abp_init
    use pmf_unit
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call abp_init_dat

    write(PMF_OUT,'(/,a)') &
        '--- [abp] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'abp') ) then
        write(PMF_OUT,10)
        return
    end if

    ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('ABP','fmode',fmode,0,1)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no abp - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    call pmf_ctrl_read_real8_wunit(prm_fin,'fhbias',EnergyUnit,fhbias,'f12.2')
    call pmf_ctrl_check_real8_wunit('ABP','fhbias',EnergyUnit,fhbias,0.0d0,CND_GE,'f12.2')

    call pmf_ctrl_read_integer(prm_fin,'feimode',feimode,'i12')
    call pmf_ctrl_check_integer_in_range('ABP','feimode',feimode,1,1)

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('ABP','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'frestart',frestart)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'i12')
    call pmf_ctrl_check_integer('ABP','frstupdate',frstupdate,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftrjsample',ftrjsample,'i12')
    call pmf_ctrl_check_integer('ABP','ftrjsample',ftrjsample,0,CND_GE)

    select case(feimode)
        case(1)
            write(PMF_OUT,20)
            call pmf_ctrl_read_integer(prm_fin,'fhramp',fhramp,'i12')
            call pmf_ctrl_check_integer('ABP','fhramp',fhramp,0,CND_GT)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Unknown extrapolation/interpolation mode!')
    end select

   ! network setup ----------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '--- [abp-walker] ---------------------------------------------------------------'

    ! try open section
    if( prmfile_open_section(prm_fin,'abp-walker') ) then

#ifdef PMFLIB_NETWORK
    if( prmfile_get_string_by_key(prm_fin,'fserverkey',fserverkey)) then
        write(PMF_OUT,110) trim(fserverkey)
        fserver_enabled = .true.
    end if

    call pmf_ctrl_read_integer(prm_fin,'fserverupdate',fserverupdate,'i12')
    call pmf_ctrl_check_integer('ABF','fserverupdate',fserverupdate,0,CND_GT)

    call pmf_ctrl_read_logical(prm_fin,'fabortonmwaerr',fabortonmwaerr)
#else
    fserver_enabled = .false.
    write(PMF_OUT,105)
#endif
    ! network setup ----------------------------------------------------------------

    else
        write(PMF_OUT,100)
    endif

    ! restart is read from server
    if( fserver_enabled .and. frestart ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABP] frestart cannot be ON if multiple-walkers approach is used!')
    end if

    abp_enabled = fmode .gt. 0

    return

 10 format (' >> Adaptive biasing force method is disabled!')
 20 format (/,'>> Linear ramp mode I')

100 format (' >> Multiple-walkers ABP method is disabled!')
#ifndef PMFLIB_NETWORK
105 format (' >> Multiple-walkers ABP method is not compiled in!')
#endif
110 format ('fserverkey                             = ',a)

end subroutine abp_control_read_abp

!===============================================================================
! Subroutine:  abp_control_read_cvs
!===============================================================================

subroutine abp_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use abp_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{ABP}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fabpdef(1:1) .eq. '{' ) then
        grpname = fabpdef(2:len_trim(fabpdef)-1)
        write(PMF_OUT,110) trim(grpname)
        ! open goup with name from abpdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            abp_enabled = .false.
            return
        end if
        call abp_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fabpdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fabpdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Unable to load file: ' // trim(fabpdef) // '!')
        end if

        call abp_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)
130 format(' >> No {ABP} group was specified - disabling ABP method!')

end subroutine abp_control_read_cvs

!===============================================================================
! Subroutine:  abp_control_read_cvs_from_group
!===============================================================================

subroutine abp_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use abp_dat
    use abp_cvs

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME)    :: resname
    character(len=PRMFILE_MAX_LINE)        :: cvname
    integer                                :: i, alloc_failed
    logical                                :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfABPCVs = prmfile_count_group(prm_fin)

    if( NumOfABPCVs .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        abp_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfABPCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(ABPCVList(NumOfABPCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABP] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfABPCVs
        call abp_cvs_reset_cv(ABPCVList(i))
    end do

    ! enumerate sections ----------------------------------------------------------
    eresult = prmfile_first_section(prm_fin)
    i = 1
    do while(eresult)
        eresult = prmfile_get_section_name(prm_fin,resname)
        write(PMF_OUT,*)
        write(PMF_OUT,130) i
        if( resname .ne. 'CV' ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[ABP] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        ABPCVList(i)%cvindx = cv_common_find_cv(cvname)
        ABPCVList(i)%cv     => CVList(ABPCVList(i)%cvindx)%cv

        ! read the rest of abp CV
        call abp_cvs_read_cv(prm_fin,ABPCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    return

100 format('>>> INFO: No CVs are defined. Adaptive Biasing Force method is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine abp_control_read_cvs_from_group

!===============================================================================

end module abp_control
