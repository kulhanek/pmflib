!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_control_read_abf
!===============================================================================

subroutine abf_control_read_abf(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use abf_dat
    use abf_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call abf_init_dat

    write(PMF_OUT,'(/,a)') '--- [abf] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'abf') ) then
        write(PMF_OUT,10)
        return
    end if

    ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('ABF','fmode',fmode,0,2)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no abf - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    call pmf_ctrl_read_integer(prm_fin,'fmask_mode',fmask_mode,'i12')
    call pmf_ctrl_check_integer_in_range('ABF','fmask_mode',fmask_mode,0,1)

    call pmf_ctrl_read_logical(prm_fin,'fapply_abf',fapply_abf)

    call pmf_ctrl_read_integer(prm_fin,'feimode',feimode,'i12')
    call pmf_ctrl_check_integer_in_range('ABF','feimode',feimode,0,3)

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('ABF','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'frestart',frestart)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'i12')
    call pmf_ctrl_check_integer('ABF','frstupdate',frstupdate,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftrjsample',ftrjsample,'i12')
    call pmf_ctrl_check_integer('ABF','ftrjsample',ftrjsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'fprint_icf',fprint_icf)
    call pmf_ctrl_read_logical(prm_fin,'fcache_icf',fcache_icf)
    call pmf_ctrl_read_logical(prm_fin,'frawicf',frawicf)

    select case(feimode)
        case(1)
            write(PMF_OUT,20)
            call pmf_ctrl_read_integer(prm_fin,'fhramp',fhramp,'i12')
            call pmf_ctrl_check_integer('ABF','fhramp',fhramp,0,CND_GT)
        case(2)
            write(PMF_OUT,30)
            call pmf_ctrl_read_integer(prm_fin,'fhramp_min',fhramp_min,'i12')
            call pmf_ctrl_check_integer('ABF','fhramp_min',fhramp_min,0,CND_GT)
            call pmf_ctrl_read_integer(prm_fin,'fhramp_max',fhramp_max,'i12')
            call pmf_ctrl_check_integer('ABF','fhramp_max',fhramp_max,0,CND_GT)
            if( fhramp_max .le. fhramp_min ) then
                call pmf_utils_exit(PMF_OUT,1,'[ABF] fhramp_max must be > fhramp_min!')
            end if
        case(3)
            write(PMF_OUT,40)
            call pmf_ctrl_read_integer(prm_fin,'fblock_size',fblock_size,'i12')
            call pmf_ctrl_check_integer('ABF','fblock_size',fblock_size,0,CND_GT)
        case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown extrapolation/interpolation mode!')
    end select

    ! network setup ----------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '--- [abf-walker] ---------------------------------------------------------------'

    ! try open section
    if( prmfile_open_section(prm_fin,'abf-walker') ) then

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
        call pmf_utils_exit(PMF_OUT,1,'[ABF] frestart cannot be ON if multiple-walkers approach is used!')
    end if

    abf_enabled = fmode .gt. 0

    return

 10 format (' >> Adaptive biasing force method is disabled!')
 20 format (/,'>> Linear ramp mode I (feimode == 1)')
 30 format (/,'>> Linear ramp mode II (feimode == 2)')
 40 format (/,'>> Block averages (feimode == 3)')

100 format (' >> Multiple-walkers ABF method is disabled!')
#ifndef PMFLIB_NETWORK
105 format (' >> Multiple-walkers ABF method is not compiled in!')
#endif
110 format ('fserverkey                             = ',a)

end subroutine abf_control_read_abf

!===============================================================================
! Subroutine:  abf_control_read_cvs
!===============================================================================

subroutine abf_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use abf_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)   :: grpname
    type(PRMFILE_TYPE)                  :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{ABF}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fabfdef(1:1) .eq. '{' ) then
        grpname = fabfdef(2:len_trim(fabfdef)-1)
         write(PMF_OUT,110) grpname
        ! open goup with name from abfdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            abf_enabled = .false.
            return
        end if
        call abf_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fabfdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fabfdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to load file: ' // trim(fabfdef) // '!')
        end if

        call abf_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)
130 format(' >> No {ABF} group was specified - disabling ABF method!')

end subroutine abf_control_read_cvs

!===============================================================================
! Subroutine:  abf_control_read_cvs_from_group
!===============================================================================

subroutine abf_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use abf_dat
    use abf_cvs

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)     :: cvname
    integer                             :: i, j, alloc_failed
    logical                             :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfABFCVs = prmfile_count_group(prm_fin)

    if( NumOfABFCVs .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        abf_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfABFCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(ABFCVList(NumOfABFCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfABFCVs
        call abf_cvs_reset_cv(ABFCVList(i))
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
                 '[ABF] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        ABFCVList(i)%cvindx = cv_common_find_cv(cvname)
        ABFCVList(i)%cv     => CVList(ABFCVList(i)%cvindx)%cv

        ! read the rest of abf CV
        call abf_cvs_read_cv(prm_fin,ABFCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfABFCVs
        do j=i+1,NumOfABFCVs
            if( ABFCVList(i)%cvindx .eq. ABFCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[ABF] Two different ABF collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No CVs are defined. Adaptive Biasing Force method is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine abf_control_read_cvs_from_group

!===============================================================================

end module abf_control
