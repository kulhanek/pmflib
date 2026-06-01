!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module stm_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  stm_control_read_stm
!===============================================================================

subroutine stm_control_read_stm(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use stm_dat
    use stm_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! -----------------------------------------------
    logical                            :: fserver_enabled
    ! --------------------------------------------------------------------------

    call stm_init_dat

    fserver_enabled = .false.

    write(PMF_OUT,'(/,a)') '--- [stm] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'stm') ) then
        write(PMF_OUT,10)
        return
    end if

    ! process options from [stm] section
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('STM','fmode',fmode,0,1)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no stm - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

#ifdef PMFLIB_NETWORK
    if( prmfile_get_string_by_key(prm_fin,'fserverkey',fserverkey)) then
        write(PMF_OUT,110) trim(fserverkey)
        fserver_enabled = .true.
    end if

    if( .not. fserver_enabled ) then
        call pmf_utils_exit(PMF_OUT,1,'The string method requires specification of the access credentials to stm-server!')
    end if

#else
    fserver_enabled = .false.
    write(PMF_OUT,105)
    call pmf_utils_exit(PMF_OUT,1,'The string method is requested but a network support is not available in PMFLib!')
#endif

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('STM','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftensor',ftensor,'i12')
    call pmf_ctrl_check_integer_in_range('STM','ftensor',ftensor,0,1)

    call pmf_ctrl_read_integer(prm_fin,'fbeadid',fbeadid,'i12')

    ! read beadid from file if fbeadid == 0
    if( fbeadid .le. 0 ) then
        call pmf_ctrl_read_stritem(prm_fin,'fbeadidfile',fbeadidfile)

        call pmf_utils_open(STM_BEADID,fbeadidfile,'O')
        read(STM_BEADID,*,end=500,err=500) fbeadid
        close(STM_BEADID)
    end if

    call pmf_ctrl_check_integer('STM','fbeadid',fbeadid,0,CND_GT)

    call pmf_ctrl_read_integer(prm_fin,'fsteadylen',fsteadylen,'i12')
    call pmf_ctrl_check_integer_in_range('STM','fsteadylen',fsteadylen,0,100)

    bead_id = fbeadid
    stm_enabled = fmode .gt. 0

    return

500 call pmf_utils_exit(PMF_OUT,1,'Unable to read bead id!')

 10 format (' >> String method is disabled!')

#ifdef PMFLIB_NETWORK
110 format ('fserverkey                            = ',a)
#endif

#ifndef PMFLIB_NETWORK
105 format (' >> The string method is not compiled in (network support is required)!')
#endif

end subroutine stm_control_read_stm

!===============================================================================
! Subroutine:  stm_control_read_cvs
!===============================================================================

subroutine stm_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use stm_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! -----------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{STM}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fstmdef(1:1) .eq. '{' ) then
        grpname = fstmdef(2:len_trim(fstmdef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open goup with name from stmdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            call pmf_utils_exit(PMF_OUT,1,'[STM] Unable to open group {' // trim(grpname) // '}!')
        end if
        call stm_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fstmdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fstmdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[STM] Unable to load file: ' // trim(fstmdef) // '!')
        end if

        call stm_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)

end subroutine stm_control_read_cvs

!===============================================================================
! Subroutine:  stm_control_read_cvs_from_group
!===============================================================================

subroutine stm_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use stm_dat
    use stm_cvs

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME)    :: resname
    character(len=PRMFILE_MAX_LINE)        :: cvname
    integer                                :: i, j, alloc_failed
    logical                                :: eresult
    ! -----------------------------------------------------------------------------

    ! count number of sections in group
    NumOfSTMCVs = prmfile_count_group(prm_fin)

    if( NumOfSTMCVs .le. 0 ) then
    ! on CV in current or specified group
    fmode = 0
    stm_enabled = .false.
    write(PMF_OUT,100)
    return
    end if

    write(PMF_OUT,110) NumOfSTMCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(STMCVList(NumOfSTMCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[STM] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfSTMCVs
        call stm_cvs_reset_cv(STMCVList(i))
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
                 '[STM] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[STM] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        STMCVList(i)%cvindx = cv_common_find_cv(cvname)
        STMCVList(i)%cv => CVList(STMCVList(i)%cvindx)%cv

        ! read the rest of stm CV
        call stm_cvs_read_cv(prm_fin,STMCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfSTMCVs
        do j=i+1,NumOfSTMCVs
            if( STMCVList(i)%cvindx .eq. STMCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[MTD] Two different STM collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No CVs are defined. Adaptive Biasing Force method is switched off!')
110 format('Number of collective variables : ',I2)
130 format('== Reading collective variable #',I2.2)
140 format('   Collective variable name : ',a)

end subroutine stm_control_read_cvs_from_group

!===============================================================================

end module stm_control
