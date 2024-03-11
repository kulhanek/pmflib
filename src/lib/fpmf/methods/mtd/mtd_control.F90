!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtd_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_control_read_mtd
!===============================================================================

subroutine mtd_control_read_mtd(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use mtd_dat
    use mtd_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call mtd_init_dat

    write(PMF_OUT,'(/,a)')   '--- [mtd] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'mtd') ) then
        write(PMF_OUT,10)
        return
    end if

   ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('MTD','fmode',fmode,0,1)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no abf - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    call pmf_ctrl_read_real8_wunit(prm_fin,'fheight',EnergyUnit,fheight,'F10.5')
    call pmf_ctrl_check_real8_wunit('MTD','fheight',EnergyUnit,fheight,0.0d0,CND_GE,'F10.5')

    call pmf_ctrl_read_integer(prm_fin,'fmetastep',fmetastep,'i12')
    call pmf_ctrl_check_integer('MTD','fmetastep',fmetastep,0,CND_GE)

    call pmf_ctrl_read_real8_wunit(prm_fin,'fmetatemp',TemperatureUnit,fmetatemp,'F10.1')
    call pmf_ctrl_check_real8_wunit('MTD','fmetatemp',TemperatureUnit,fmetatemp,0.0d0,CND_GE,'F10.1')

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('ABF','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'frestart',frestart)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'i12')
    call pmf_ctrl_check_integer('ABF','frstupdate',frstupdate,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftrjsample',ftrjsample,'i12')
    call pmf_ctrl_check_integer('ABF','ftrjsample',ftrjsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'fwritehills',fwritehills)

    call pmf_ctrl_read_logical(prm_fin,'fswitch2zero',fswitch2zero)

    ! network setup ----------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '--- [mtd-walker] -----------------------------------------------------------'

    ! try open section
    if( prmfile_open_section(prm_fin,'mtd-walker') ) then

#ifdef PMFLIB_NETWORK
    if( prmfile_get_string_by_key(prm_fin,'fserverkey',fserverkey)) then
        write(PMF_OUT,110) trim(fserverkey)
        fserver_enabled = .true.
    end if

    call pmf_ctrl_read_integer(prm_fin,'fserverupdate',fserverupdate,'i12')
    call pmf_ctrl_check_integer('MTD','fserverupdate',fserverupdate,0,CND_GT)

    call pmf_ctrl_read_logical(prm_fin,'fabortonmwaerr',fabortonmwaerr)

#else
    fserver_enabled = .false.
    fserver_key_enabled = .false.
    write(PMF_OUT,105)
#endif
    ! network setup ----------------------------------------------------------------

    else
        write(PMF_OUT,100)
    end if

    ! restart is read from server
    if( fserver_enabled .and. frestart ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] frestart cannot be on if multiple-walker approach is used!')
    end if

    mtd_enabled = fmode .gt. 0

    return

 10 format (' >> Metadynamics is disabled!')

100 format (' >> Multiple-walkers metadynamics is disabled!')
#ifndef PMFLIB_NETWORK
105 format (' >> Multiple-walkers metadynamics is not compiled in!')
#else
110 format ('fserverkey                             = ',a)
#endif

end subroutine mtd_control_read_mtd

!===============================================================================
! Subroutine:  mtd_control_read_cvs
!===============================================================================

subroutine mtd_control_read_cvs(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use mtd_dat

    implicit none
    type(PRMFILE_TYPE),intent(inout)        :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_GROUP_NAME)   :: grpname
    type(PRMFILE_TYPE)                      :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{MTD}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fmtddef(1:1) .eq. '{' ) then
        grpname = fmtddef(2:len_trim(fmtddef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open group with name from metadef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            mtd_enabled = .false.
            return
        end if
        call mtd_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fmtddef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fmtddef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to load file: ' // trim(fmtddef) // '!')
        end if

        call mtd_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Metadynamics collective variables are read from group: ',A)
120 format('Metadynamics collective variables are read from file : ',A)
130 format(' >> No {MTD} group was specified - disabling metadynamics!')

end subroutine mtd_control_read_cvs

!===============================================================================
! Subroutine:  mtd_control_read_cvs_from_group
!===============================================================================

subroutine mtd_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use mtd_dat
    use mtd_cvs_mod

    implicit none
    type(PRMFILE_TYPE),intent(inout)        :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)         :: cvname
    integer                                 :: i, j, alloc_failed
    logical                                 :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfMTDCVs = prmfile_count_group(prm_fin)

    if( NumOfMTDCVs .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        mtd_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfMTDCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(MTDCVList(NumOfMTDCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to allocate memory for coordinate data!')
    end if

    ! enumerate sections ----------------------------------------------------------
    eresult = prmfile_first_section(prm_fin)
    i = 1

    do while(eresult)
        eresult = prmfile_get_section_name(prm_fin,resname)
        write(PMF_OUT,*)
        write(PMF_OUT,130) i
        if( resname .ne. 'CV' ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[MTD] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[MTD] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        MTDCVList(i)%cvindx = cv_common_find_cv(cvname)
        MTDCVList(i)%cv => CVList(MTDCVList(i)%cvindx)%cv

        ! read the rest of mtd CV
        call mtd_cvs_read_cv(prm_fin,MTDCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfMTDCVs
        do j=i+1,NumOfMTDCVs
            if( MTDCVList(i)%cvindx .eq. MTDCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[MTD] Two different MTD collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No collective variables are defined. Metadynamics is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine mtd_control_read_cvs_from_group

!===============================================================================

end module mtd_control
