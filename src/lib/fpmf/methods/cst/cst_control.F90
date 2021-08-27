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

module cst_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_control_read_con
!===============================================================================

subroutine cst_control_read_con(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cst_dat
    use cst_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call cst_init_dat

    write(PMF_OUT,'(/,a)') '--- [cst] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'cst') ) then
        write(PMF_OUT,10)
        return
    end if

    ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'I12')
    call pmf_ctrl_check_integer_in_range('CST','fmode',fmode,0,1)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no cst - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'I12')
    call pmf_ctrl_check_integer('CST','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'fplevel',fplevel,'I12')
    call pmf_ctrl_check_integer_in_range('CST','fplevel',fplevel,0,1)

    call pmf_ctrl_read_logical(prm_fin,'frestart',frestart)

    call pmf_ctrl_read_integer(prm_fin,'faccurst',faccurst,'I12')
    call pmf_ctrl_check_integer('CST','faccurst',faccurst,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'I12')
    call pmf_ctrl_check_integer('CST','frstupdate',frstupdate,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftrjsample',ftrjsample,'I12')
    call pmf_ctrl_check_integer('CST','ftrjsample',ftrjsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'fenthalpy',fenthalpy)
    call pmf_ctrl_read_logical(prm_fin,'fentropy',fentropy)

    call pmf_ctrl_read_real8_wunit(prm_fin,'fepotoffset',EnergyUnit,fepotoffset,'F10.1')
    call pmf_ctrl_read_real8_wunit(prm_fin,'fekinoffset',EnergyUnit,fekinoffset,'F10.1')

    call pmf_ctrl_read_integer(prm_fin,'flambdasolver',flambdasolver,'I12')
    call pmf_ctrl_check_integer_in_range('CST','flambdasolver',flambdasolver,0,2)

    if( flambdasolver .eq. 2 ) then
        call pmf_ctrl_read_real8(prm_fin,'frcond',frcond,'E12.4')
    end if

    call pmf_ctrl_read_real8(prm_fin,'flambdatol',fepotoffset,'E12.4')
    call pmf_ctrl_read_integer(prm_fin,'fmaxiter',fmaxiter,'I12')

    cst_enabled = fmode .gt. 0

    return

 10 format (' >> Constrained dynamics is disabled!')

end subroutine cst_control_read_con

!===============================================================================
! Subroutine:  cst_control_read_cvs
!===============================================================================

subroutine cst_control_read_cvs(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cst_dat
    use cst_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprm_fin
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{CST}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fcstdef(1:1) .eq. '{' ) then
        grpname = fcstdef(2:len_trim(fcstdef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open goup with name from metadef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            cst_enabled = .false.
            return
        end if
        call cst_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fcstdef)

        call prmfile_init(locprm_fin)

        if( .not. prmfile_read(locprm_fin,fcstdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to load file: ' // trim(fcstdef) // '!')
        end if

        call cst_control_read_cvs_from_group(locprm_fin)

        call prmfile_clear(locprm_fin)
    end if

    ! init constrained atom arrays for CST
    call cst_init_cst_atoms

    return

110 format('Constraints are read from group: ',A)
120 format('Constraints are read from file : ',A)
130 format(' >> No {CST} group was specified - disabling constrained dynamics!')

end subroutine cst_control_read_cvs

!===============================================================================
! Subroutine:  cst_control_read_cvs_from_group
!===============================================================================

subroutine cst_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use cst_dat
    use cst_constraints

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)     :: cvname
    integer                             :: i, j, alloc_failed
    logical                             :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfCONs = prmfile_count_group(prm_fin)

    if( NumOfCONs .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        cst_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfCONs

    ! allocate constraint list ----------------------
    allocate(CONList(NumOfCONs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to allocate memory for constraint data!')
    end if

    do i=1,NumOfCONs
        call cst_constraints_reset_con(CONList(i))
    end do

    ! enumerate sections ----------------------------
    eresult = prmfile_first_section(prm_fin)
    i = 1

    do while(eresult)
        eresult = prmfile_get_section_name(prm_fin,resname)
        write(PMF_OUT,*)
        write(PMF_OUT,130) i
        if( resname .ne. 'CV' ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[CST] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        CONList(i)%cvindx = cv_common_find_cv(cvname)
        CONList(i)%cv     => CVList(CONList(i)%cvindx)%cv

        call cst_constraints_read_con(prm_fin,CONList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfCONs
        do j=i+1,NumOfCONs
            if( CONList(i)%cvindx .eq. CONList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1,'[CST] Two different constraints share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No constraints are defined. Blue moon method is switched off!')
110 format('Number of constraints : ',I4)
130 format('== Reading constraint #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine cst_control_read_cvs_from_group

!===============================================================================

end module cst_control

