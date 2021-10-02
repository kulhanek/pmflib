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

module usabf_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  usabf_control_read_abf
!===============================================================================

subroutine usabf_control_read_abf(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use usabf_dat
    use usabf_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call usabf_init_dat

    write(PMF_OUT,'(/,a)') '--- [us-abf] -------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'us-abf') ) then
        write(PMF_OUT,10)
        return
    end if

    ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('US-ABF','fmode',fmode,0,5)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no abf - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if


    call pmf_ctrl_read_logical(prm_fin,'fcontbias',fcontbias)
    call pmf_ctrl_read_logical(prm_fin,'falignbias',falignbias)
    call pmf_ctrl_read_logical(prm_fin,'ftdsbiased',ftdsbiased)

    call pmf_ctrl_read_logical(prm_fin,'frestart',frestart)

    call pmf_ctrl_read_integer(prm_fin,'faccurst',faccurst,'i12')
    call pmf_ctrl_check_integer('US-ABF','faccurst',faccurst,-1,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('US-ABF','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'i12')
    call pmf_ctrl_check_integer('US-ABF','frstupdate',frstupdate,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftrjsample',ftrjsample,'i12')
    call pmf_ctrl_check_integer('US-ABF','ftrjsample',ftrjsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'fenthalpy',fenthalpy)
    call pmf_ctrl_read_logical(prm_fin,'fentropy',fentropy)

    call pmf_ctrl_read_logical(prm_fin,'fsmoothetot',fsmoothetot)

    call pmf_ctrl_read_real8_wunit(prm_fin,'fepotaverage',EnergyUnit,fepotaverage,'F10.1')
    call pmf_ctrl_read_real8_wunit(prm_fin,'fekinaverage',EnergyUnit,fekinaverage,'F10.1')

    if( (fmode .eq. 4) .or. (fmode .eq. 5) ) then
        call pmf_ctrl_read_integer(prm_fin,'gpr_len',gpr_len,'i12')
        call pmf_ctrl_read_real8(prm_fin,'gpr_width',gpr_width,'F12.3')
        call pmf_ctrl_read_real8(prm_fin,'gpr_noise',gpr_noise,'F12.3')
    end if

    usabf_enabled = fmode .gt. 0

    return

 10 format (' >> Umbrella Sampling / Adaptive Biasing Force (US-ABF) method is disabled!')

end subroutine usabf_control_read_abf

!===============================================================================
! Subroutine:  usabf_control_read_cvs
!===============================================================================

subroutine usabf_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use usabf_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)   :: grpname
    type(PRMFILE_TYPE)                  :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{US-ABF}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fusabfdef(1:1) .eq. '{' ) then
        grpname = fusabfdef(2:len_trim(fusabfdef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open goup with name from abfdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            usabf_enabled = .false.
            return
        end if
        call usabf_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fusabfdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fusabfdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to load file: ' // trim(fusabfdef) // '!')
        end if

        call usabf_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)
130 format(' >> No {US-ABF} group was specified - disabling ABF method!')

end subroutine usabf_control_read_cvs

!===============================================================================
! Subroutine:  usabf_control_read_cvs_from_group
!===============================================================================

subroutine usabf_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use usabf_dat
    use usabf_cvs

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)     :: cvname
    integer                             :: i, j, alloc_failed
    logical                             :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfUSABFCVs = prmfile_count_group(prm_fin)

    if( NumOfUSABFCVs .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        usabf_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfUSABFCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(USABFCVList(NumOfUSABFCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfUSABFCVs
        call usabf_cvs_reset_cv(USABFCVList(i))
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
                 '[US-ABF] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[US-ABF] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        USABFCVList(i)%cvindx = cv_common_find_cv(cvname)
        USABFCVList(i)%cv     => CVList(USABFCVList(i)%cvindx)%cv

        ! read the rest of abf CV
        call usabf_cvs_read_cv(prm_fin,USABFCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfUSABFCVs
        do j=i+1,NumOfUSABFCVs
            if( USABFCVList(i)%cvindx .eq. USABFCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[US-ABF] Two different ABF collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No CVs are defined. US-ABF method is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine usabf_control_read_cvs_from_group

!===============================================================================

end module usabf_control
