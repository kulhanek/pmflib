!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtc_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtc_control_read_mtc
! load only [mtc] section
!===============================================================================

subroutine mtc_control_read_mtc(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use mtc_dat
    use mtc_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    call mtc_init_dat

    write(PMF_OUT,'(/,a)')   '--- [mtc] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'mtc') ) then
        write(PMF_OUT,5)
        return
    end if

    ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('MTC','fmode',fmode,0,1)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no mtc - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('MTC','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'frestart',frestart)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'I12')
    call pmf_ctrl_check_integer('MTC','frstupdate',frstupdate,0,CND_GE)

    mtc_enabled = fmode .gt. 0

    return

  5 format (' >> Metric Tensor Correction is disabled!')

end subroutine mtc_control_read_mtc

!===============================================================================
! Subroutine:  mtc_control_read_cvs
!===============================================================================

subroutine mtc_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use mtc_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_GROUP_NAME)  :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{MTC}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fmtcdef(1:1) .eq. '{' ) then
        grpname = fmtcdef(2:len_trim(fmtcdef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open goup with name from mtcdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            call pmf_utils_exit(PMF_OUT,1,'[MTC] Unable to open group {' // trim(grpname) // '}!')
        end if
        call mtc_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fmtcdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fmtcdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[MTC] Unable to load file: ' // trim(fmtcdef) // '!')
        end if

        call mtc_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)

end subroutine mtc_control_read_cvs

!===============================================================================
! Subroutine:  mtc_control_read_cvs_from_group
!===============================================================================

subroutine mtc_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use mtc_dat
    use mtc_cvs

    implicit none
    type(PRMFILE_TYPE),intent(inout)        :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)         :: cvname
    integer                                 :: i, alloc_failed
    logical                                 :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfMTCCVs = prmfile_count_group(prm_fin)

    if( NumOfMTCCVs .le. 0 ) then
        ! no CV in current or specified group
        fmode = 0
        mtc_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfMTCCVs

    ! allocate list of CVs indexes ------------------------------------------------
    allocate(MTCCVList(NumOfMTCCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTC] Unable to allocate memory for coordinate data!')
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
                 '[MTC] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[MTC] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)
        MTCCVList(i)%cvindx = cv_common_find_cv(cvname)
        MTCCVList(i)%cv => CVList(MTCCVList(i)%cvindx)%cv

        ! read the rest of MTC CV
        call mtc_cvs_read_cv(prm_fin,MTCCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    return

100 format('>>> INFO: No collective variables are defined. Monitoring is switched off!')
110 format('Number of collective variables : ',I2)
130 format('== Reading collective variable #',I2.2)
140 format('   Collective variable name : ',a)

end subroutine mtc_control_read_cvs_from_group

!===============================================================================

end module mtc_control
