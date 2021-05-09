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

module tabf_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  tabf_control_read_abf
!===============================================================================

subroutine tabf_control_read_abf(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use tabf_dat
    use tabf_init
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call tabf_init_dat

    write(PMF_OUT,'(/,a)') '--- [tabf] ---------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,10)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'tabf') ) then
        write(PMF_OUT,10)
        return
    end if

    ! read configuration
    call pmf_ctrl_read_integer(prm_fin,'fmode',fmode,'i12')
    call pmf_ctrl_check_integer_in_range('ABF','fmode',fmode,0,3)

    if( fmode .eq. 0 ) then
        write(PMF_OUT,10)
        ! no abf - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    call pmf_ctrl_read_logical(prm_fin,'fapply_abf',fapply_abf)

    call pmf_ctrl_read_integer(prm_fin,'feimode',feimode,'i12')
    call pmf_ctrl_check_integer_in_range('ABF','feimode',feimode,1,1)

    call pmf_ctrl_read_integer(prm_fin,'fsample',fsample,'i12')
    call pmf_ctrl_check_integer('ABF','fsample',fsample,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'frstupdate',frstupdate,'i12')
    call pmf_ctrl_check_integer('ABF','frstupdate',frstupdate,0,CND_GE)

    call pmf_ctrl_read_integer(prm_fin,'ftrjsample',ftrjsample,'i12')
    call pmf_ctrl_check_integer('ABF','ftrjsample',ftrjsample,0,CND_GE)

    call pmf_ctrl_read_logical(prm_fin,'fprint_icf',fprint_icf)

    call pmf_ctrl_read_logical(prm_fin,'fenthalpy',fenthalpy)
    call pmf_ctrl_read_logical(prm_fin,'fentropy',fentropy)
    call pmf_ctrl_read_real8(prm_fin,'fepotoffset',fepotoffset,'F10.1')
    call pmf_ctrl_read_real8(prm_fin,'fekinoffset',fekinoffset,'F10.1')

    select case(feimode)
        case(1)
            write(PMF_OUT,20)
            call pmf_ctrl_read_integer(prm_fin,'fhramp_min',fhramp_min,'i12')
            call pmf_ctrl_check_integer('ABF','fhramp_min',fhramp_min,0,CND_GE)
            call pmf_ctrl_read_integer(prm_fin,'fhramp_max',fhramp_max,'i12')
            call pmf_ctrl_check_integer('ABF','fhramp_max',fhramp_max,0,CND_GE)
            if( fhramp_max .lt. fhramp_min ) then
                call pmf_utils_exit(PMF_OUT,1,'[TABF] fhramp_max must be >= fhramp_min!')
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Unknown extrapolation/interpolation mode!')
    end select

    tabf_enabled = fmode .gt. 0

    return

 10 format (' >> Adaptive biasing force method is disabled!')
 20 format (/,'>> Linear ramp mode I (feimode == 1)')

end subroutine tabf_control_read_abf

!===============================================================================
! Subroutine:  tabf_control_read_cvs
!===============================================================================

subroutine tabf_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use tabf_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)   :: grpname
    type(PRMFILE_TYPE)                  :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{TABF}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( ftabfdef(1:1) .eq. '{' ) then
        grpname = ftabfdef(2:len_trim(ftabfdef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open goup with name from abfdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            tabf_enabled = .false.
            return
        end if
        call tabf_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(ftabfdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,ftabfdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Unable to load file: ' // trim(ftabfdef) // '!')
        end if

        call tabf_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)
130 format(' >> No {TABF} group was specified - disabling ABF method!')

end subroutine tabf_control_read_cvs

!===============================================================================
! Subroutine:  tabf_control_read_cvs_from_group
!===============================================================================

subroutine tabf_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use tabf_dat
    use tabf_cvs

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
        tabf_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfABFCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(ABFCVList(NumOfABFCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[TABF] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfABFCVs
        call tabf_cvs_reset_cv(ABFCVList(i))
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
                 '[TABF] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[TABF] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        ABFCVList(i)%cvindx = cv_common_find_cv(cvname)
        ABFCVList(i)%cv     => CVList(ABFCVList(i)%cvindx)%cv

        ! read the rest of abf CV
        call tabf_cvs_read_cv(prm_fin,ABFCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfABFCVs
        do j=i+1,NumOfABFCVs
            if( ABFCVList(i)%cvindx .eq. ABFCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[TABF] Two different ABF collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No CVs are defined. Adaptive Biasing Force method is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine tabf_control_read_cvs_from_group

!===============================================================================

end module tabf_control
