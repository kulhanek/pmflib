!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module gap_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  gap_control_read_gap
! load only [gap] section
!===============================================================================

subroutine gap_control_read_gap(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use gap_dat
    use gap_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    call gap_init_dat

    write(PMF_OUT,'(/,a)')   '--- [gap] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'gap') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [gap] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[GAP] fmode item is mandatory in this section')
    else
        write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0 or 1'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no gap - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fsample', fsample)) then
        write(PMF_OUT,50) fsample
    else
        write(PMF_OUT,55) fsample
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fplevel', fplevel)) then
        write(PMF_OUT,60) fplevel
    else
        write(PMF_OUT,65) fplevel
    end if

    gap_enabled = fmode .gt. 0

    return

  5 format (' >> GAP is disabled!')
 10 format ('fmode                                  = ',i12)
 15 format ('fmode                                  = ',i12,'                  (default)')
 50 format ('fsample                                = ',i12)
 55 format ('fsample                                = ',i12,'                  (default)')
 60 format ('fplevel                                = ',i12)
 65 format ('fplevel                                = ',i12,'                  (default)')

end subroutine gap_control_read_gap

!===============================================================================
! Subroutine:  gap_control_read_cvs
!===============================================================================

subroutine gap_control_read_cvs(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use gap_dat

    implicit none
    type(PRMFILE_TYPE),intent(inout)        :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_GROUP_NAME)   :: grpname
    type(PRMFILE_TYPE)                      :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{GAP}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fgapdef(1:1) .eq. '{' ) then
        grpname = fgapdef(2:len_trim(fgapdef)-1)
         write(PMF_OUT,110) grpname
        ! open group with name from gapdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            gap_enabled = .false.
            return
        end if
        call gap_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fgapdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fgapdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[GAP] Unable to load file: ' // trim(fgapdef) // '!')
        end if

        call gap_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('GAP collective variables are read from group: ',A)
120 format('GAP collective variables are read from file : ',A)
130 format(' >> No {GAP} group was specified - disabling GAP!')

end subroutine gap_control_read_cvs

!===============================================================================
! Subroutine:  gap_control_read_cvs_from_group
!===============================================================================

subroutine gap_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use gap_dat
    use gap_cvs

    implicit none
    type(PRMFILE_TYPE),intent(inout)        :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)         :: cvname
    integer                                 :: i, alloc_failed
    logical                                 :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfGAPCVs = prmfile_count_group(prm_fin)

    if( NumOfGAPCVs .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        gap_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfGAPCVs

    ! allocate gap cv list --------------------------------------------------------
    allocate(GAPCVList(NumOfGAPCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[GAP] Unable to allocate memory for coordinate data!')
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
                 '[GAP] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[GAP] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        GAPCVList(i)%cvindx = cv_common_find_cv(cvname)
        GAPCVList(i)%cv => CVList(GAPCVList(i)%cvindx)%cv

        ! read the rest of gap CV
        call gap_cvs_read_cv(prm_fin,GAPCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    return

100 format('>>> INFO: No collective variables are defined. GAP is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine gap_control_read_cvs_from_group

!===============================================================================

end module gap_control
