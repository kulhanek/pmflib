!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf2_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf2_control_read_abf
!===============================================================================

subroutine abf2_control_read_abf(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use abf2_dat
    use abf2_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call abf2_init_dat

    write(PMF_OUT,'(/,a)') '--- [abf2] ---------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'abf2') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [abf] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF2] fmode item is mandatory in this section')
    else
     write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1 .and. fmode .ne. 2 ) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0, 1, 2'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no abf - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fmask_mode', fmask_mode)) then
        write(PMF_OUT,16) fmask_mode
    else
        write(PMF_OUT,19) fmask_mode
    end if

    if (fmask_mode .ne. 0 .and. fmask_mode .ne. 1) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmask_mode (', fmask_mode, ') must be 0 or 1'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_logical_by_key(prm_fin,'fapply_abf', fapply_abf)) then
        write(PMF_OUT,21) prmfile_onoff(fapply_abf)
    else
        write(PMF_OUT,26) prmfile_onoff(fapply_abf)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fsample', fsample)) then
        write(PMF_OUT,50) fsample
    else
        write(PMF_OUT,55) fsample
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fblock_size', fblock_size)) then
        write(PMF_OUT,90) fblock_size
    else
        write(PMF_OUT,94) fblock_size
    end if

    if( fblock_size .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'fblock_size has to be > 0')
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fintrpl', fintrpl)) then
        write(PMF_OUT,96) fintrpl
    else
        write(PMF_OUT,98) fintrpl
    end if

    if (fintrpl .ne. 0 .and. fintrpl .ne. 1) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fintrpl (', fintrpl, ') must be 0 or 1'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_logical_by_key(prm_fin,'frestart', frestart)) then
        write(PMF_OUT,70) prmfile_onoff(frestart)
    else
        write(PMF_OUT,75) prmfile_onoff(frestart)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'frstupdate', frstupdate)) then
        write(PMF_OUT,76) frstupdate
    else
        write(PMF_OUT,77) frstupdate
    end if

    if( frstupdate .lt. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'frstupdate has to be >=0')
    end if

    if(prmfile_get_integer_by_key(prm_fin,'ftrjsample', ftrjsample)) then
        write(PMF_OUT,80) ftrjsample
    else
        write(PMF_OUT,85) ftrjsample
    end if

    if( ftrjsample .lt. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'ftrjsample has to be >=0')
    end if

    ! network setup ----------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '--- [abf-walker] ---------------------------------------------------------------'

    ! try open section
    if( prmfile_open_section(prm_fin,'abf-walker') ) then

#ifdef PMFLIB_NETWORK
    if( prmfile_get_string_by_key(prm_fin,'fserverkey',fserverkey)) then
        write(PMF_OUT,110) trim(fserverkey)
        fserver_enabled = .true.
        use_key = .true.
    end if

    if( .not. fserver_enabled ) then
        if( prmfile_get_string_by_key(prm_fin,'fserver', fserver) ) then
            write(PMF_OUT,120) fserver
            fserver_enabled = .true.
            use_key = .false.
        else
            call pmf_utils_exit(PMF_OUT,1,'fserver is required when [abf-walker] is specified')
        end if
        if( prmfile_get_string_by_key(prm_fin,'fpassword', fpassword) ) then
            write(PMF_OUT,130) fpassword
        else
            call pmf_utils_exit(PMF_OUT,1,'fpassword is required when [abf-walker] is specified')
        end if
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fserverupdate', fserverupdate)) then
        write(PMF_OUT,140) fserverupdate
    else
        write(PMF_OUT,145) fserverupdate
    end if

    if(prmfile_get_logical_by_key(prm_fin,'fabortonmwaerr', fabortonmwaerr)) then
        write(PMF_OUT,150) prmfile_onoff(fabortonmwaerr)
    else
        write(PMF_OUT,155) prmfile_onoff(fabortonmwaerr)
    end if

#else
    fserver_enabled = .false.
    use_key = .false.
    write(PMF_OUT,105)
#endif
    ! network setup ----------------------------------------------------------------

    else
        write(PMF_OUT,100)
    endif

    ! restart is read from server
    if( fserver_enabled .and. frestart ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF2] frestart cannot be ON if multiple-walkers approach is used!')
    end if

    abf2_enabled = fmode .gt. 0

    return

  5 format (' >> Adaptive biasing force 2 method is disabled!')
 10 format ('fmode                                  = ',i12)
 16 format ('fmask_mode                             = ',i12)
 19 format ('fmask_mode                             = ',i12,'                  (default)')
 21 format ('fapply_abf                             = ',a12)
 26 format ('fapply_abf                             = ',a12,'                  (default)')
 50 format ('fsample                                = ',i12)
 55 format ('fsample                                = ',i12,'                  (default)')
 90 format ('fblock_size                            = ',i12)
 94 format ('fblock_size                            = ',i12,'                  (default)')
 96 format ('fintrpl                                = ',i12)
 98 format ('fintrpl                                = ',i12,'                  (default)')
 70 format ('frestart                               = ',a12)
 75 format ('frestart                               = ',a12,'                  (default)')
 76 format ('frstupdate                             = ',i12)
 77 format ('frstupdate                             = ',i12,'                  (default)')
 80 format ('ftrjsample                             = ',i12)
 85 format ('ftrjsample                             = ',i12,'                  (default)')

100 format (' >> Multiple-walkers ABF2 method is disabled!')
#ifndef PMFLIB_NETWORK
105 format (' >> Multiple-walkers ABF2 method is not compiled in!')
#endif
110 format ('fserverkey                             = ',a)
120 format ('fserver                                = ',a)
130 format ('fpassword                              = ',a16)
140 format ('fserverupdate                          = ',i12)
145 format ('fserverupdate                          = ',i12,'                  (default)')
150 format ('fabortonmwaerr                         = ',a12)
155 format ('fabortonmwaerr                         = ',a12,'                  (default)')

end subroutine abf2_control_read_abf

!===============================================================================
! Subroutine:  abf2_control_read_cvs
!===============================================================================

subroutine abf2_control_read_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use abf2_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)   :: grpname
    type(PRMFILE_TYPE)                  :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{ABF2}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fabfdef(1:1) .eq. '{' ) then
        grpname = fabfdef(2:len_trim(fabfdef)-1)
         write(PMF_OUT,110) grpname
        ! open goup with name from abfdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            abf2_enabled = .false.
            return
        end if
        call abf2_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fabfdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fabfdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF2] Unable to load file: ' // trim(fabfdef) // '!')
        end if

        call abf2_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)
130 format(' >> No {ABF} group was specified - disabling ABF2 method!')

end subroutine abf2_control_read_cvs

!===============================================================================
! Subroutine:  abf2_control_read_cvs_from_group
!===============================================================================

subroutine abf2_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use abf2_dat
    use abf2_cvs

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
        abf2_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfABFCVs

    ! allocate constraint list ----------------------------------------------------
    allocate(ABFCVList(NumOfABFCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF2] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfABFCVs
        call abf2_cvs_reset_cv(ABFCVList(i))
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
                 '[ABF2] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF2] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        ABFCVList(i)%cvindx = cv_common_find_cv(cvname)
        ABFCVList(i)%cv     => CVList(ABFCVList(i)%cvindx)%cv

        ! read the rest of abf CV
        call abf2_cvs_read_cv(prm_fin,ABFCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfABFCVs
        do j=i+1,NumOfABFCVs
            if( ABFCVList(i)%cvindx .eq. ABFCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[ABF2] Two different ABF collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No CVs are defined. Adaptive Biasing Force method is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine abf2_control_read_cvs_from_group

!===============================================================================

end module abf2_control
