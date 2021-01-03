!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'stm') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [stm] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[STM] fmode item is mandatory in this section')
    else
        write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0 or 1'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no stm - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

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
            call pmf_utils_exit(PMF_OUT,1,'fserver is required when STM method is requested')
        end if
        if( prmfile_get_string_by_key(prm_fin,'fpassword', fpassword) ) then
            write(PMF_OUT,130) fpassword
        else
            call pmf_utils_exit(PMF_OUT,1,'fpassword is required when STM method is requested')
        end if
    end if

    if( .not. fserver_enabled ) then
        call pmf_utils_exit(PMF_OUT,1,'String method requires specification of access to stm-server!')
    end if

#else
    fserver_enabled = .false.
    use_key = .false.
    write(PMF_OUT,105)
    call pmf_utils_exit(PMF_OUT,1,'String method is requested but network support is not available in PMFLib!')
#endif

    if(prmfile_get_integer_by_key(prm_fin,'fsample', fsample)) then
        write(PMF_OUT,50) fsample
    else
        write(PMF_OUT,55) fsample
    end if

    if(prmfile_get_integer_by_key(prm_fin,'ftensor', ftensor)) then
        write(PMF_OUT,60) ftensor
    else
        write(PMF_OUT,65) ftensor
    end if

    if( ftensor .ne. 0 .and. ftensor .ne. 1 .and. ftensor .ne. 2 ) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'ftensor (', ftensor, ') must be 0, 1, or 2'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_string_by_key(prm_fin,'fbeadidfile', fbeadidfile)) then
        write(PMF_OUT,140) fbeadidfile
    else
        write(PMF_OUT,145) fbeadidfile
    end if

    stm_enabled = (fmode .gt. 0) .and. (fserver_enabled .eqv. .true.)

    if( .not. stm_enabled ) then
        write(PMF_OUT,5)
    end if

    ! read bead id ------------------------------
    write(PMF_OUT,*)
    write(PMF_OUT,150) trim(fbeadidfile)

    call pmf_utils_open(STM_BEADID,fbeadidfile,'O')
    read(STM_BEADID,*,end=500,err=500) bead_id
    close(STM_BEADID)
    write(PMF_OUT,155) bead_id

    if( bead_id .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Bead ID must be greater than zero!')
    end if

    return

500 call pmf_utils_exit(PMF_OUT,1,'Unable to read bead id!')

  5 format (' >> String method is disabled!')
 10 format ('fmode                                  = ',i12)
 50 format ('fsample                                = ',i12)
 55 format ('fsample                                = ',i12,'                  (default)')
 60 format ('ftensor                                = ',i12)
 65 format ('ftensor                                = ',i12,'                  (default)')
140 format ('fbeadidfile                            = ',a)
145 format ('fbeadidfile                            = ',a12,'                  (default)')

110 format ('fserverkey                             = ',a)
120 format ('fserver                                = ',a)
130 format ('fpassword                              = ',a16)

150 format ('Reading bead ID from file              = ',a)
155 format ('Bead ID                                = ',I3)

#ifndef PMFLIB_NETWORK
105 format (' >> String method is not compiled in (network support is required)!')
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
