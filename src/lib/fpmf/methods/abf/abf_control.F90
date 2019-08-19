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

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call abf_init_dat

    write(PMF_OUT,'(/,a)') '--- [abf] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'abf') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [abf] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] fmode item is mandatory in this section')
    else
     write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1 .and. fmode .ne. 2) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0, 1 or 2'
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

    if(prmfile_get_integer_by_key(prm_fin,'feimode', feimode)) then
        write(PMF_OUT,20) feimode
    else
        write(PMF_OUT,25) feimode
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fsample', fsample)) then
        write(PMF_OUT,50) fsample
    else
        write(PMF_OUT,55) fsample
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

    if(prmfile_get_logical_by_key(prm_fin,'fprint_icf', fprint_icf)) then
        write(PMF_OUT,90) prmfile_onoff(fprint_icf)
    else
        write(PMF_OUT,95) prmfile_onoff(fprint_icf)
    end if

    if(prmfile_get_logical_by_key(prm_fin,'fcache_icf', fcache_icf)) then
        write(PMF_OUT,400) prmfile_onoff(fcache_icf)
    else
        write(PMF_OUT,405) prmfile_onoff(fcache_icf)
    end if

    if(prmfile_get_logical_by_key(prm_fin,'frawicf', frawicf)) then
        write(PMF_OUT,410) prmfile_onoff(frawicf)
    else
        write(PMF_OUT,415) prmfile_onoff(frawicf)
    end if

    select case(feimode)
        case(1)
            write(PMF_OUT,190)
            if(prmfile_get_integer_by_key(prm_fin,'fhramp', fhramp)) then
                write(PMF_OUT,192) fhramp
            else
                write(PMF_OUT,195) fhramp
            end if
            if( fhramp .le. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[ABF] fhramp must be > 0!')
            end if
        case(2)
            write(PMF_OUT,200)
            if(prmfile_get_integer_by_key(prm_fin,'fhramp_min', fhramp_min)) then
                write(PMF_OUT,202) fhramp_min
            else
                write(PMF_OUT,205) fhramp_min
            end if
            if(prmfile_get_integer_by_key(prm_fin,'fhramp_max', fhramp_max)) then
                write(PMF_OUT,206) fhramp_max
            else
                write(PMF_OUT,207) fhramp_max
            end if
            if( fhramp_min .le. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[ABF] fhramp_min must be > 0!')
            end if
            if( fhramp_max .le. fhramp_min ) then
                call pmf_utils_exit(PMF_OUT,1,'[ABF] fhramp_max must be > fhramp_min!')
            end if
        case(3)
            write(PMF_OUT,300)
            if(prmfile_get_integer_by_key(prm_fin,'fgpmin_samples', fgpmin_samples)) then
                write(PMF_OUT,310) fgpmin_samples
            else
                write(PMF_OUT,315) fgpmin_samples
            end if
            if(prmfile_get_integer_by_key(prm_fin,'fgpmodel_update', fgpmodel_update)) then
                write(PMF_OUT,320) fgpmodel_update
            else
                write(PMF_OUT,325) fgpmodel_update
            end if
            if(prmfile_get_integer_by_key(prm_fin,'fgpprint_period', fgpprint_period)) then
                write(PMF_OUT,330) fgpprint_period
            else
                write(PMF_OUT,335) fgpprint_period
            end if
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
        call pmf_utils_exit(PMF_OUT,1,'[ABF] frestart cannot be ON if multiple-walkers aproach is used!')
    end if

    abf_enabled = fmode .gt. 0

    return

  5 format (' >> Adaptive biasing force method is disabled!')
 10 format ('fmode                                  = ',i12)
 16 format ('fmask_mode                             = ',i12)
 19 format ('fmask_mode                             = ',i12,'                  (default)')
 21 format ('fapply_abf                             = ',a12)
 26 format ('fapply_abf                             = ',a12,'                  (default)')
 20 format ('feimode                                = ',i12)
 25 format ('feimode                                = ',i12,'                  (default)')
 50 format ('fsample                                = ',i12)
 55 format ('fsample                                = ',i12,'                  (default)')
 70 format ('frestart                               = ',a12)
 75 format ('frestart                               = ',a12,'                  (default)')
 76 format ('frstupdate                             = ',i12)
 77 format ('frstupdate                             = ',i12,'                  (default)')
 80 format ('ftrjsample                             = ',i12)
 85 format ('ftrjsample                             = ',i12,'                  (default)')
 90 format ('fprint_icf                             = ',a12)
 95 format ('fprint_icf                             = ',a12,'                  (default)')
400 format ('fcache_icf                             = ',a12)
405 format ('fcache_icf                             = ',a12,'                  (default)')
410 format ('frawicf                                = ',a12)
415 format ('frawicf                                = ',a12,'                  (default)')

100 format (' >> Multiple-walkers ABF method is disabled!')
#ifndef PMFLIB_NETWORK
105 format (' >> Multiple-walkers ABF method is not compiled in!')
#endif
110 format ('fserverkey                             = ',a)
120 format ('fserver                                = ',a)
130 format ('fpassword                              = ',a16)
140 format ('fserverupdate                          = ',i12)
145 format ('fserverupdate                          = ',i12,'                  (default)')
150 format ('fabortonmwaerr                          = ',a12)
155 format ('fabortonmwaerr                          = ',a12,'                  (default)')

190 format (/,'>> Linear ramp mode I (feimode == 1)')
192 format ('   fhramp                              = ',i12)
195 format ('   fhramp                              = ',i12,'                  (default)')

200 format (/,'>> Linear ramp mode II (feimode == 2)')
202 format ('   fhramp_min                          = ',i12)
205 format ('   fhramp_min                          = ',i12,'                  (default)')
206 format ('   fhramp_max                          = ',i12)
207 format ('   fhramp_max                          = ',i12,'                  (default)')

300 format (/,'>> Gaussian process (feimode == 3)')
310 format ('   fgpmin_samples                      = ',i12)
315 format ('   fgpmin_samples                      = ',i12,'                  (default)')
320 format ('   fgpmodel_update                     = ',i12)
325 format ('   fgpmodel_update                     = ',i12,'                  (default)')
330 format ('   fgpprint_period                     = ',i12)
335 format ('   fgpprint_period                     = ',i12,'                  (default)')

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
