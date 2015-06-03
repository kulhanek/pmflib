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

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call mtd_init_dat

    write(PMF_OUT,'(/,a)')   '--- [mtd] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'mtd') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [mtd] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] fmode item is mandatory in this section')
    else
        write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1 .and. fmode .ne. 2 .and. fmode .ne. 3) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0, 1, 2 or 3'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no metadynamics - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    if(prmfile_get_real8_by_key(prm_fin,'fheight', fheight) ) then
        write(PMF_OUT,30) fheight
    else
        write(PMF_OUT,35) fheight
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fmetastep', fmetastep)) then
        write(PMF_OUT,40) fmetastep
    else
        write(PMF_OUT,45) fmetastep
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fpsample', fpsample)) then
        write(PMF_OUT,50) fpsample
    else
        write(PMF_OUT,55) fpsample
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fplevel', fplevel)) then
        write(PMF_OUT,60) fplevel
    else
        write(PMF_OUT,65) fplevel
    end if

    if(prmfile_get_logical_by_key(prm_fin,'frestart', frestart)) then
        write(PMF_OUT,70) prmfile_onoff(frestart)
    else
        write(PMF_OUT,75) prmfile_onoff(frestart)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fextout', fextout)) then
        write(PMF_OUT,80) fextout
    else
        write(PMF_OUT,85) fextout
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fbuffersize', fbuffersize)) then
        write(PMF_OUT,90) fbuffersize
    else
        write(PMF_OUT,95) fbuffersize
    end if

    if( fbuffersize .le. 0 ) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fbuffersize (', fbuffersize, ') must be grater then zero.'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_real8_by_key(prm_fin,'fmetatemp', fmetatemp)) then
        fmetavary = 1
        write(PMF_OUT,100) fmetatemp
    else
        write(PMF_OUT,105) fmetatemp
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fmetavary', fmetavary)) then
        meta_next_fstep = fmetastep
        write(PMF_OUT,110) fmetavary
    else
        write(PMF_OUT,115) fmetavary
    end if

    if( (fmetavary .lt. 0) .or. (fmetavary .gt. 2) ) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmetavary (', fmetavary, ') must be 0, 1 or 2.'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fscaling', fscaling)) then
        write(PMF_OUT,120) fscaling
    else
        write(PMF_OUT,125) fscaling
    end if

    if( (fscaling .lt. 0) .or. (fscaling .gt. 2) ) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fscaling (', fscaling, ') must be 0, 1 or 2.'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fdelaying', fdelaying)) then
        ftransition = fdelaying
        write(PMF_OUT,130) fdelaying
    else
        write(PMF_OUT,135) fdelaying
    end if

    if( fdelaying .lt. 0 ) then
        write(PMF_OUT, '(/2x,a,i5,a)') 'fdelaying (', fdelaying, ') must be grater than zero'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( (fscaling .eq. 0) .and. (fdelaying .ne. 0) ) then
        write(PMF_OUT, '(/2x,a)') 'fscaling must be also specified when fdelaying is specified'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'ftransition', ftransition)) then
        write(PMF_OUT,140) ftransition
    else
        write(PMF_OUT,145) ftransition
    end if 

    if( ftransition .lt. 0 ) then
        write(PMF_OUT, '(/2x,a,i5,a)') 'ftransition (', ftransition, ') must be grater than zero'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( (fscaling .eq. 0) .and. (ftransition .ne. 0) ) then
        write(PMF_OUT, '(/2x,a)') 'fscaling must be also specified when ftransition is specified'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( ftransition .lt. fdelaying ) then
        write(PMF_OUT, '(/2x,a)') 'ftransition must be grater or equal than fdelaying'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 3 ) then
            write(PMF_OUT,410)
            if( .not. prmfile_get_integer_by_key(prm_fin,'fgpcovfreq', fgpcovfreq)) then
               call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpcovfreq item is mandatory when GP-MTD is on!')
            else
                write(PMF_OUT,415) fgpcovfreq
            end if
            if( fgpcovfreq .le. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpcovfreq must be > 0!')
            end if
            fbuffersize = fgpcovfreq / fmetastep
            if( .not. prmfile_get_integer_by_key(prm_fin,'fgpsparse', fgpsparse)) then
                call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpsparse item is mandatory when GP-MTD is on!')
            else
                write(PMF_OUT,425) fgpsparse
            end if
            if(prmfile_get_real8_by_key(prm_fin,'fgpdelta', fgpdelta)) then
                write(PMF_OUT,430) fgpdelta
            else
                write(PMF_OUT,435) fgpdelta
            end if
            if( fgpdelta .le. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpdelta must be > 0!')
            end if
            if(prmfile_get_real8_by_key(prm_fin,'fgpjitter', fgpjitter)) then
                write(PMF_OUT,440) fgpjitter
            else
                write(PMF_OUT,445) fgpjitter
            end if
            if( fgpjitter .le. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpjitter must be > 0!')
            end if
            if(prmfile_get_integer_by_key(prm_fin,'fgpsparsification', fgpsparsification)) then
                write(PMF_OUT,450) fgpsparsification
            else
                write(PMF_OUT,455) fgpsparsification
            end if
            if ( (fgpsparsification .le. 0) .or. (fgpsparsification .gt. 2) ) then
                  call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpsparsification must be 1 or 2!')
            end if
            if(prmfile_get_integer_by_key(prm_fin,'fgpsparsefreq', fgpsparsefreq)) then
                  write(PMF_OUT,460) fgpsparsefreq
            else
                  fgpsparsefreq = fgpcovfreq
                  write(PMF_OUT,465) fgpsparsefreq
            end if
            if( fgpsparsefreq .le. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpsparsefreq must be > 0!')
            end if
            select case (fgpsparsification)
                case(1)
                    if( .not. prmfile_get_integer_by_key(prm_fin,'fgpclusterfreq', fgpclusterfreq)) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpclusterfreq item is mandatory when fgpclustering = 1!')
                    else
                        write(PMF_OUT,470) fgpclusterfreq
                    end if
                    if( fgpclusterfreq .le. 0 ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] fgpclusterevery must be > 0!')
                    end if
            end select
            if(prmfile_get_integer_by_key(prm_fin,'fgpprint', fgpprint)) then
                write(PMF_OUT,480) fgpprint
            else
                write(PMF_OUT,485) fgpprint
            end if
            if( fgpprint .lt. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[ABF] fgpprint must be => 0!')
            end if
    end if

    ! network setup ----------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '--- [mtd-walker] -----------------------------------------------------------'

    ! try open section
    if( prmfile_open_section(prm_fin,'mtd-walker') ) then

#ifdef PMFLIB_NETWORK
    if( prmfile_get_string_by_key(prm_fin,'fserverkey',fserverkey)) then
        write(PMF_OUT,490) trim(fserver)
        fserver_enabled = .true.
        fserver_key_enabled = .true.
    end if

    if( .not. fserver_enabled ) then
        if( prmfile_get_string_by_key(prm_fin,'fserver', fserver) ) then
            write(PMF_OUT,500) fserver
            fserver_enabled = .true.
            fserver_key_enabled = .false.
        else
            call pmf_utils_exit(PMF_OUT,1,'fserver is required when [abf-walker] is specified')
        end if
        if( prmfile_get_string_by_key(prm_fin,'fpassword', fpassword) ) then
            write(PMF_OUT,510) fpassword
        else
            call pmf_utils_exit(PMF_OUT,1,'fpassword is required when [abf-walker] is specified')
        end if
    end if

#else
    fserver_enabled = .false.
    fserver_key_enabled = .false.
    write(PMF_OUT,7)
#endif
    ! network setup ----------------------------------------------------------------

    else
        write(PMF_OUT,6)
    end if

    ! restart is read from server
    if( fserver_enabled .and. frestart ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] frestart cannot be on if multiple-walker aproach is used!')
    end if

    if( fserver_enabled .and. fextout .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] fextout has to zero if multiple-walker aproach is used!')
    end if
    mtd_enabled = fmode .gt. 0

    return

  5 format (' >> Metadynamics is disabled!')
 10 format ('fmode                                  = ',i12)
 15 format ('fmode                                  = ',i12,'                  (default)')
 30 format ('fheight                                = ',E12.4)
 35 format ('fheight                                = ',E12.4,'                  (default)')
 40 format ('fmetastep                              = ',i12)
 45 format ('fmetastep                              = ',i12,'                  (default)')
 50 format ('fsample                                = ',i12)
 55 format ('fsample                                = ',i12,'                  (default)')
 60 format ('fplevel                                = ',i12)
 65 format ('fplevel                                = ',i12,'                  (default)')
 70 format ('frestart                               = ',a12)
 75 format ('frestart                               = ',a12,'                  (default)')
 80 format ('fextout                                = ',i12)
 85 format ('fextout                                = ',i12,'                  (default)')
 90 format ('fbuffersize                            = ',i12)
 95 format ('fbuffersize                            = ',i12,'                  (default)')
100 format ('fmetatemp                              = ',E12.4)
105 format ('fmetatemp                              = ',E12.4,'                  (default)')
110 format ('fmetavary                              = ',i12)
115 format ('fmetavary                              = ',i12,'                  (default)')
120 format ('fscaling                               = ',i12)
125 format ('fscaling                               = ',i12,'                  (default)')
130 format ('fdelaying                              = ',i12)
135 format ('fdelaying                              = ',i12,'                  (default)')
140 format ('ftransition                            = ',i12)
145 format ('ftransition                            = ',i12,'                  (default)')

410 format (/,'>> GP-ABF mode IV')
415 format ('   fgpcovfreq                          = ',i12)
425 format ('   fgpsparse                           = ',i12)
430 format ('   fgpdelta                            = ',e12.4)
435 format ('   fgpdelta                            = ',e12.4,'                  (default)')
440 format ('   fgpjitter                           = ',e12.4)
445 format ('   fgpjitter                           = ',e12.4,'                  (default)')
450 format ('   fgpsparsification                   = ',i12)
455 format ('   fgpsparsification                   = ',i12,'                  (default)')
460 format ('   fgpsparsefreq                       = ',i12)
465 format ('   fgpsparsefreq                       = ',i12,'                  (default)')
470 format ('   fgpclusterfreq                      = ',i12)
480 format ('   fgpprint                            = ',i12)
485 format ('   fgpprint                            = ',i12,'                  (default)')

  6 format (' >> Multiple-walkers metadynamics is disabled!')
  7 format (' >> Multiple-walkers metadynamics is not compiled in!')
490 format ('fserverkey                             = ',a)
500 format ('fserver                                = ',a)
510 format ('fpassword                              = ',a12)

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
         write(PMF_OUT,110) grpname
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
    use mtd_cvs

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
