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

module con_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  con_control_read_con
!===============================================================================

subroutine con_control_read_con(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use con_dat
    use con_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call con_init_dat

    write(PMF_OUT,'(/,a)') '--- [con] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'con') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [metadyn] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[CON] fmode item is mandatory in this section')
    else
        write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0 or 1!'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no umbrealla - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fsample', fsample)) then
        write(PMF_OUT,30) fsample
    else
        write(PMF_OUT,35) fsample
    end if

    if(prmfile_get_integer_by_key(prm_fin,'faccurst', faccurst)) then
        write(PMF_OUT,40) faccurst
    else
        write(PMF_OUT,45) faccurst
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fplevel', fplevel)) then
        write(PMF_OUT,50) fplevel
    else
        write(PMF_OUT,55) fplevel
    end if

    if(prmfile_get_logical_by_key(prm_fin,'frestart', frestart)) then
        write(PMF_OUT,70) prmfile_onoff(frestart)
    else
        write(PMF_OUT,75) prmfile_onoff(frestart)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'flambdasolver', flambdasolver) ) then
        write(PMF_OUT,76) flambdasolver, trim(con_init_get_lsolver_name(flambdasolver))
    else
        write(PMF_OUT,77) flambdasolver, trim(con_init_get_lsolver_name(flambdasolver))
    end if

    if(prmfile_get_real8_by_key(prm_fin,'flambdatol', flambdatol) ) then
        write(PMF_OUT,80) flambdatol
    else
        write(PMF_OUT,85) flambdatol
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fmaxiter', fmaxiter)) then
        write(PMF_OUT,90) fmaxiter
    else
        write(PMF_OUT,95) fmaxiter
    end if

    ! check input parameters ------------------------------------------------------
    if( flambdasolver .ne. 0 .and. flambdasolver .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'[CON] flambdasolver has to be 0 or 1!')
    end if

    con_enabled = fmode .gt. 0

    return

  5 format (' >> Constrained dynamics is disabled!')
 10 format ('fmode         = ',i12)
 15 format ('fmode         = ',i12,'                                           (default)')
 30 format ('fsample       = ',i12)
 35 format ('fsample       = ',i12,'                                           (default)')
 40 format ('faccurst      = ',i12)
 45 format ('faccurst      = ',i12,'                                           (default)')
 50 format ('fplevel       = ',i12)
 55 format ('fplevel       = ',i12,'                                           (default)')
 70 format ('frestart      = ',a12)
 75 format ('frestart      = ',a12,'                                           (default)')
 76 format ('flambdasolver = ',i12,1x,a)
 77 format ('flambdasolver = ',i12,1x,a25,'                 (default)')
 80 format ('flambdatol    = ',E12.4)
 85 format ('flambdatol    = ',E12.4,'                                           (default)')
 90 format ('fmaxiter      = ',i12)
 95 format ('fmaxiter      = ',i12,'                                           (default)')

end subroutine con_control_read_con

!===============================================================================
! Subroutine:  con_control_read_cvs
!===============================================================================

subroutine con_control_read_cvs(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use con_dat
    use con_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprm_fin
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{CON}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fcondef(1:1) .eq. '{' ) then
        grpname = fcondef(2:len_trim(fcondef)-1)
         write(PMF_OUT,110) grpname
        ! open goup with name from metadef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            con_enabled = .false.
            return
        end if
        call con_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fcondef)

        call prmfile_init(locprm_fin)

        if( .not. prmfile_read(locprm_fin,fcondef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[CON] Unable to load file: ' // trim(fcondef) // '!')
        end if

        call con_control_read_cvs_from_group(locprm_fin)

        call prmfile_clear(locprm_fin)
    end if

    ! init constrained atom arrays for CON
    call con_init_con_atoms

    return

110 format('Constraints are read from group: ',A)
120 format('Constraints are read from file : ',A)
130 format(' >> No {CON} group was specified - disabling constrained dynamics!')

end subroutine con_control_read_cvs

!===============================================================================
! Subroutine:  con_control_read_cvs_from_group
!===============================================================================

subroutine con_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use con_dat
    use con_constraints

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
        con_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfCONs

    ! allocate constraint list ----------------------
    allocate(CONList(NumOfCONs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[CON] Unable to allocate memory for constraint data!')
    end if

    do i=1,NumOfCONs
        call con_constraints_reset_con(CONList(i))
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
                 '[CON] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[CON] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        CONList(i)%cvindx = cv_common_find_cv(cvname)
        CONList(i)%cv     => CVList(CONList(i)%cvindx)%cv

        call con_constraints_read_con(prm_fin,CONList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfCONs
        do j=i+1,NumOfCONs
            if( CONList(i)%cvindx .eq. CONList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1,'[CON] Two different constraints share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No constraints are defined. Blue moon method is switched off!')
110 format('Number of constraints : ',I4)
130 format('== Reading constraint #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine con_control_read_cvs_from_group

!===============================================================================

end module con_control

