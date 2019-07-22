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

module rst_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  rst_control_read_rst
!===============================================================================

subroutine rst_control_read_rst(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use rst_dat
    use rst_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    call rst_init_dat

    write(PMF_OUT,'(/,a)') '--- [rst] ----------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'rst') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [metadyn] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[RST] fmode item is mandatory in this section')
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

    if(prmfile_get_integer_by_key(prm_fin,'fplevel', fplevel)) then
        write(PMF_OUT,40) fplevel
    else
        write(PMF_OUT,45) fplevel
    end if

    if(prmfile_get_logical_by_key(prm_fin,'frestart', frestart)) then
        write(PMF_OUT,70) prmfile_onoff(frestart)
    else
        write(PMF_OUT,75) prmfile_onoff(frestart)
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fhistupdate', fhistupdate)) then
        write(PMF_OUT,76) fhistupdate
    else
        write(PMF_OUT,77) fhistupdate
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fhistclear', fhistclear)) then
        write(PMF_OUT,86) fhistclear
    else
        write(PMF_OUT,87) fhistclear
    end if

    if(prmfile_get_real8_by_key(prm_fin,'fwarnlevel', fwarnlevel)) then
        call pmf_unit_conv_to_ivalue(EnergyUnit,fwarnlevel)
        write(PMF_OUT,90) pmf_unit_get_rvalue(EnergyUnit,fwarnlevel), pmf_unit_label(EnergyUnit)
    else
        write(PMF_OUT,95) pmf_unit_get_rvalue(EnergyUnit,fwarnlevel), pmf_unit_label(EnergyUnit)
    end if

    rst_enabled = fmode .gt. 0

    return

  5 format (' >> Restrained dynamics is disabled!')
 10 format ('fmode                                  = ',i12)
 30 format ('fsample                                = ',i12)
 35 format ('fsample                                = ',i12,'                  (default)')
 40 format ('fplevel                                = ',i12)
 45 format ('fplevel                                = ',i12,'                  (default)')
 70 format ('frestart                               = ',a12)
 75 format ('frestart                               = ',a12,'                  (default)')
 76 format ('fhistupdate                            = ',i12)
 77 format ('fhistupdate                            = ',i12,'                  (default)')
 86 format ('fhistclear                             = ',i12)
 87 format ('fhistclear                             = ',i12,'                  (default)')
 90 format ('fwarnlevel                             = ',f12.3,1X,A)
 95 format ('fwarnlevel                             = ',f12.3,1X,A15'  (default)')

end subroutine rst_control_read_rst

!===============================================================================
! Subroutine:  rst_control_read_cvs
!===============================================================================

subroutine rst_control_read_cvs(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use rst_dat

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{RST}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( frstdef(1:1) .eq. '{' ) then
        grpname = frstdef(2:len_trim(frstdef)-1)
         write(PMF_OUT,110) grpname
        ! open goup with name from rstdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130)
            rst_enabled = .false.
            return
        end if
        call rst_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(frstdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,frstdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to load file: ' // trim(frstdef) // '!')
        end if

        call rst_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Restraints are read from group: ',A)
120 format('Restraints are read from file : ',A)
130 format(' >> No {RST} group was specified - disabling restrained dynamics!')

end subroutine rst_control_read_cvs

!===============================================================================
! Subroutine:  rst_control_read_cvs_from_group
!===============================================================================

subroutine rst_control_read_cvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common
    use rst_dat
    use rst_restraints

    implicit none
    type(PRMFILE_TYPE),intent(inout)           :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_SECTION_NAME)    :: resname
    character(len=PRMFILE_MAX_LINE)            :: cvname
    integer                                    :: i, j, alloc_failed
    logical                                    :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfRSTItems = prmfile_count_group(prm_fin)

    if( NumOfRSTItems .le. 0 ) then
        ! on CV in current or specified group
        fmode = 0
        rst_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfRSTItems

    ! allocate constraint list ----------------------------------------------------
    allocate(RSTCVList(NumOfRSTItems), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to allocate memory for coordinate data!')
    end if

    do i=1,NumOfRSTItems
        call rst_restraints_reset_rst(RSTCVList(i))
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
                 '[RST] Illegal section name ['//trim(resname)//'] - only [CV] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',cvname)) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] CV name is not provided!')
        end if
        write(PMF_OUT,140) trim(cvname)

        RSTCVList(i)%cvindx = cv_common_find_cv(cvname)
        RSTCVList(i)%cv => CVList(RSTCVList(i)%cvindx)%cv

        ! read the rest of restraint
        call rst_restraints_read_rst(prm_fin,RSTCVList(i))

        eresult = prmfile_next_section(prm_fin)

        i = i + 1
    end do

    ! check if there is CV overlap
    do i=1,NumOfRSTItems
        do j=i+1,NumOfRSTItems
            if( RSTCVList(i)%cvindx .eq. RSTCVList(j)%cvindx ) then
                call pmf_utils_exit(PMF_OUT,1,'[RST] Two different restraints share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No restraints are defined. Umbrella sampling method is switched off!')
110 format('Number of restraints : ',I4)
130 format('== Reading restraint variable #',I4.4)
140 format('   Collective variable name : ',a)

end subroutine rst_control_read_cvs_from_group

!===============================================================================

end module rst_control
