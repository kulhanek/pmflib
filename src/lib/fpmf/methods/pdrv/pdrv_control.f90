!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pdrv_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  pdrv_control_read_pdrv
! load only [pdrv] section
!===============================================================================

subroutine pdrv_control_read_pdrv(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pdrv_dat
    use pdrv_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    call pdrv_init_dat

    write(PMF_OUT,'(/,a)')   '--- [pdrv] ---------------------------------------------------------------------'

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,5)
        return
    end if

    ! try open section
    if( .not. prmfile_open_section(prm_fin,'pdrv') ) then
        write(PMF_OUT,5)
        return
    end if

    ! process options from [pdrv] section
    if( .not. prmfile_get_integer_by_key(prm_fin,'fmode',fmode) ) then
        call pmf_utils_exit(PMF_OUT,1,'[PDRV] fmode item is mandatory in this section')
    else
        write(PMF_OUT,10) fmode
    end if

    if (fmode .ne. 0 .and. fmode .ne. 1) then
        write(PMF_OUT, '(/2x,a,i3,a)') 'fmode (', fmode, ') must be 0 or 1'
        call pmf_utils_exit(PMF_OUT,1)
    end if

    if( fmode .eq. 0 ) then
        write(PMF_OUT,5)
        ! no mon - rest of section is skipped
        call prmfile_set_sec_as_processed(prm_fin)
        return
    end if

    if(prmfile_get_integer_by_key(prm_fin,'fsample', fsample)) then
        write(PMF_OUT,50) fsample
    else
        write(PMF_OUT,55) fsample
    end if

    pdrv_enabled = fmode .gt. 0

    return

  5 format (' >> Path driving is disabled!')
 10 format ('fmode                                  = ',i12)
 50 format ('fsample                                = ',i12)
 55 format ('fsample                                = ',i12,'                  (default)')

end subroutine pdrv_control_read_pdrv

!===============================================================================
! Subroutine:  pdrv_control_read_pdrvs
!===============================================================================

subroutine pdrv_control_read_pdrvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use pdrv_dat
    use prmfile

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_GROUP_NAME)  :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{PDRV}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fpdrvdef(1:1) .eq. '{' ) then
        grpname = fpdrvdef(2:len_trim(fpdrvdef)-1)
         write(PMF_OUT,110) trim(grpname)
        ! open goup with name from mondef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            call pmf_utils_exit(PMF_OUT,1,'[PDRV] Unable to open group {' // trim(grpname) // '}!')
        end if
        call pdrv_control_read_pdrvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fpdrvdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fpdrvdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[PDRV] Unable to load file: ' // trim(fpdrvdef) // '!')
        end if

        call pdrv_control_read_pdrvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Driven paths are read from group: ',A)
120 format('Driven paths are read from file : ',A)

end subroutine pdrv_control_read_pdrvs

!===============================================================================
! Subroutine:  pdrv_control_read_pdrvs_from_group
!===============================================================================

subroutine pdrv_control_read_pdrvs_from_group(prm_fin)

    use prmfile
    use pmf_dat
    use pdrv_dat
    use pmf_utils
    use pdrv_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)        :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_SECTION_NAME) :: resname
    character(len=PRMFILE_MAX_LINE)         :: pathname
    integer                                 :: i, j, alloc_failed
    logical                                 :: eresult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfPDRVItems = prmfile_count_group(prm_fin)

    if( NumOfPDRVItems .le. 0 ) then
        ! no PATH in current or specified group
        fmode = 0
        pdrv_enabled = .false.
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfPDRVItems

    ! allocate list of CVs indexes ------------------------------------------------
    allocate(PDRVCVList(NumOfPDRVItems), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PDRV] Unable to allocate memory for driven paths!')
    end if

    ! enumerate sections ----------------------------------------------------------
    eresult = prmfile_first_section(prm_fin)
    i = 1
    do while(eresult)
        eresult = prmfile_get_section_name(prm_fin,resname)
        write(PMF_OUT,*)
        write(PMF_OUT,130) i
        if( resname .ne. 'PATH' ) then
            call pmf_utils_exit(PMF_OUT,1, &
                 '[PDRV] Illegal section name ['//trim(resname)//'] - only [PATH] is allowed!')
        end if
        if( .not. prmfile_get_string_by_key(prm_fin,'name',pathname)) then
            call pmf_utils_exit(PMF_OUT,1,'[PDRV] PATH name is not provided!')
        end if
        write(PMF_OUT,140) trim(pathname)

        PDRVCVList(i)%pathindx = pmf_paths_find_path(pathname)
        PDRVCVList(i)%path => PathList(PDRVCVList(i)%pathindx)%path
        PDRVCVList(i)%path%driven_mode = .true.

        call pdrv_paths_read_pdrv(prm_fin,PDRVCVList(i))

        eresult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    ! check if there is PATH overlap
    do i=1,NumOfPDRVItems
        do j=i+1,NumOfPDRVItems
            if( PDRVCVList(i)%pathindx .eq. PDRVCVList(j)%pathindx ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     '[MTD] Two different STM collective variables share the same general collective variable!')
            end if
        end do
    end do

    return

100 format('>>> INFO: No driven paths are defined. The path driving is switched off!')
110 format('Number of driven paths : ',I2)
130 format('== Reading driven path #',I2.2)
140 format('   Driven path name : ',a)

end subroutine pdrv_control_read_pdrvs_from_group

!===============================================================================

end module pdrv_control
