!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
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

module pmf_control

implicit none
contains

!===============================================================================
! Subroutine:  pmf_control_read_pmflib_group
!===============================================================================

subroutine pmf_control_read_pmflib_group(prm_fin)

    use prmfile
    use pmf_constants
    use pmf_utils
    use mon_control
    use rst_control
    use abf_control
    use tabf_control
    use usabf_control
    use abp_control
    use mtd_control
    use cst_control
    use stm_control
    use pdrv_control
    use pmf_dat
    use pmf_mask
    use pmf_init

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    ! read PMFLib setup from {PMFLIB} group ---------
    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{PMFLIB}',':')

    ! try open group
    if( .not. prmfile_open_group(prm_fin,'PMFLIB') ) then
        write(PMF_OUT,110)
        return
    end if

    call pmf_control_read_control(prm_fin)
    call pmf_control_read_units(prm_fin)

    ! load fake topology - if it is provided
    if( (len_trim(ftopology) .gt. 0) .and. (ftopology .ne. '-system-') ) then
        call pmf_mask_topo_load(ftopology)
    end if

    ! read method setup
    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Methods','=')

    call abf_control_read_abf(prm_fin)
    call abp_control_read_abp(prm_fin)
    call mtd_control_read_mtd(prm_fin)
    call stm_control_read_stm(prm_fin)
    call cst_control_read_con(prm_fin)
    call rst_control_read_rst(prm_fin)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Methods (Testing)','=')
    call tabf_control_read_abf(prm_fin)
    call usabf_control_read_abf(prm_fin)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Extensions','=')
    call pdrv_control_read_pdrv(prm_fin)
    call mon_control_read_mon(prm_fin)

    ! read filenames
    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Files','=')
    call pmf_control_read_files(prm_fin)

    ! print initial topology
    if( fprint_inpcrds ) then
        write(PMF_OUT,*)
        call pmf_utils_heading(PMF_OUT,'PMFLIB TOPOLOGY','=')
        write(PMF_OUT,*)
        call pmf_mask_topo_print
    end if

    ! load collective variables
    call pmf_control_read_cvs(prm_fin)

    ! load paths
    call pmf_control_read_paths(prm_fin)

    return

110 format (' >> No PMFLIB group was found - PMFLib is disabled!')

end subroutine pmf_control_read_pmflib_group

!===============================================================================
! Subroutine:  pmf_control_read_control
!===============================================================================

subroutine pmf_control_read_control(prm_fin)

    use pmf_constants
    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '--- [control] ------------------------------------------------------------------'

    ftopology = '-system-'

    if(.not. prmfile_open_section(prm_fin,'control')) then
        call pmf_ctrl_print_default_stritem('ftopology',ftopology)
        call pmf_ctrl_print_default_logical('fdebug',fdebug)
        call pmf_ctrl_print_default_logical('frepmpifrag',frepmpifrag)
        call pmf_ctrl_print_default_logical('fprint_inpcrds',fprint_inpcrds)
        call pmf_ctrl_print_default_logical('fprint_masks',fprint_masks)
        call pmf_ctrl_print_default_logical('fenable_pbc',fenable_pbc)
        return
    end if

    call pmf_ctrl_read_stritem(prm_fin,'ftopology',ftopology)
    call pmf_ctrl_read_logical(prm_fin,'fdebug',fdebug)
    call pmf_ctrl_read_logical(prm_fin,'frepmpifrag',frepmpifrag)
    call pmf_ctrl_read_logical(prm_fin,'fprint_inpcrds',fprint_inpcrds)
    call pmf_ctrl_read_logical(prm_fin,'fprint_masks',fprint_masks)
    call pmf_ctrl_read_logical(prm_fin,'fenable_pbc',fenable_pbc)

end subroutine pmf_control_read_control

!===============================================================================
! Subroutine:  pmf_control_read_units
!===============================================================================

subroutine pmf_control_read_units(prm_fin)

    use pmf_constants
    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_unit

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! -----------------------------------------------
    character(PMF_MAX_SUNIT)           :: buffer
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') &
          '--- [units] --------------------------------------------------------------------'

    if(.not. prmfile_open_section(prm_fin,'units')) then
        write(PMF_OUT,15) pmf_unit_label(MassUnit)
        write(PMF_OUT,25) pmf_unit_label(TimeUnit)
        write(PMF_OUT,35) pmf_unit_label(LengthUnit)
        write(PMF_OUT,45) pmf_unit_label(AngleUnit)
        write(PMF_OUT,55) pmf_unit_label(EnergyUnit)
        write(PMF_OUT,65) pmf_unit_label(TemperatureUnit)
        return
    end if

    if(.not. prmfile_get_string_by_key(prm_fin,'mass', buffer)) then
        write(PMF_OUT,15) pmf_unit_label(MassUnit)
    else
        call pmf_unit_decode_mass_unit(buffer,MassUnit)
        write(PMF_OUT,10) pmf_unit_label(MassUnit)
    end if

    if(.not. prmfile_get_string_by_key(prm_fin,'time', buffer)) then
        write(PMF_OUT,25) pmf_unit_label(TimeUnit)
    else
        call pmf_unit_decode_time_unit(buffer,TimeUnit)
        write(PMF_OUT,20) pmf_unit_label(TimeUnit)
    end if

    if(.not. prmfile_get_string_by_key(prm_fin,'length', buffer)) then
        write(PMF_OUT,35) pmf_unit_label(LengthUnit)
    else
        call pmf_unit_decode_length_unit(buffer,LengthUnit)
        write(PMF_OUT,30) pmf_unit_label(LengthUnit)
    end if

    if(.not. prmfile_get_string_by_key(prm_fin,'angle', buffer)) then
        write(PMF_OUT,45) pmf_unit_label(AngleUnit)
    else
        call pmf_unit_decode_angle_unit(buffer,AngleUnit)
        write(PMF_OUT,40) pmf_unit_label(AngleUnit)
    end if

    if(.not. prmfile_get_string_by_key(prm_fin,'energy', buffer)) then
        write(PMF_OUT,55) pmf_unit_label(EnergyUnit)
    else
        call pmf_unit_decode_energy_unit(buffer,EnergyUnit)
        write(PMF_OUT,50) pmf_unit_label(EnergyUnit)
    end if

    if(.not. prmfile_get_string_by_key(prm_fin,'temperature', buffer)) then
        write(PMF_OUT,65) pmf_unit_label(TemperatureUnit)
    else
        call pmf_unit_decode_temp_unit(buffer,TemperatureUnit)
        write(PMF_OUT,60) pmf_unit_label(TemperatureUnit)
    end if

 10 format('mass                                   = ',a)
 15 format('mass                                   = ',a29,' (default)')
 20 format('time                                   = ',a)
 25 format('time                                   = ',a29,' (default)')
 30 format('length                                 = ',a)
 35 format('length                                 = ',a29,' (default)')
 40 format('angle                                  = ',a)
 45 format('angle                                  = ',a29,' (default)')
 50 format('energy                                 = ',a)
 55 format('energy                                 = ',a29,' (default)')
 60 format('temperature                            = ',a)
 65 format('temperature                            = ',a29,' (default)')

end subroutine pmf_control_read_units

!===============================================================================
! Subroutine:  pmf_control_read_files
!===============================================================================

subroutine pmf_control_read_files(prm_fin)

    use pmf_constants
    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') &
          '--- [files] --------------------------------------------------------------------'

    if( .not. prmfile_open_section(prm_fin,'files') ) then
        write(PMF_OUT,10)
        call pmf_ctrl_print_default_stritem('fcvsdef',fcvsdef)
        call pmf_ctrl_print_default_stritem('fpathsdef',fpathsdef)
        if( abf_enabled ) then
            write(PMF_OUT,400)
            call pmf_ctrl_print_default_stritem('fabfdef',fabfdef)
            call pmf_ctrl_print_default_stritem('fabfmask',fabfmask)
            call pmf_ctrl_print_default_stritem('fabfout',fabfout)
            call pmf_ctrl_print_default_stritem('fabfrst',fabfrst)
        end if
        if( tabf_enabled ) then
            write(PMF_OUT,410)
            call pmf_ctrl_print_default_stritem('ftabfdef',ftabfdef)
            call pmf_ctrl_print_default_stritem('ftabfout',ftabfout)
            call pmf_ctrl_print_default_stritem('ftabfrst',ftabfrst)
            call pmf_ctrl_print_default_stritem('ftabficf',ftabficf)
        end if
        if( usabf_enabled ) then
            write(PMF_OUT,420)
            call pmf_ctrl_print_default_stritem('fusabfdef',fusabfdef)
            call pmf_ctrl_print_default_stritem('fusabfout',fusabfout)
            call pmf_ctrl_print_default_stritem('fusabfrst',fusabfrst)
        end if
        if( mtd_enabled ) then
            write(PMF_OUT,300)
            call pmf_ctrl_print_default_stritem('fmtddef',fmtddef)
            call pmf_ctrl_print_default_stritem('fmtdout',fmtdout)
            call pmf_ctrl_print_default_stritem('fmtdrst',fmtdrst)
            call pmf_ctrl_print_default_stritem('fmtdtrj',fmtdtrj)
            call pmf_ctrl_print_default_stritem('fmtdhills',fmtdhills)
        end if
        if( stm_enabled ) then
            write(PMF_OUT,700)
            call pmf_ctrl_print_default_stritem('fstmdef',fstmdef)
            call pmf_ctrl_print_default_stritem('fstmout',fstmout)
        end if
        if( cst_enabled ) then
            write(PMF_OUT,100)
            call pmf_ctrl_print_default_stritem('fcstdef',fcstdef)
            call pmf_ctrl_print_default_stritem('fcstout',fcstout)
            call pmf_ctrl_print_default_stritem('fcstrst',fcstrst)
            call pmf_ctrl_print_default_stritem('fcstfrst',fcstfrst)
            call pmf_ctrl_print_default_stritem('fcsttrj',fcsttrj)
        end if
        if( rst_enabled ) then
            write(PMF_OUT,200)
            call pmf_ctrl_print_default_stritem('frstdef',frstdef)
            call pmf_ctrl_print_default_stritem('frstout',frstout)
            call pmf_ctrl_print_default_stritem('frsthist',frsthist)
        end if
        if( pdrv_enabled ) then
            write(PMF_OUT,800)
            call pmf_ctrl_print_default_stritem('fpdrvdef',fpdrvdef)
            call pmf_ctrl_print_default_stritem('fpdrvout',fpdrvout)
        end if
        if( mon_enabled ) then
            write(PMF_OUT,500)
            call pmf_ctrl_print_default_stritem('fmondef',fmondef)
            call pmf_ctrl_print_default_stritem('fmonout',fmonout)
        end if
        return
    end if

    write(PMF_OUT,10)
    call pmf_ctrl_print_default_stritem('fcvsdef',fcvsdef)
    call  pmf_ctrl_read_stritem(prm_fin,'fpathsdef',fpathsdef)
    if( abf_enabled ) then
        write(PMF_OUT,400)
        call  pmf_ctrl_read_stritem(prm_fin,'fabfdef',fabfdef)
        call  pmf_ctrl_read_stritem(prm_fin,'fabfmask',fabfmask)
        call  pmf_ctrl_read_stritem(prm_fin,'fabfout',fabfout)
        call  pmf_ctrl_read_stritem(prm_fin,'fabfrst',fabfrst)
    end if
    if( tabf_enabled ) then
        write(PMF_OUT,410)
        call  pmf_ctrl_read_stritem(prm_fin,'ftabfdef',ftabfdef)
        call  pmf_ctrl_read_stritem(prm_fin,'ftabfout',ftabfout)
        call  pmf_ctrl_read_stritem(prm_fin,'ftabfrst',ftabfrst)
        call  pmf_ctrl_read_stritem(prm_fin,'ftabficf',ftabficf)
    end if
    if( usabf_enabled ) then
        write(PMF_OUT,420)
        call  pmf_ctrl_read_stritem(prm_fin,'fusabfdef',fusabfdef)
        call  pmf_ctrl_read_stritem(prm_fin,'fusabfout',fusabfout)
        call  pmf_ctrl_read_stritem(prm_fin,'fusabfrst',fusabfrst)
    end if
    if( mtd_enabled ) then
        write(PMF_OUT,300)
        call  pmf_ctrl_read_stritem(prm_fin,'fmtddef',fmtddef)
        call  pmf_ctrl_read_stritem(prm_fin,'fmtdout',fmtdout)
        call  pmf_ctrl_read_stritem(prm_fin,'fmtdrst',fmtdrst)
        call  pmf_ctrl_read_stritem(prm_fin,'fmtdtrj',fmtdtrj)
        call  pmf_ctrl_read_stritem(prm_fin,'fmtdhills',fmtdhills)
    end if
    if( stm_enabled ) then
        write(PMF_OUT,700)
        call  pmf_ctrl_read_stritem(prm_fin,'fstmdef',fstmdef)
        call  pmf_ctrl_read_stritem(prm_fin,'fstmout',fstmout)
    end if
    if( cst_enabled ) then
        write(PMF_OUT,100)
        call  pmf_ctrl_read_stritem(prm_fin,'fcstdef',fcstdef)
        call  pmf_ctrl_read_stritem(prm_fin,'fcstout',fcstout)
        call  pmf_ctrl_read_stritem(prm_fin,'fcstrst',fcstrst)
        call  pmf_ctrl_read_stritem(prm_fin,'fcstfrst',fcstfrst)
        call  pmf_ctrl_read_stritem(prm_fin,'fcsttrj',fcsttrj)
    end if
    if( rst_enabled ) then
        write(PMF_OUT,200)
        call  pmf_ctrl_read_stritem(prm_fin,'frstdef',frstdef)
        call  pmf_ctrl_read_stritem(prm_fin,'frstout',frstout)
        call  pmf_ctrl_read_stritem(prm_fin,'frsthist',frsthist)
    end if
    if( pdrv_enabled ) then
        write(PMF_OUT,800)
        call  pmf_ctrl_read_stritem(prm_fin,'fpdrvdef',fpdrvdef)
        call  pmf_ctrl_read_stritem(prm_fin,'fpdrvout',fpdrvout)
    end if
    if( mon_enabled ) then
        write(PMF_OUT,500)
        call  pmf_ctrl_read_stritem(prm_fin,'fmondef',fmondef)
        call  pmf_ctrl_read_stritem(prm_fin,'fmonout',fmonout)
    end if

    return

 10 format('# -------- General Setup')
100 format('# -------- Constrained Dynamics (CST)')
200 format('# -------- Restrained Dynamics (RST)')
300 format('# -------- Metadynamics (MTD)')
400 format('# -------- Adaptive Biasing Force Method (ABF)')
410 format('# -------- Adaptive Biasing Force Method (Testing - TABF)')
420 format('# -------- Umbrella Sampling / Adaptive Biasing Force Method (US-ABF)')
500 format('# -------- Monitoring (MON)')
700 format('# -------- String Method (STM)')
800 format('# -------- Path Driving (PDRV)')


end subroutine pmf_control_read_files

!===============================================================================
! Subroutine:  pmf_control_read_cvs
!===============================================================================

subroutine pmf_control_read_cvs(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{CVS}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fcvsdef(1:1) .eq. '{' ) then
        grpname = fcvsdef(2:len_trim(fcvsdef)-1)
        write(PMF_OUT,110) trim(grpname)
        ! open goup with name from cvsdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130) trim(grpname)
            pmf_enabled = .false.
            return
        end if
        call pmf_control_read_cvs_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fcvsdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fcvsdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to load file: ' // trim(fcvsdef) // '!')
        end if

        call pmf_control_read_cvs_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    call pmf_control_list_periodic_cvs()

    return

110 format('Collective variables are read from group: ',A)
120 format('Collective variables are read from file : ',A)
130 format(' >> No ',A,' group was specified - disabling PMFLib!')

end subroutine pmf_control_read_cvs

!===============================================================================
! Subroutine:  pmf_control_read_cvs_from_group
!===============================================================================

subroutine pmf_control_read_cvs_from_group(prm_fin)

    use pmf_dat
    use pmf_utils
    use prmfile
    use pmf_cvs
    use pmf_alloc_cv

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME)    :: cv_type
    integer                                :: i, alloc_failed
    logical                                :: erstult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfCVs = prmfile_count_group(prm_fin)

    if( NumOfCVs .le. 0 ) then
        ! on CV in current or specified group
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfCVs

    ! allocate constraint list -------------------------------------------------
    allocate(CVList(NumOfCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLib] Unable to allocate memory for collective variables!')
    end if

    ! enumerate sections -------------------------------------------------------
    erstult = prmfile_first_section(prm_fin)
    i = 1

    do while(erstult)
        erstult = prmfile_get_section_name(prm_fin,cv_type)
        write(PMF_OUT,*)
        write(PMF_OUT,130) i,trim(cv_type)

        call pmf_alloc_cv_allocate(cv_type,CVList(i)%cv)

        ! init and load CV data
        call CVList(i)%cv%reset_cv()
        CVList(i)%cv%idx = i
        call CVList(i)%cv%load_cv(prm_fin)

        erstult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    return

100 format('>>> INFO: No collective variables are defined. PMFLib is switched off!')
110 format('Number of collective variables : ',I4)
130 format('== Reading collective variable #',I4.4,' of type "',A,'"')

end subroutine pmf_control_read_cvs_from_group

!===============================================================================
! Subroutine:  pmf_control_read_paths
!===============================================================================

subroutine pmf_control_read_paths(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_GROUP_NAME)      :: grpname
    type(PRMFILE_TYPE)                     :: locprmfile
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{PATHS}',':')
    write(PMF_OUT,*)

    ! get name of group
    if( fpathsdef(1:1) .eq. '{' ) then
        grpname = fpathsdef(2:len_trim(fpathsdef)-1)
        write(PMF_OUT,110) trim(grpname)
        ! open goup with name from fpathsdef
        if( .not. prmfile_open_group(prm_fin,trim(grpname)) ) then
            write(PMF_OUT,130) trim(grpname)
            return
        end if
        call pmf_control_read_paths_from_group(prm_fin)
    else
        write(PMF_OUT,120) trim(fpathsdef)

        call prmfile_init(locprmfile)

        if( .not. prmfile_read(locprmfile,fpathsdef) ) then
            call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to load file: ' // trim(fpathsdef) // '!')
        end if

        call pmf_control_read_paths_from_group(locprmfile)

        call prmfile_clear(locprmfile)
    end if

    return

110 format('Paths are read from group: ',A)
120 format('Paths are read from file : ',A)
130 format(' >> No ',A,' group was specified.')

end subroutine pmf_control_read_paths

!===============================================================================
! Subroutine:  pmf_control_read_paths_from_group
!===============================================================================

subroutine pmf_control_read_paths_from_group(prm_fin)

    use pmf_dat
    use pmf_utils
    use prmfile
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME) :: section
    integer                             :: i, alloc_failed
    logical                             :: erstult
    ! --------------------------------------------------------------------------

    ! count number of sections in group
    NumOfPaths = prmfile_count_group(prm_fin)

    if( NumOfPaths .le. 0 ) then
        ! on PATH in current or specified group
        write(PMF_OUT,100)
        return
    end if

    write(PMF_OUT,110) NumOfPaths

    ! allocate constraint list -------------------------------------------------
    allocate(PathList(NumOfPaths), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLib] Unable to allocate memory for paths!')
    end if

    ! enumerate sections -------------------------------------------------------
    erstult = prmfile_first_section(prm_fin)
    i = 1
    do while(erstult)
        erstult = prmfile_get_section_name(prm_fin,section)
        if( trim(section) .ne. 'PATH' ) then
            call pmf_utils_exit(PMF_OUT,1,'[PMFLib] Unsupported section ''' // trim(section) // ''' in {PATHS} group!')
        end if
        write(PMF_OUT,*)
        write(PMF_OUT,130) i

        allocate(PathList(i)%path, stat = alloc_failed)
        if ( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[PMFLib] Unable to allocate memory for the path!')
        end if

        PathList(i)%path%idx = i

        call pmf_paths_load_path(prm_fin,PathList(i)%path)

        erstult = prmfile_next_section(prm_fin)
        i = i + 1
    end do

    return

100 format('>>> INFO: No paths are defined.')
110 format('Number of paths : ',I4)
130 format('== Reading path #',I4.4)

end subroutine pmf_control_read_paths_from_group

!===============================================================================
! Subroutine:  pmf_control_read_paths_gen_fake_cvs
!===============================================================================

subroutine pmf_control_read_paths_gen_fake_cvs(prm_fin)

    use pmf_dat
    use pmf_utils
    use prmfile
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! -----------------------------------------------
    character(PRMFILE_MAX_SECTION_NAME)    :: section
    integer                                :: ncvs, alloc_failed
    logical                                :: erstult
    ! --------------------------------------------------------------------------

    ! CVs among paths cannot overlap

    NumOfCVs = 0
    NumOfFakeCVs = 0

    ! get number of CVS  -------------------------------------------------------
    erstult = prmfile_first_section(prm_fin)
    do while(erstult)
        erstult = prmfile_get_section_name(prm_fin,section)
        if( trim(section) .ne. 'PATH' ) then
            cycle
        end if
        if( prmfile_get_integer_by_key(prm_fin,'ncvs',ncvs) ) then
            NumOfCVs = NumOfCVs + ncvs
        end if

        erstult = prmfile_next_section(prm_fin)
    end do

    write(PMF_OUT,110) NumOfCVs

    ! allocate constraint list -------------------------------------------------
    allocate(CVList(NumOfCVs), stat = alloc_failed)

    if ( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLib] Unable to allocate memory for collective variables!')
    end if

    ! now create individual CVs -------------------------------------------------
    erstult = prmfile_first_section(prm_fin)
    do while(erstult)
        erstult = prmfile_get_section_name(prm_fin,section)
        if( trim(section) .ne. 'PATH' ) then
            cycle
        end if
        call pmf_paths_load_fake_cvs(prm_fin)

        erstult = prmfile_next_section(prm_fin)
    end do

    return

110 format('Number of collective variables : ',I4)

end subroutine pmf_control_read_paths_gen_fake_cvs

!===============================================================================
! Subroutine:  pmf_control_list_periodic_cvs
!===============================================================================

subroutine pmf_control_list_periodic_cvs()

    use pmf_utils
    use pmf_constants
    use pmf_dat

    implicit none
    integer     :: i,pcvs
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Periodic collective variables','-')
    write(PMF_OUT,*)

    pcvs = 0
    do i=1,NumOfCVs
        if( CVList(i)%cv%is_periodic_cv() ) then
            pcvs = pcvs + 1
        end if
    end do

    write(PMF_OUT,10) pcvs

    if( pcvs .eq. 0 ) return

    write(PMF_OUT,*)
    write(PMF_OUT,20)
    write(PMF_OUT,30)

    do i=1,NumOfCVs
        if( CVList(i)%cv%is_periodic_cv() ) then
            write(PMF_OUT,40) i, CVList(i)%cv%ctype, trim(CVList(i)%cv%name), &
                              CVList(i)%cv%get_rvalue(CVList(i)%cv%get_min_cv_value()), &
                              CVList(i)%cv%get_rvalue(CVList(i)%cv%get_max_cv_value()), &
                              CVList(i)%cv%get_rvalue(CVList(i)%cv%get_period_cv_value()), &
                              '['//trim(CVList(i)%cv%get_ulabel())//']'
        end if
    end do

    return

10 format('Number of periodic CVs : ',I2)
20 format('ID   Type             Name            Min Value        Max Value       Periodicity   ')
30 format('-- ---------- -------------------- ---------------- ---------------- ----------------')

40 format(I2,1X,A10,1X,A20,1X,F16.6,1X,F16.6,1X,F16.6,1X,A)

end subroutine pmf_control_list_periodic_cvs

!===============================================================================
! Subroutine:  pmf_control_read_method_cvs_and_paths
!===============================================================================

subroutine pmf_control_read_method_cvs_and_paths(prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use pmf_init
    use mon_control
    use rst_control
    use abf_control
    use tabf_control
    use usabf_control
    use abp_control
    use mtd_control
    use cst_control
    use stm_control
    use pdrv_control

    implicit none
    type(PRMFILE_TYPE),intent(inout)       :: prm_fin
    ! --------------------------------------------------------------------------

    ! this must be read first as it influences other CV readers
    if( pdrv_enabled ) then
        call pdrv_control_read_pdrvs(prm_fin)
    end if

    ! read method CVs setup -------------------------
    if( abf_enabled ) then
        call abf_control_read_cvs(prm_fin)
    end if
    if( tabf_enabled ) then
        call tabf_control_read_cvs(prm_fin)
    end if
    if( usabf_enabled ) then
        call usabf_control_read_cvs(prm_fin)
    end if
    if( abp_enabled ) then
        call abp_control_read_cvs(prm_fin)
    end if
    if( mtd_enabled ) then
        call mtd_control_read_cvs(prm_fin)
    end if
    if( stm_enabled ) then
        call stm_control_read_cvs(prm_fin)
    end if
    if( cst_enabled ) then
        call cst_control_read_cvs(prm_fin)
    end if
    if( rst_enabled ) then
        call rst_control_read_cvs(prm_fin)
    end if
    if( mon_enabled ) then
        call mon_control_read_cvs(prm_fin)
    end if

end subroutine pmf_control_read_method_cvs_and_paths

!===============================================================================

end module pmf_control
