! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

module test_coords_control

use pmf_sizes
use pmf_constants
use pmf_core

implicit none
contains

!===============================================================================
! Subroutine:  process_control
!===============================================================================

subroutine process_control

 use prmfile
 use pmf_utils
 use pmf_constants
 use pmf_init
 use pmf_pbc
 use pmf_mask
 use test_coords_dat

 implicit none
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'Reading control file', '-')
 write(PMF_OUT,'(a,a)') 'Control file name : ',trim(ControlFile)

 call prmfile_init(ControlPrmfile)

 if( .not. prmfile_read(ControlPrmfile,ControlFile) ) then
    call pmf_utils_exit(PMF_OUT,1,'specified control file cannot be opened!')
 end if

 if( .not. prmfile_open_group(ControlPrmfile,'MAIN') ) then
    call pmf_utils_exit(PMF_OUT,1,'Specified control file does not contain {MAIN} group!')
 end if

 ! read user specificaton ------------------------
 call read_main
 call read_files

 ! now we check if everything was understood from control file
 if( prmfile_count_ulines(ControlPrmfile,'MAIN') .gt. 0 ) then
    write(PMF_OUT,'(/,a,/)') '=== [unprocessed options in {MAIN} group] ======================================' 
    call prmfile_dump_group(ControlPrmfile,PMF_OUT,'MAIN',.true.)
    call pmf_utils_exit(PMF_OUT,1,'Some items in the control file were not understood (see above)!')
 end if

end subroutine process_control

!===============================================================================
! Subroutine:  read_main
!===============================================================================

subroutine read_main

 use prmfile
 use pmf_utils
 use pmf_dat
 use pmf_core
 use test_coords_dat

 implicit none
 character(80)      :: string
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,'(/,a,/)') '=== [test-cvs] ==============================================================' 

 ! try open section
 if( .not. prmfile_open_section(ControlPrmfile,'test-cvs') ) then
    call pmf_utils_exit(PMF_OUT,1,'Section [test-cvs] is mandatory!')
 end if

 if( prmfile_get_string_by_key(ControlPrmfile,'test_type', string)) then
    select case(trim(string))
        case('deriv')
            test_type = TEST_DERIV
            write(PMF_OUT,1) 'deriv'
        case('identity')
            test_type = TEST_IDENTITY
            write(PMF_OUT,1) 'identity'
        case('value')
            test_type = TEST_VALUE
            write(PMF_OUT,1) 'value'
        case('combined')
            test_type = TEST_COMBINED
            write(PMF_OUT,1) 'combined'
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported test type!')
    end select
 else
    call pmf_utils_exit(PMF_OUT,1,'Missing test_type!')
 end if

 max_atoms = 10

 if(prmfile_get_integer_by_key(ControlPrmfile,'max_atoms', max_atoms)) then
    write(PMF_OUT,10) max_atoms
 else
    write(PMF_OUT,15) max_atoms
 end if
 if(prmfile_get_real8_by_key(ControlPrmfile,'max_radius', max_radius)) then
    write(PMF_OUT,20) max_radius
 else
    write(PMF_OUT,25) max_radius
 end if
 if(prmfile_get_real8_by_key(ControlPrmfile,'max_mass', max_mass)) then
    write(PMF_OUT,30) max_mass
 else
    write(PMF_OUT,35) max_mass
 end if
 if(prmfile_get_real8_by_key(ControlPrmfile,'uniform_mass', uniform_mass)) then
    write(PMF_OUT,36) uniform_mass
 else
    write(PMF_OUT,37) uniform_mass
 end if
 if(prmfile_get_integer_by_key(ControlPrmfile,'num_of_tests', num_of_tests)) then
    write(PMF_OUT,40) num_of_tests
 else
    write(PMF_OUT,45) num_of_tests
 end if
 if(prmfile_get_real8_by_key(ControlPrmfile,'num_diff', num_diff)) then
    write(PMF_OUT,50) num_diff
 else
    write(PMF_OUT,55) num_diff
 end if
 if(prmfile_get_real8_by_key(ControlPrmfile,'alarm_treshold', alarm_treshold)) then
    write(PMF_OUT,60) alarm_treshold
 else
    write(PMF_OUT,65) alarm_treshold
 end if
 if(prmfile_get_logical_by_key(ControlPrmfile,'stop_when_error', stop_when_error)) then
    write(PMF_OUT,70) prmfile_onoff(stop_when_error)
 else
    write(PMF_OUT,75) prmfile_onoff(stop_when_error)
 end if
 if(prmfile_get_logical_by_key(ControlPrmfile,'verbose_print', verbose_print)) then
    write(PMF_OUT,80) prmfile_onoff(verbose_print)
 else
    write(PMF_OUT,85) prmfile_onoff(verbose_print)
 end if

 return

  1 format ('test_type                              = ',a12)
 10 format ('max_atoms                              = ',i12)
 15 format ('max_atoms                              = ',i12,'                  (default)')
 20 format ('max_radius                             = ',F12.4)
 25 format ('max_radius                             = ',F12.4,'                  (default)')
 30 format ('max_mass                               = ',F12.4)
 35 format ('max_mass                               = ',F12.4,'                  (default)')
 36 format ('uniform_mass                           = ',F12.4)
 37 format ('uniform_mass                           = ',F12.4,'                  (default)')
 40 format ('num_of_tests                           = ',i12)
 45 format ('num_of_tests                           = ',i12,'                  (default)')
 50 format ('num_diff                               = ',E12.4)
 55 format ('num_diff                               = ',E12.4,'                  (default)')
 60 format ('alarm_treshold                         = ',E12.4)
 65 format ('alarm_treshold                         = ',E12.4,'                  (default)')
 70 format ('stop_when_error                        = ',a12)
 75 format ('stop_when_error                        = ',a12,'                  (default)')
 80 format ('verbose_print                          = ',a12)
 85 format ('verbose_print                          = ',a12,'                  (default)')

end subroutine read_main

!===============================================================================
! Subroutine:  read_files
!===============================================================================

subroutine read_files

 use prmfile
 use pmf_utils
 use pmf_dat
 use pmf_core
 use test_coords_dat

 implicit none
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,'(/,a,/)') '=== [files] ===================================================================='

 ! try open section
 if( .not. prmfile_open_section(ControlPrmfile,'files') ) then
    write(PMF_OUT,30)
    return
 end if

 if(prmfile_get_string_by_key(ControlPrmfile,'snapshots', SnapshotFile)) then
    write(PMF_OUT,20) trim(SnapshotFile)
    UseExternalSnapshosts = .true.
 else
    write(PMF_OUT,30)
 end if

 return

 20 format ('snapshots                              = ',a)
 30 format ('snapshots                              = ','-not used-')

end subroutine read_files

!===============================================================================
! Subroutine:  load_coordinates
!===============================================================================

subroutine load_coordinates

 use prmfile
 use pmf_utils
 use pmf_dat
 use pmf_core
 use pmf_control
 use test_coords_dat

 implicit none
 character(PRMFILE_MAX_GROUP_NAME)      :: grpname
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,*)
 call pmf_utils_heading(PMF_OUT,'TESTED COORDINATES', ':')

 ! get name of group
 if( CoordFile(1:1) .eq. '{' ) then
    grpname = CoordFile(2:len_trim(CoordFile)-1)
     write(PMF_OUT,110) grpname
    ! open goup with name from abfdef
    if( .not. prmfile_open_group(ControlPrmfile,trim(grpname)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to open group {' // trim(grpname) // '}!')
    end if
    call pmf_control_read_cvs_from_group(ControlPrmfile)
 else
    write(PMF_OUT,120) trim(CoordFile)

    call prmfile_init(CoordPrmfile)

    if( .not. prmfile_read(CoordPrmfile,CoordFile) ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to load file: ' // trim(CoordFile) // '!')
    end if

    call pmf_control_read_cvs_from_group(CoordPrmfile)

    call prmfile_clear(CoordPrmfile)
 end if

 return

110 format('Coordinates are read from group: ',A)
120 format('Coordinates are read from file : ',A)

end subroutine load_coordinates

!===============================================================================

end module test_coords_control
