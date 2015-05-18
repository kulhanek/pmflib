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

program pmf_path_derivatives

    use pmf_sizes
    use pmf_constants
    use pmf_utils
    use pmf_control
    use prmfile
    use pmf_paths

    character(len=PMF_MAX_PATH) :: inname       ! input file name
    character(len=PMF_MAX_PATH) :: outname      ! output file name
    ! --------------------------------------------
    type(PRMFILE_TYPE)          :: locprmfile
    integer                     :: i
    logical                     :: res
    ! --------------------------------------------------------------------------

    call pmf_utils_header('Path derivatives')

    ! test number of input arguments
    if( command_argument_count() .ne. 2 ) then
        call print_usage()
        call pmf_utils_exit(PMF_OUT,1,'Incorrect number of arguments was specified (two expected)!')
    end if

    call get_command_argument(1, inname)
    call get_command_argument(2, outname)

    ! load paths ---------------------------------------------------------------
    write(PMF_OUT,*)
    write(PMF_OUT,100) trim(inname)

    call prmfile_init(locprmfile)

    if( .not. prmfile_read(locprmfile,inname) ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFLIB] Unable to load file: ' // trim(inname) // '!')
    end if

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Generating fake CVS',':')
    write(PMF_OUT,*)
    if( .not. prmfile_open_group(locprmfile,'PATHS') ) then
        res = prmfile_open_group(locprmfile,'MAIN')
    end if
    call pmf_control_read_paths_gen_fake_cvs(locprmfile)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{PATHS}',':')
    write(PMF_OUT,*)
    if( .not. prmfile_open_group(locprmfile,'PATHS') ) then
        res = prmfile_open_group(locprmfile,'MAIN')
    end if
    call pmf_control_read_paths_from_group(locprmfile)

    call prmfile_clear(locprmfile)

    ! save paths ---------------------------------------------------------------
    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'Saving derivatives',':')
    write(PMF_OUT,*)
    write(PMF_OUT,110) trim(outname)

    call pmf_utils_open(PDRV_OUT,outname,'R')

    do i=1,NumOfPaths
        call pmf_paths_write_path_derivatives(PDRV_OUT,PathList(i)%path)
    end do

    close(PDRV_OUT)

    call pmf_utils_footer('Path derivatives')

100 format('Paths are read from file : ',A)
110 format('Derivatives of paths are written to file : ',A)
contains

!===============================================================================
! subroutine:  print_usage
!===============================================================================

subroutine print_usage()

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(PMF_OUT,10)
    write(PMF_OUT,*)

    return

10 format('    pmf-path-derivatives <input> <output>')

end subroutine print_usage

!===============================================================================

end program pmf_path_derivatives
