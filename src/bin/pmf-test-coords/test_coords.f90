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

program test_coords

    use pmf_sizes
    use pmf_utils
    use test_coords_dat
    use test_coords_control
    use test_coords_utils
    use pmf_timers

    implicit none
    integer             ::  command_argument_count
    ! --------------------------------------------------------------------------

    call pmf_utils_header('Test CVs')

    ! checkk number of arguments ---------------------
    if( command_argument_count() .eq. 0 ) then
        call print_usage
        write(PMF_OUT,*)
        stop
    end if

    if( command_argument_count() .ne. 2 ) then
        call print_usage
        call pmf_utils_exit(PMF_OUT,1,'Incorrect number of arguments was specified (two expected)!')
    end if

    ! get arguments and print -----------------------
    call get_command_argument(1, ControlFile)
    call get_command_argument(2, CoordFile)

    ! init timers
    call pmf_timers_init_top()

    ! read control ----------------------------------
    call process_control

    ! initialization --------------------------------
    call initialization

    ! load coordinates and test them ----------------
    call load_coordinates

    ! test coordinates
    call test_coordinates

    ! print final results ---------------------------
    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'FINAL RESULTS', ':')
    write(PMF_OUT,*)
    write(PMF_OUT,10) total_num_of_tests,num_of_failed_tests

    ! print header ----------------------------------
    call pmf_utils_footer('Test CVs')

    stop

10 format('Total number of tests: ',I6,'                     Number of failed tests: ',I6)

!===============================================================================

end program test_coords


