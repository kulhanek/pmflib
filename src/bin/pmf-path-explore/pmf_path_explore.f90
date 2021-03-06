! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

program PMFPathExplore

    use explore
    use pmf_constants
    use pmf_system
    use pmf_utils
    use smf_profiling

    implicit none
    character(90),parameter        :: PROGNAME= 'pmf-path-explore: path explorer by string method'
    ! --------------------------------------------------------------------------

    call pmf_utils_header(PROGNAME)

    call init_profiling(PMF_OUT,50)
    call init_timers

    call start_timer(TOTAL_TIMER)

    call exp_startup

    call exp_run

    call exp_shutdown

    call stop_timer(TOTAL_TIMER)

    call write_timing
    call finalize_profiling

    call pmf_utils_footer(PROGNAME)

    stop

end program PMFPathExplore

!===============================================================================

