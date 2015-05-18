! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
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

module pmfdyn_system_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine process_control

    use pmf_utils
    use smf_utils
    use prmfile
    use pmfdyn_system_dat
    use pmfdyn_thermostat_control
    use pmfdyn_restraints_control

    implicit none
    integer    ::  nrst_args
    integer    ::  command_argument_count
    ! --------------------------------------------------------------------------

    ! get name of control file from the command line
    nrst_args = command_argument_count()

    ! check if control file was provided
    if( nrst_args .eq. 0 ) then
        call print_usage
        call pmf_utils_exit(PMF_OUT,1,'No input file specified on the command line!')
    end if
    if( nrst_args .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'More then one argument is provided!')
    end if

    ! process options
    call get_command_argument(nrst_args, ControlFile)   ! argument is name of control file

    ! write header ----------------------------------------------------------
    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Reading control file', '-')
    write(PMF_OUT,'(a,a)') 'Control file name : ',trim(ControlFile)

    call prmfile_init(ControlPrmfile)

    if( .not. prmfile_read(ControlPrmfile,ControlFile) ) then
        call pmf_utils_exit(PMF_OUT,1,'specified control file cannot be opened!')
    end if

    if( .not. prmfile_open_group(ControlPrmfile,'MAIN') ) then
        call pmf_utils_exit(PMF_OUT,1,'Specified control file does not contain {MAIN} group!')
    end if

    ! read user specificaton ------------------------------------------------
    call read_dynamics
    call read_thermostat
    call read_intervals
    call read_files
    call read_restraints

    ! now we check if everything was understood from control file
    if( prmfile_count_ulines(ControlPrmfile,'MAIN') .gt. 0 ) then
        write(PMF_OUT,'(/,a,/)') '=== [unprocessed options in {MAIN} group] ======================================'
        call prmfile_dump_group(ControlPrmfile,PMF_OUT,'MAIN',.true.)
        call pmf_utils_exit(PMF_OUT,1,'Some items in the control file were not understood (see above)!')
    end if

    return

end subroutine process_control

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine print_usage

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a,/)') '=== [usage] ===================================================================='

    write(PMF_OUT,10)
    return

    10 format('    pmf-dyn <control>')

end subroutine print_usage

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_dynamics

    use prmfile
    use pmf_utils
    use pmfdyn_system_dat


    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [dynamics] ================================================================='

    ! open first section
    if( .not. prmfile_open_section(ControlPrmfile,'dynamics') ) then
        call pmf_utils_exit(PMF_OUT,1,'[dynamics] section has to be specified in control file!')
    end if

    ! timings
    if(.not. prmfile_get_integer_by_key(ControlPrmfile,'steps', nsteps)) then
        call pmf_utils_exit(PMF_OUT,1,'Steps are not specified!')
    end if

    if(.not. prmfile_get_real8_by_key(ControlPrmfile,'stepsize', stepsize)) then
        call pmf_utils_exit(PMF_OUT,1,'Stepsize is not specified!')
    end if
    write(PMF_OUT,10) nsteps
    write(PMF_OUT,20) stepsize

    if( prmfile_get_logical_by_key(ControlPrmfile,'onlyinit', OnlyInit)) then
        write(PMF_OUT,30) prmfile_onoff(OnlyInit)
    else
        write(PMF_OUT,35) prmfile_onoff(OnlyInit)
    end if

    return

 10  format ('Number of MD steps (steps)                = ',i12)
 20  format ('Stepsize (stepsize)                       = ',f16.3,' [fs]')
 30  format ('Only init engine (onlyinit)               = ',a12)
 35  format ('Only init engine (onlyinit)               = ',a12,'               (default)')

end subroutine read_dynamics

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_files

    use pmf_utils
    use pmf_dat
    use pmfdyn_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [files] ===================================================================='

    if(.not. prmfile_open_section(ControlPrmfile,'files')) then
        call pmf_utils_exit(PMF_OUT,1,'[files] section not found.')
    end if

    if(.not. prmfile_get_string_by_key(ControlPrmfile,'initial', InputRstFile)) then
        call pmf_utils_exit(PMF_OUT,1,'Initial coordinates (initial) are not specified!')
    end if
    write (PMF_OUT,10) trim(InputRstFile)

    ! final coordinates
    if(.not. prmfile_get_string_by_key(ControlPrmfile,'final', OutputRstFile)) then
        call pmf_utils_exit(PMF_OUT,1,'Final coordinate file (file) not specified!')
    end if
    write (PMF_OUT,30) trim(OutputRstFile)

    ! trajectory file
    if(.not. prmfile_get_string_by_key(ControlPrmfile,'trajectory', OutputTrjFile)) then
        if(traj_freq > 0) then
            call pmf_utils_exit(PMF_OUT,1,'Trajectory file name (trajectory) required to write trajectory!')
        end if
    else
        write (PMF_OUT,50) trim(OutputTrjFile)
        if(traj_freq < 0) then
            call pmf_utils_exit(PMF_OUT,1,'>>> Error: Trajectory file is given but no output interval!')
        end if
    end if

    return

 10  format ('Input restart file (initial)              = ',a)
 30  format ('Output restart file (final)               = ',a)
 50  format ('Trajectory file (trajectory)              = ',a)

end subroutine read_files

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_intervals

    use prmfile
    use pmf_utils
    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat

    implicit none
    ! -----------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [intervals] ================================================================'

    if(.not. prmfile_open_section(ControlPrmfile,'intervals')) then
        write(PMF_OUT,21) out_freq
        write(PMF_OUT,41) traj_freq
        write(PMF_OUT,61) nscm_freq
        return
    end if

    if(.not. prmfile_get_integer_by_key(ControlPrmfile,'output', out_freq)) then
        write(PMF_OUT,21) out_freq
    else
        write(PMF_OUT,20) out_freq
    end if

    if(.not. prmfile_get_integer_by_key(ControlPrmfile,'trajectory', traj_freq)) then
        write(PMF_OUT,41) traj_freq
    else
        write(PMF_OUT,40) traj_freq
    end if

    if(.not. prmfile_get_integer_by_key(ControlPrmfile,'nscm', nscm_freq)) then
        write(PMF_OUT,61) nscm_freq
    else
        write(PMF_OUT,60) nscm_freq
    end if

    ! avoid division by zero if interval is zero

    if(out_freq .le. 0) then
        out_freq = -1
    end if

    if(traj_freq .le. 0) then
        traj_freq = -1
    end if

    if(nscm_freq .le. 0) then
        nscm_freq = -1
    end if

    if( ThermostatType .eq. THERMOSTAT_LANGEVIN .and. nscm_freq .gt. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'nscm has to be zero if Langevin thermostat is used!')
    end if

    return

 20  format ('Dynamics summary interval (output)        = ',i12)
 21  format ('Dynamics summary interval (output)        = ',i12,'               (default)')

 40  format ('Trajectory write interval (trajectory)    = ',i12)
 41  format ('Trajectory write interval (trajectory)    = ',i12,'               (default)')

 60  format ('COM moment removal interval (nscm)        = ',i12)
 61  format ('COM moment removal interval (nscm)        = ',i12,'               (default)')

end subroutine read_intervals

! ------------------------------------------------------------------------------

end module pmfdyn_system_control
