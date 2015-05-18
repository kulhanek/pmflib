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

module pmfdyn_thermostat_control

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_thermostat

    use prmfile
    use pmf_utils
    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat

    implicit none
    character(80)      :: string
    ! --------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [thermostat] ==============================================================='

    NeedRestart = .true.

    ! open first section
    if( .not. prmfile_open_section(ControlPrmfile,'thermostat') ) then
        write(PMF_OUT,15) 'berendsen'
        ThermostatType = THERMOSTAT_BERENDSEN
        write(PMF_OUT,25) Temp0
        write(PMF_OUT,75) prmfile_onoff(MonitorHotAtoms)
        if( MonitorHotAtoms ) then
            write(PMF_OUT,85) HotAtomTreshold
            write(PMF_OUT,95) prmfile_onoff(RemoveHotAtoms)
        end if
        write(PMF_OUT,105) Iseed
        return
    end if

    if( prmfile_get_string_by_key(ControlPrmfile,'type', string)) then
        select case(trim(string))
            case('none')
                write(PMF_OUT,10) 'none'
                ThermostatType = THERMOSTAT_NONE
            case('berendsen')
                write(PMF_OUT,10) 'berendsen'
                ThermostatType = THERMOSTAT_BERENDSEN
            case('cheatham')
                write(PMF_OUT,10) 'cheatham'
                ThermostatType = THERMOSTAT_CHEATHAM
            case('langevin')
                write(PMF_OUT,10) 'langevin'
                ThermostatType = THERMOSTAT_LANGEVIN
            case default
                call pmf_utils_exit(PMF_OUT,1,'Unsupported thermostat type ' // trim(string) // ' !' )
        end select
    else
        write(PMF_OUT,15) 'berendsen'
        ThermostatType = THERMOSTAT_BERENDSEN
    end if

    ! check some possible conflicts ------------------------------
    select case(ThermostatType)
        case(THERMOSTAT_NONE)
            if( prmfile_get_real8_by_key(ControlPrmfile,'target_temperature', Temp1) ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     'Target temperature (target_temperature) cannot be specified if none thermostat is used!')
            end if
            if( prmfile_get_real8_by_key(ControlPrmfile,'temperature', Temp1) ) then
                call pmf_utils_exit(PMF_OUT,1, &
                     'Temperature (temperature) cannot be specified if none thermostat is used!')
            end if
    end select

    if(.not. prmfile_get_real8_by_key(ControlPrmfile,'initial_temperature', Tmaxw)) then
        NeedRestart = .true.
        Tmaxw = -1.0d0
    else
        NeedRestart = .false.
        write(PMF_OUT,60) Tmaxw
    end if

    Temp1 = -1.0
    select case(ThermostatType)
        case(THERMOSTAT_NONE)
        case default
            if( prmfile_get_real8_by_key(ControlPrmfile,'temperature', Temp0)) then
                write(PMF_OUT,20) Temp0
            else
                if( Tmaxw .gt. -1.0d0 ) then
                    Temp0 = Tmaxw
                end if
                write(PMF_OUT,25) Temp0
            end if

            if( prmfile_get_real8_by_key(ControlPrmfile,'target_temperature', Temp1)) then
                write(PMF_OUT,30) Temp1
            end if
    end select

    if( prmfile_get_logical_by_key(ControlPrmfile,'monitor_hot_atoms', MonitorHotAtoms) ) then
        write(PMF_OUT,70) prmfile_onoff(MonitorHotAtoms)
    else
        write(PMF_OUT,75) prmfile_onoff(MonitorHotAtoms)
    end if

    if( prmfile_get_real8_by_key(ControlPrmfile,'hot_atom_treshold', HotAtomTreshold) ) then
        write(PMF_OUT,80) HotAtomTreshold
    else
        write(PMF_OUT,85) HotAtomTreshold
    end if

    if( prmfile_get_logical_by_key(ControlPrmfile,'remove_hot_atoms', RemoveHotAtoms) ) then
        write(PMF_OUT,90) prmfile_onoff(RemoveHotAtoms)
    else
        write(PMF_OUT,95) prmfile_onoff(RemoveHotAtoms)
    end if

    if( (RemoveHotAtoms .eqv. .true.) .and. (MonitorHotAtoms .eqv. .false.) ) then
        call pmf_utils_exit(PMF_OUT,1,'Hot atoms can be removed only if they are monitored!')
    end if

    if( prmfile_get_integer_by_key(ControlPrmfile,'random_seed', Iseed) ) then
        write(PMF_OUT,100) Iseed
    else
        write(PMF_OUT,105) Iseed
    end if

    ! read setup of individual thermostat
    call read_thermostat_spec

    return

 10  format ('Thermostat (type)                         = ',a12)
 15  format ('Thermostat (type)                         = ',a12,'               (default)')
 20  format ('Temperature (temperature)                 = ',f14.1,' [K]')
 25  format ('Temperature (temperature)                 = ',f14.1,' [K]         (default)')
 30  format ('Target temperature (target_temperature)   = ',f14.1,' [K]')
 60  format ('Initial temperature (initial_temperature) = ',f14.1)
 70  format ('Monitor hot atoms (monitor_hot_atoms)     = ',a12)
 75  format ('Monitor hot atoms (monitor_hot_atoms)     = ',a12,'               (default)')
 80  format ('Hot atom treshold (hot_atom_treshold)     = ',f14.1,' [K]')
 85  format ('Hot atom treshold (hot_atom_treshold)     = ',f14.1,' [K]         (default)')
 90  format ('Remove hot atoms (remove_hot_atoms)       = ',a12)
 95  format ('Remove hot atoms (remove_hot_atoms)       = ',a12,'               (default)')
100  format ('Random seed (random_seed)                 = ',i12)
105  format ('Random seed (random_seed)                 = ',i12,'               (default)')

end subroutine read_thermostat

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_thermostat_spec

    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat

    implicit none
    ! --------------------------------------------------------------------------

    select case(ThermostatType)
        case(THERMOSTAT_NONE)
            ! nothing to read
        case(THERMOSTAT_BERENDSEN)
            call read_thermostat_berendsen
        case(THERMOSTAT_CHEATHAM)
            call read_thermostat_cheatham
        case(THERMOSTAT_LANGEVIN)
            call read_thermostat_langevin
    end select

end subroutine read_thermostat_spec

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_thermostat_berendsen

    use prmfile
    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [berendsen_thermostat] ====================================================='

    ! open first section
    if( .not. prmfile_open_section(ControlPrmfile,'berendsen_thermostat') ) then
        write(PMF_OUT,45) tau_T
        return
    end if

    if( prmfile_get_real8_by_key(ControlPrmfile,'bath_coupling', Tau_T)) then
        write(PMF_OUT,40) tau_T
    else
        write(PMF_OUT,45) tau_T
    end if

    return

 40  format ('Bath coupling (bath_coupling)             = ',f16.3,' [fs]')
 45  format ('Bath coupling (bath_coupling)             = ',f16.3,' [fs]        (default)')

end subroutine read_thermostat_berendsen

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_thermostat_cheatham

     use prmfile
     use pmfdyn_system_dat
     use pmfdyn_thermostat_dat

     implicit none
     ! -------------------------------------------------------------------------

     write(PMF_OUT,'(/,a)') '=== [cheatham_thermostat] ======================================================'

     ! open first section
     if( .not. prmfile_open_section(ControlPrmfile,'cheatham_thermostat') ) then
        write(PMF_OUT,45) tau_T
        return
     end if

     if( prmfile_get_real8_by_key(ControlPrmfile,'bath_coupling', Tau_T)) then
        write(PMF_OUT,40) tau_T
     else
        write(PMF_OUT,45) tau_T
     end if

     return

 40  format ('Bath coupling (bath_coupling)             = ',f16.3,' [fs]')
 45  format ('Bath coupling (bath_coupling)             = ',f16.3,' [fs]        (default)')

end subroutine read_thermostat_cheatham

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_thermostat_langevin

    use prmfile
    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a)') '=== [langevin_thermostat] ======================================================'

    ! open first section
    if( .not. prmfile_open_section(ControlPrmfile,'langevin_thermostat') ) then
        write(PMF_OUT,55) Gamma_T
        return
    end if

    if( prmfile_get_real8_by_key(ControlPrmfile,'friction_coefficient', Gamma_T)) then
        write(PMF_OUT,50) Gamma_T
    else
        write(PMF_OUT,55) Gamma_T
    end if

    return

50  format ('Friction coefficient (friction_coefficient) = ',f16.3,' [fs^-1]')
55  format ('Friction coefficient (friction_coefficient) = ',f16.3,' [fs^-1]    (default)')

end subroutine read_thermostat_langevin

!===============================================================================

end module pmfdyn_thermostat_control
