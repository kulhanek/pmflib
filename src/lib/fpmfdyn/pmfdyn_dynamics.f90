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

module pmfdyn_dynamics

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine md_startup

    use smf_utils
    use smf_profiling
    use pmfdyn_system_control
    use pmfdyn_system
    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat
    use pmfdyn_thermostat
    use pmfdyn_restraints

    implicit none
    ! --------------------------------------------------------------------------

    call start_timer(INITIALIZATION_TIMER)

    ! initialize subsystems -------------------------
    call process_control       ! read control data
    call read_restart          ! and read cvs and velocities
    call open_output_files     ! open output files

    dt  = stepsize * PMF_DT2VDT
    idt = 1.0d0 / dt

    call pmf_pmfdyn_init(ControlPrmfile,natoms,nsteps,stepsize,Temp0,mass,xtop)    ! pmflib - part I

    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Initialising dynamics', ':')

    call init_thermostat_subsystem(md_x,md_v)
    call init_restraints

    write(PMF_OUT,*)

    call stop_timer(INITIALIZATION_TIMER)

    return

end subroutine md_startup

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmfdyn_init(ctrlprmfile,natoms,nsteps,stepsize,temp0,amass,ax)

    use prmfile
    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_utils
    use pmf_pbc
    use pmf_core_lf
    use pmf_mask
    use pmf_control

    implicit none
    type(PRMFILE_TYPE) :: ctrlprmfile
    integer            :: natoms
    integer            :: nsteps
    real(PMFDP)        :: stepsize
    real(PMFDP)        :: temp0
    real(PMFDP)        :: amass(:)
    real(PMFDP)        :: ax(:,:)
    ! -----------------------------------------------
    integer            :: i,has_box
    real(PMFDP)        :: cbox(3)
    ! ----------------------------------------------------------------------------

    ! setup conversion factors
    MassConv         = 1.0d0               ! g/mol -> g/mol
    LengthConv       = 1.0d0               ! A -> A
    AngleConv        = 1.0d0               ! rad -> rad
    TimeConv         = 1.0d0               ! fs -> fs
    VelocityConv     = 1.0d0               ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0               ! kcal/mol -> kcal/mol
    ForceConv        = -1.0d0              ! derivatives kcal/mol/A -> kcal/mol/A

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_LEAP_FROG,natoms,0,nsteps,stepsize,0.0d0,temp0)

    ! init mask subsystem
    call pmf_pbc_get_cbox(has_box,cbox)
    call pmf_mask_topo_init(natoms,1,has_box,cbox(1),cbox(2),cbox(3))

    ! set only one residue
    call pmf_mask_set_topo_residue(1,"ALL",1)

    ! init mask topology atom masses and positions
    do i=1,natoms
        call pmf_mask_set_topo_atom(i,"","")
        call pmf_mask_set_topo_atom_mcrd(i,amass(i),ax(1,i),ax(2,i),ax(3,i))
    end do

    ! finalize mask
    call pmf_mask_topo_finalize()

    ! print PMFLib header with system description
    call pmf_init_title('PMFDYN')

    ! load method setups and CVs definition
    call pmf_control_read_pmflib_group(ctrlprmfile)

    ! read method CV setup
    call pmf_control_read_method_cvs_and_paths(ctrlprmfile)

    ! write footer
    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'{END}',':')
    write(PMF_OUT,*)

    ! now we check if everything was understood from control file
    if( prmfile_count_ulines(ctrlprmfile,'PMFLIB') .gt. 0 ) then
        write(PMF_OUT,'(/,a,/)') '=== [unprocessed options in {PMFLIB} group] ======================================'
        call prmfile_dump_group(ctrlprmfile,PMF_OUT,'PMFLIB',.true.)
        call pmf_utils_exit(PMF_OUT,1,'Some items in the control file were not understood (see above)!')
    end if

    ! init all methods
    call pmf_init_all(amass,ax)

    ! init methods
    call pmf_init_pmf_methods

    return

end subroutine pmf_pmfdyn_init

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine md_shutdown

    use smf_profiling
    use pmfdyn_system
    use pmfdyn_system_dat
    use pmf_utils
    use pmf_finalize
    use pmfdyn_system
    use smf_utils

    implicit none
    ! --------------------------------------------------------------------------

    call start_timer(FINALIZATION_TIMER)

    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Finalizing dynamics', ':')

    call pmf_finalize_all(.false.)

    ! write final restart file
    call write_restart

    ! close output files
    call close_output_files

    ! deallocate memory
    call system_shutdown

    call stop_timer(FINALIZATION_TIMER)

    return

end subroutine md_shutdown

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine md_run

    use smf_utils
    use smf_profiling
    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat
    use pmfdyn_thermostat
    use pmfdyn_system
    use pmf_dat
    use pmf_core_lf

    implicit none
    integer            :: i
    ! -----------------------------------------------------------------------------

    ! only init subsystems
    if( OnlyInit ) return

    call start_timer(CORE_TIMER)

    call centered_heading(PMF_OUT,'Molecular dynamics', ':')
    call write_mdout_header

    do istep = 1, nsteps

        ! update time position
        time = time + stepsize

        !===============================================================================
        ! get potential energy and derivatives
        call pot_energy()

        !===============================================================================
        ! update velocities from accelerations,

        ! backup old position
        md_old_x(:,:) = md_x(:,:)

        ! backup old velocities
        md_old_v(:,:) = md_v(:,:)

        ! do leap-frog step for
        select case(ThermostatType)
            case(THERMOSTAT_NONE,THERMOSTAT_BERENDSEN,THERMOSTAT_CHEATHAM)
                do i=1,natoms
                    md_v(:,i) = md_v(:,i) - md_d(:,i)*winv(i)*dt
                    md_x(:,i) = md_x(:,i) + md_v(:,i)*dt
                end do
            case(THERMOSTAT_LANGEVIN)
                call langevin_thermostat(md_v,md_d)
                do i=1,natoms
                    md_x(:,i) = md_x(:,i) + md_v(:,i)*dt
                end do
        end select

        !===============================================================================
        ! shake and bluemoon constraints
        if( cst_enabled ) then
            call pmf_core_lf_shake(md_x)
            md_v(:,:) = (md_x(:,:) - md_old_x(:,:)) * idt
        end if

        ! remove COM tran and rot
        if( (nscm_freq .gt. 0) .and. (mod(istep,nscm_freq) .eq. 0) ) call stop_com(md_x,md_v)

        !===============================================================================
        ! update temperature controls and remove hot atoms
        call adjust_thermostat_controls(md_v)

        ! calculate temperature
        call get_temperature(md_old_v,md_v,Ekin,EkinH)


        Etot = Epot+Ekin+Erst+Epmf

        !===============================================================================
        ! thermostat
        select case(ThermostatType)
            case(THERMOSTAT_BERENDSEN)
                call berendsen_thermostat(md_v)
            case(THERMOSTAT_CHEATHAM)
                call cheatham_thermostat(md_old_v,md_v)
        end select

        !===============================================================================
        ! print [intermediate] results (master node only)
        call write_results
    end do

    call stop_timer(CORE_TIMER)

    return

end subroutine md_run

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmfdyn_dynamics




