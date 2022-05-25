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

module pmfdyn_system

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine init_timers

    use smf_profiling_dat
    use smf_profiling
    use pmfdyn_system_dat
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    ! add standard timers --------------------------------
    INITIALIZATION_TIMER   = add_timer(TOTAL_TIMER,'Initialization')
    CORE_TIMER             = add_timer(TOTAL_TIMER,'Program core')
        FORCES_TIMER           = add_timer(CORE_TIMER,'Forces')
            EXTERNAL_TIMER          = add_timer(FORCES_TIMER,'External forces')
            RE_TIMER                = add_timer(FORCES_TIMER,'Restraints')
        PMFLIB_TOTAL_TIMER     = add_timer(CORE_TIMER,'PMFLib')
    FINALIZATION_TIMER     = add_timer(TOTAL_TIMER,'Finalization')

end subroutine init_timers

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_restart

    use smf_xyzfile_type
    use smf_xyzfile
    use smf_periodic_table
    use pmfdyn_system_dat
    use pmf_utils

    implicit none
    type(XYZFILE_TYPE)     :: restart_xyz
    integer                :: alloc_failed, i
    ! --------------------------------------------------------------------------

    ! init xyz file structure
    call init_xyz(restart_xyz)

    ! open xyz file
    call open_xyz(IO_RST,InputRstFile,restart_xyz,'OLD')

    ! read coordinates
    call read_xyz(IO_RST,restart_xyz)

    natoms = restart_xyz%natoms

    if( natoms .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'System does not contain any atom!')
    end if

    ! allocate all arrays
    allocate(  xtop(3,natoms), &
            md_x(3,natoms), &
            md_old_x(3,natoms), &
            md_d(3,natoms), &
            md_v(3,natoms), &
            md_old_v(3,natoms), &
            symbols(natoms), &
            mass(natoms), &
            winv(natoms), &
            heavy(natoms), &
            stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for arrays used in pmf-dyn!')
    end if

    xtop(:,:)      = restart_xyz%cvs(:,:)
    md_x(:,:)      = xtop(:,:)
    md_old_x(:,:)  = 0.0d0
    md_d(:,:)      = 0.0d0
    md_v(:,:)      = 0.0d0
    md_old_v(:,:)  = 0.0d0
    mass(:)        = 0.0d0
    winv(:)        = 0.0d0

    ! init masses
    do i=1,natoms
        mass(i) = SearchMassBySymbol(restart_xyz%symbols(i))
        if( mass(i) .eq. 0.0 ) then
            call pmf_utils_exit(PMF_OUT,1,'Some atom has zero mass!')
        end if
        if( mass(i) .gt. 1.5d0 ) then
            heavy(i) = .true.
        else
            heavy(i) = .false.
        end if
        winv(i) = 1.0d0/mass(i)
        symbols(i) = restart_xyz%symbols(i)
    end do

    if( NeedRestart .eqv. .false. ) then
        ! close restart file
        call close_xyz(IO_RST,restart_xyz)

        return
    end if

    ! try to read velocities
    if( is_next_xyz_record(IO_RST,restart_xyz) .eqv. .false. ) then
        call pmf_utils_exit(PMF_OUT,1,'Velocities are required but they were not found!')
    end if

    ! read velocities
    call read_xyz(IO_RST,restart_xyz)

    md_v(:,:)      = restart_xyz%cvs(:,:)

    Restart = .true.

    ! close restart file
    call close_xyz(IO_RST,restart_xyz)

    return

end subroutine read_restart

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_restart

    use smf_xyzfile_type
    use smf_xyzfile
    use smf_periodic_table
    use pmf_utils
    use pmfdyn_system_dat

    implicit none
    type(XYZFILE_TYPE)     :: restart_xyz
    ! --------------------------------------------------------------------------

    ! allocate arrays
    call allocate_xyz(restart_xyz,natoms)

    ! open xyz file
    call open_xyz(IO_RST,OutputRstFile,restart_xyz,'UNKNOWN')

    ! copy symbols and coordinates
    restart_xyz%comment = 'cvs'
    restart_xyz%symbols(:) = symbols(:)
    restart_xyz%cvs(:,:) = md_x(:,:)

    ! write coordinates
    call write_xyz(IO_RST,restart_xyz)

    restart_xyz%comment = 'velocities'
    restart_xyz%symbols(:) = symbols(:)
    restart_xyz%cvs(:,:) = md_v(:,:)

    ! write velocities
    call write_xyz(IO_RST,restart_xyz)

    ! close restart file
    call close_xyz(IO_RST,restart_xyz)

    return

end subroutine write_restart

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine open_output_files

    use smf_xyzfile
    use pmfdyn_system_dat

    implicit none
    ! -----------------------------------------------------------------------------

    ! trajectory
    if( traj_freq .gt. 0 ) then
        ! init xyz file structure
        call allocate_xyz(TrajectoryCrd,natoms)

        ! open xyz file
        call open_xyz(IO_TRJ,OutputTrjFile,TrajectoryCrd,'UNKNOWN')
    end if

    return

end subroutine open_output_files

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine close_output_files

    use smf_xyzfile
    use pmfdyn_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! trajectory
    if( traj_freq .gt. 0 ) then
        ! close file
        call close_xyz(IO_TRJ,TrajectoryCrd)
    end if

    return

end subroutine close_output_files

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_traj_snapshot

    use smf_xyzfile_type
    use smf_xyzfile
    use pmfdyn_system_dat
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    ! copy symbols and coordinates
    TrajectoryCrd%comment = 'trajectory'
    TrajectoryCrd%symbols(:) = symbols(:)
    TrajectoryCrd%cvs(:,:) = md_x(:,:)

    ! write coordinates
    call write_xyz(IO_TRJ,TrajectoryCrd)

    return

end subroutine write_traj_snapshot

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_results

    use pmfdyn_system_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! print [intermediate] results (master node only)

    if( traj_freq .gt. 0 .and. mod(istep,traj_freq) == 0 ) then
        call write_traj_snapshot()    ! trajectory
    end if

    if(  (out_freq .gt. 0 .and. (mod(istep,out_freq) == 0)) .or. istep .eq. nsteps .or. istep .eq. 1 ) then
        call write_mdout() ! results to stdout
        ! call md_timing(istep,nsteps)
    end if

    if( rst_freq .gt. 0 .and. mod(istep,rst_freq) .eq. 0 ) then
        call write_restart           ! restart
    end if

    return

end subroutine write_results

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_mdout_header

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,10)
    write(PMF_OUT,20)

    return

10 format('   Step      E(tot)     E(kin)     E(pot)     E(rst)     E(pmf)     Temp    ')
20 format('---------- ---------- ---------- ---------- ---------- ---------- ----------')

end subroutine write_mdout_header

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine write_mdout

    use pmfdyn_system_dat
    use pmfdyn_thermostat_dat

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,10) istep,Etot,Ekin,Epot,Erst,Epmf,TempA

    return

    10 format(I10,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2)

end subroutine write_mdout

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine system_shutdown

 return

end subroutine system_shutdown

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pot_energy

    use smf_profiling
    use pmfdyn_system_dat
    use pmfdyn_restraints
    use pmf_core_lf
    use pmfdyn_thermostat_dat

    implicit none
    ! --------------------------------------------------------------------------

    call start_timer(FORCES_TIMER)

    Epot = 0.0d0
    Erst = 0.0d0
    Epmf = 0.0d0
    md_d(:,:) = 0.0d0

    call pmf_core_lf_update_step

    ! calculate external forces ---------------------
    if( associated(pot_ext_energy_pts) ) then
        call pot_ext_energy_pts(md_x,md_d,Epot)
    end if

    ! calculate restraints --------------------------
    call potene_restraints(md_x,md_d,Erst)

    ! pmf force  -------------------------------
    !md_d(:,:) = - md_d(:,:)
    ! FIXME - EkinH,0.0
    call pmf_core_lf_force(md_x,md_v,md_d,Epot+Erst,Ekin,EkinH,0.0d0,Epmf)
    !md_d(:,:) = - md_d(:,:)

    call stop_timer(FORCES_TIMER)

    return

end subroutine pot_energy

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmfdyn_system
