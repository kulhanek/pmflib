!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2010,2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

!===============================================================================
! subroutine pmf_cpmd_init_control_name
!===============================================================================

subroutine pmf_cpmd_init_control_name(mdin)

 use pmf_cpmd_dat

 implicit none
 character(*)   :: mdin
 ! -----------------------------------------------------------------------------

 ControlFileName = mdin

 return

end subroutine pmf_cpmd_init_control_name

!===============================================================================
! subroutine pmf_cpmd_init
!===============================================================================

subroutine pmf_cpmd_init(NSP,NAX,NSX,NA,PMA0,TAU0,NOMORE,DELT_IONS,TEMPW,IBRAV,CELLDM)

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_utils
    use pmf_cpmd_dat
    use pmf_cpmd_control
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers
    use pmf_cpmd_common

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: PMA0(NSX)            ! masses
    real(PMFDP)    :: TAU0(3,NAX,NSX)      ! coordinates
    integer        :: NOMORE               ! number of MD steps
    real(PMFDP)    :: DELT_IONS            ! time step in atomic unit
    real(PMFDP)    :: TEMPW                ! temperature in K
    integer        :: IBRAV                ! symmetry
    real(PMFDP)    :: CELLDM(6)            ! box dimensions
    ! -----------------------------------------------
    integer        :: ntb,latom,is,alloc_failed
    real(PMFDP)    :: lbox(6)
    integer        :: has_box
    real(PMFDP)    :: cbox(3)
    integer        :: anres                        ! number of residues
    integer        :: i
    ! --------------------------------------------------------------------------

    ! setup conversion factors ----------------------
    MassConv         = 1.0d0                           ! g/mol -> g/mol
    LengthConv       = PMF_AU2A                        ! a.u. -> A
    AngleConv        = 1.0d0                           ! rad -> rad
    TimeConv         = PMF_AU2FS                       ! a.u. -> fs
    VelocityConv     = PMF_AU2A/PMF_AU2FS*PMF_VDT2DT ! a.u./a.u. -> pmflib velocity
    EnergyConv       = PMF_HARTREE2KCL                 ! hartree -> kcal/mol
    ForceConv        = PMF_HARTREE2KCL/PMF_AU2A        ! hartree/A -> kcal/mol/A

    ! setup box -------------------------------------
    lbox(1) = CELLDM(1)
    lbox(2) = lbox(1)*CELLDM(2)
    lbox(3) = lbox(1)*CELLDM(3)
    lbox(4:6) = acos(CELLDM(4:6))

    ntb = SYS_NT
    if( IBRAV .ne. 0 ) ntb = SYS_NTV

    ! init pmf core ---------------------------------
    ! calculate the number of atoms
    latom=0
    do is=1,NSP
        latom=latom+NA(IS)
    end do

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_VEL_VERLET,latom,ntb,NOMORE,DELT_IONS,0.0d0,TEMPW)
    call pmf_pbc_set_box(lbox(1),lbox(2),lbox(3),lbox(4),lbox(5),lbox(6))

    ! remap cpmd arrays to libpmf arrays ------------
    ! allocate local arrays
    allocate(  tmp_mass(fnatoms),   &
            tmp_x(3,fnatoms),    &
            tmp_y(3,fnatoms),    &
            tmp_v(3,fnatoms),    &
            tmp_f(3,fnatoms), stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory used in pmf_init')
    end if

    ! copy coordinates to tmp_x array
    call remap_coords_in(NSP,NAX,NSX,NA,TAU0,tmp_x)
    call remap_masses_in(NSP,NSX,NA,PMA0,tmp_mass)

    ! init mask subsystem
    anres = 1
    call pmf_pbc_get_cbox(has_box,cbox)
    call pmf_mask_topo_init(fnatoms,anres,has_box,cbox(1),cbox(2),cbox(3))

    ! init mask topology atom masses and positions
    do i=1,fnatoms
        call pmf_mask_set_topo_atom_mcrd(i,tmp_mass(i),tmp_x(1,i),tmp_x(2,i),tmp_x(3,i))
    end do

    ! finalize mask
    call pmf_mask_topo_finalize()

    ! print PMFLib header with system description
    call pmf_init_title('CPMD')

    ! load method setups
    call pmf_cpmd_process_control_begin()

    ! print main header
    call pmf_init_all(tmp_mass,tmp_x)

    ! finalize input file processing
    call pmf_cpmd_process_control_end()

    ! end of pmflib initialization
    call pmf_init_pmf_methods()

    ! rewrite the setup
    fcanexmdloop     = .true.       ! the client is able to terminate md loop

    return

end subroutine pmf_cpmd_init

!===============================================================================
! subroutine pmf_cpmd_shouldexit
!===============================================================================

subroutine pmf_cpmd_shouldexit(exitcode)

    use pmf_constants
    use pmf_sizes
    use pmf_dat

    implicit none
    integer        :: exitcode       ! MD loop exit code
    ! --------------------------------------------------------------------------

    exitcode = fexit_mdloop

    return

end subroutine pmf_cpmd_shouldexit

!===============================================================================
! subroutine pmf_cpmd_shouldexit_mpi
!===============================================================================
#ifdef MPI
subroutine pmf_cpmd_shouldexit_mpi(exitcode)

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_utils
    use mpi

    implicit none
    integer        :: exitcode       ! MD loop exit code
    ! --------------------------------------------
    integer        :: ierr
    ! --------------------------------------------------------------------------

    ! this can be call from all processes, fexit_mdloop must be redistributed
    ! to them

    ! broadcast MD exit status
    call mpi_bcast(fexit_mdloop, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast fexit_mdloop variable!')
    end if

    exitcode = fexit_mdloop

    return

end subroutine pmf_cpmd_shouldexit_mpi
#endif

!===============================================================================
! subroutine pmf_cpmd_shake
!===============================================================================

subroutine pmf_cpmd_shake(NSP,NAX,NSX,NA,TAUP,CVELP)

    use pmf_dat
    use pmf_sizes
    use pmf_cpmd_common
    use pmf_core_vv
    use pmf_cpmd_dat

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates r_u(t+dt)->R(t+dt)
    real(PMFDP)    :: CVELP(3,NAX,NSX)     ! velocities  v_u(t')->??
    ! --------------------------------------------------------------------------

    if( .not. con_enabled ) return

    ! remap cpmd arrays to pmflib arrays
    call remap_coords_in(NSP,NAX,NSX,NA,TAUP,tmp_x)
    call remap_coords_in(NSP,NAX,NSX,NA,CVELP,tmp_v)

    ! do shake
    call pmf_core_vv_shake_SFR(tmp_x,tmp_v)

    ! remap arrays back
    call remap_coords_back(NSP,NAX,NSX,NA,TAUP,tmp_x)
    call remap_coords_back(NSP,NAX,NSX,NA,CVELP,tmp_v)

    return

end subroutine pmf_cpmd_shake

!===============================================================================
! subroutine pmf_cpmd_force
!===============================================================================

subroutine pmf_cpmd_force(NSP,NAX,NSX,NA,TAUP,FION)

    use pmf_sizes
    use pmf_core_vv
    use pmf_cpmd_dat
    use pmf_cpmd_common

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: FION(3,NAX,NSX)      ! forces
    real(PMFDP)    :: ENE
    ! -----------------------------------------------
    real(PMFDP)    :: epmf
    ! --------------------------------------------------------------------------

    ! copy coordinates to tmp_x array
    call remap_coords_in(NSP,NAX,NSX,NA,TAUP,tmp_x)

    ! it is necessary to clean tmp_f because forces are added
    tmp_f(:,:) = 0.0d0
    call remap_force_in_add(NSP,NAX,NSX,NA,FION,tmp_f)

    epmf = 0.0d0

    call pmf_core_vv_force_SRF(tmp_x,tmp_f,ENE,epmf)

    ! convert forces back ----------------------------------
    call remap_force_back(NSP,NAX,NSX,NA,FION,tmp_f)

    return

end subroutine pmf_cpmd_force

!===============================================================================
! subroutine pmf_cpmd_rattle
!===============================================================================

subroutine pmf_cpmd_rattle(NSP,NAX,NSX,NA,CVELP)

    use pmf_dat
    use pmf_sizes
    use pmf_cpmd_common
    use pmf_core_vv
    use pmf_cpmd_dat

    implicit none
    integer        :: NSP                  ! number of types
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: NA(NSX)              ! number of atoms of given type
    real(PMFDP)    :: CVELP(3,NAX,NSX)     ! velocities
    ! --------------------------------------------------------------------------

    if( .not. con_enabled ) return

    ! remap cpmd arrays to pmflib arrays
    call remap_coords_in(NSP,NAX,NSX,NA,CVELP,tmp_v)

    ! do rattle
    call pmf_core_vv_rattle_SFR(tmp_v)

    ! remap arrays back
    call remap_coords_back(NSP,NAX,NSX,NA,CVELP,tmp_v)

    return

end subroutine pmf_cpmd_rattle

!===============================================================================
! subroutine pmf_cpmd_finalize
!===============================================================================

subroutine pmf_cpmd_finalize

    use pmf_constants
    use pmf_utils
    use pmf_finalize
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Library Finalization', '-')

    call pmf_timers_stop_timer(PMFLIB_TOTAL_TIMER)
    call pmf_finalize_all(.true.)

end subroutine pmf_cpmd_finalize

!===============================================================================




