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

module pmf_cp2k

! CALLING SEQUENCE
!    -> pmf_cp2k_init_control_name
!    -> pmf_cp2k_init
!    -> pmf_cp2k_rattle
!    -> pmf_cp2k_force
!    LOOP
!       -> pmf_cp2k_shake
!       -> pmf_cp2k_force
!       -> pmf_cp2k_rattle
!    END LOOP
!    -> pmf_cp2k_finalize

contains

!===============================================================================
! subroutine pmf_cp2k_init_control_name
!===============================================================================

subroutine pmf_cp2k_init_control_name(mdin)

    use pmf_cp2k_dat

    implicit none
    character(*)   :: mdin
    ! --------------------------------------------------------------------------

    ControlFileName = mdin

    return

end subroutine pmf_cp2k_init_control_name

!===============================================================================
! subroutine pmf_cp2k_init
!===============================================================================

subroutine pmf_cp2k_init(latoms,lnst,ldt,ltemp,lbox,lx,lm)

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_utils
    use pmf_cp2k_dat
    use pmf_cp2k_control
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers

    implicit none
    integer        :: latoms               ! number of atoms
    integer        :: lnst                 ! number of steps
    real(PMFDP)    :: ldt                  ! time step
    real(PMFDP)    :: ltemp                ! temperature
    real(PMFDP)    :: lbox(3,3)            ! cell box
    real(PMFDP)    :: lx(3,latoms)         ! coordinates
    real(PMFDP)    :: lm(latoms)           ! masses
    ! -----------------------------------------------
    integer        :: ntb
    integer        :: has_box
    real(PMFDP)    :: cbox_a,cbox_b,cbox_c
    integer        :: anres                        ! number of residues
    integer        :: i
    ! --------------------------------------------------------------------------

    ! setup conversion factors ----------------------
    MassConv         = PMF_AU2AMU                      ! a.u. -> g/mol (amu)
    LengthConv       = PMF_AU2A                        ! a.u. -> A
    AngleConv        = 1.0d0                           ! rad -> rad
    TimeConv         = PMF_AU2FS                       ! a.u. -> fs
    VelocityConv     = PMF_AU2A/PMF_AU2FS*PMF_VDT2DT   ! a.u./a.u. -> pmflib velocity
    EnergyConv       = PMF_HARTREE2KCL                 ! hartree -> kcal/mol
    ForceConv        = PMF_HARTREE2KCL/PMF_AU2A        ! hartree/au -> kcal/mol/A

    ntb = SYS_NTV

    ! init pmf core ---------------------------------

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_VEL_VERLET,latoms,ntb,lnst,ldt,0.0d0,ltemp)
    call pmf_pbc_set_box_from_lvectors(lbox)

    ! init mask subsystem
    call pmf_pbc_get_cbox(has_box,cbox_a,cbox_b,cbox_c)
    anres = 1
    call pmf_mask_topo_init(fnatoms,anres,has_box,cbox_a,cbox_b,cbox_c)

    ! init mask topology atom masses and positions
    do i=1,fnatoms
        call pmf_mask_set_topo_atom_mcrd(i,lm(i),lx(1,i),lx(2,i),lx(3,i))
    end do

    ! finalize mask
    call pmf_mask_topo_finalize()

    ! print PMFLib header with system description
    call pmf_init_title('CP2K')

    ! load method setups and CVs definition
    call pmf_cp2k_process_control

    ! print main header
    call pmf_init_all(lm,lx,ControlPrmfile)

    return

end subroutine pmf_cp2k_init

!!===============================================================================
!! subroutine pmf_cp2k_shake
!!===============================================================================

!subroutine pmf_cp2k_shake(NSP,NAX,NSX,NA,TAUP,CVELP)

! use pmf_dat
! use pmf_sizes
! use pmf_cp2k_common
! use pmf_core_vv
! use pmf_cp2k_dat

! implicit none
! integer        :: NSP                  ! number of types
! integer        :: NAX                  ! max number of atom of same type
! integer        :: NSX                  ! max number of atom types
! integer        :: NA(NSX)              ! number of atoms of given type
! real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates r_u(t+dt)->R(t+dt)
! real(PMFDP)    :: CVELP(3,NAX,NSX)      ! velocities  v_u(t')->??
! ! -----------------------------------------------------------------------------

! if( .not. con_enabled ) return

! ! remap cp2k arrays to pmflib arrays
! call remap_coords_in(NSP,NAX,NSX,NA,TAUP,tmp_x)
! call remap_coords_in(NSP,NAX,NSX,NA,CVELP,tmp_v)

! ! do shake
! call pmf_core_vv_shake_SFR(tmp_x,tmp_v)

! ! remap arrays back
! call remap_coords_back(NSP,NAX,NSX,NA,TAUP,tmp_x)
! call remap_coords_back(NSP,NAX,NSX,NA,CVELP,tmp_v)

! return

!end subroutine pmf_cp2k_shake

!===============================================================================
! subroutine pmf_cp2k_force
!===============================================================================

subroutine pmf_cp2k_force(lx,lf)

    use pmf_sizes
    use pmf_core_vv
    use pmf_cp2k_dat

    implicit none
    real(PMFDP)    :: lx(:,:)         ! coordinates
    real(PMFDP)    :: lf(:,:)         ! force
    ! -----------------------------------------------
    real(PMFDP)    :: ene, epmf
    ! --------------------------------------------------------------------------

    ene = 0.0d0
    epmf = 0.0d0

    call pmf_core_vv_force_SRF(lx,lf,ene,epmf)

    return

end subroutine pmf_cp2k_force

!!===============================================================================
!! subroutine pmf_cp2k_rattle
!!===============================================================================

!subroutine pmf_cp2k_rattle(NSP,NAX,NSX,NA,CVELP)

! use pmf_dat
! use pmf_sizes
! use pmf_cp2k_common
! use pmf_core_vv
! use pmf_cp2k_dat

! implicit none
! integer        :: NSP                  ! number of types
! integer        :: NAX                  ! max number of atom of same type
! integer        :: NSX                  ! max number of atom types
! integer        :: NA(NSX)              ! number of atoms of given type
! real(PMFDP)    :: CVELP(3,NAX,NSX)      ! velocities
! ! -----------------------------------------------------------------------------

! if( .not. con_enabled ) return

! ! remap cp2k arrays to pmflib arrays
! call remap_coords_in(NSP,NAX,NSX,NA,CVELP,tmp_v)

! ! do rattle
! call pmf_core_vv_rattle_SFR(tmp_v)

! ! remap arrays back
! call remap_coords_back(NSP,NAX,NSX,NA,CVELP,tmp_v)

! return

!end subroutine pmf_cp2k_rattle

!===============================================================================
! subroutine pmf_cp2k_finalize
!===============================================================================

subroutine pmf_cp2k_finalize

    use pmf_constants
    use pmf_utils
    use pmf_finalize
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Library Finalization', '-')

    call pmf_timers_stop_timer(PMFLIB_TOTAL_TIMER)
    call pmf_finalize_all

end subroutine pmf_cp2k_finalize

!===============================================================================

end module pmf_cp2k


