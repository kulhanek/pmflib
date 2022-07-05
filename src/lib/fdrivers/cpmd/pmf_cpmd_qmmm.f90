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
! subroutine pmf_cpmd_init_qmmm
!===============================================================================

subroutine pmf_cpmd_init_qmmm(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,PMA0,TAU0,NOMORE, &
                         DELT_IONS,TEMPW,NTB,BOX,BETA)

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
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)            ! atom indexes
    integer        :: CPSP(NTAT)            ! atom types
    real(PMFDP)    :: PMA0(NSX)            ! masses
    real(PMFDP)    :: TAU0(3,NAX,NSX)      ! coordinates
    integer        :: NOMORE               ! number of MD steps
    real(PMFDP)    :: DELT_IONS            ! time step in atomic unit
    real(PMFDP)    :: TEMPW                ! temperature in K
    integer        :: NTB                  ! box type
    real(PMFDP)    :: BOX(3)               ! box dimensions
    real(PMFDP)    :: BETA                 ! box angle
    ! -----------------------------------------------
    integer        :: alloc_failed
    real(PMFDP)    :: lbox(6)
    integer        :: has_box
    real(PMFDP)    :: cbox(3)
    integer        :: anres                        ! number of residues
    integer        :: i
    integer        :: stype
    real(PMFDP)    :: GroLengthConv,GroAngleConv
    ! --------------------------------------------------------------------------

    ! setup conversion factors ----------------------
    MassConv         = 1.0d0                           ! g/mol -> g/mol
    LengthConv       = PMF_AU2A                        ! a.u. -> A
    AngleConv        = 1.0d0                           ! rad -> rad
    TimeConv         = PMF_AU2FS                       ! a.u. -> fs
    VelocityConv     = PMF_AU2A/PMF_AU2FS*PMF_VDT2DT ! a.u./a.u. -> pmflib velocity
    EnergyConv       = PMF_HARTREE2KCL                 ! hartree -> kcal/mol
    ForceConv        = PMF_HARTREE2KCL/PMF_AU2A        ! hartree/A -> kcal/mol/A

    GroLengthConv       = 10.0d0                          ! nm -> A
    GroAngleConv        = PMF_D2R                         ! deg -> rad

    ! setup box -------------------------------------
    lbox(1) = BOX(1)*GroLengthConv
    lbox(2) = BOX(2)*GroLengthConv
    lbox(3) = BOX(3)*GroLengthConv
    lbox(4:6) = BETA*GroAngleConv

    if( NTB .lt. 0 ) then
    call pmf_utils_exit(PMF_OUT,1,'Unsupported box type!')
    end if

    stype = SYS_NTV
    if( NTB .eq. 0 )  stype = SYS_NT

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_VEL_VERLET,NTAT-NDAT,ntb,NOMORE,DELT_IONS,0.0d0,TEMPW)
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
    call remap_qmmm_coords_in(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAU0,tmp_x)
    call remap_qmmm_masses_in(NTAT,NSAT,NDAT,NSX,CPSP,PMA0,tmp_mass)

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
    call pmf_init_title('CPMD-QMMM')

    ! load method setups
    call pmf_cpmd_process_control_begin()

    ! print main header
    call pmf_init_all(tmp_mass,tmp_x)

    if( NDAT .ne. 0 ) then
        write(PMF_OUT,'(A)')    '#'
        write(PMF_OUT,'(A)')    '# INFO: Dummy atoms detected (probably due to capping).'
        write(PMF_OUT,'(A)')    '#'
        write(PMF_OUT,'(A,i6)') '# Number of solute atoms       : ',NSAT
        write(PMF_OUT,'(A,i6)') '# Number of solvent atoms      : ',NTAT-NDAT-NSAT
        write(PMF_OUT,'(A,i6)') '# Total number of atoms        : ',NTAT-NDAT
        write(PMF_OUT,'(A)')    '# --------------------------------------'
        write(PMF_OUT,'(A,i6)') '# Number of dummy atoms        : ',NDAT
        write(PMF_OUT,'(A,i6)') '# Total number of system atoms : ',NTAT
    else
        write(PMF_OUT,'(A)')    '#'
        write(PMF_OUT,'(A)')    '# INFO: No dummy atoms detected.'
        write(PMF_OUT,'(A)')    '#'
        write(PMF_OUT,'(A,i6)') '# Number of solute atoms       : ',NSAT
        write(PMF_OUT,'(A,i6)') '# Number of solvent atoms      : ',NTAT-NSAT
        write(PMF_OUT,'(A,i6)') '# Total number of atoms        : ',NTAT
    end if
    write(PMF_OUT,'(A)')

    ! finalize input file processing
    call pmf_cpmd_process_control_end()

    ! end of pmflib initialization
    call pmf_init_pmf_methods()

    ! rewrite the setup
    fcanexmdloop     = .true.       ! the client is able to terminate md loop

    return

end subroutine pmf_cpmd_init_qmmm

!===============================================================================
! subroutine pmf_cpmd_shake_qmmm
!===============================================================================

subroutine pmf_cpmd_shake_qmmm(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,CVELP)

    use pmf_dat
    use pmf_sizes
    use pmf_cpmd_common
    use pmf_core_vv
    use pmf_cpmd_dat

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)           ! atom indexes
    integer        :: CPSP(NTAT)           ! atom types
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: CVELP(3,NAX,NSX)     ! velocities
    ! --------------------------------------------------------------------------

    if( .not. cst_enabled ) return

    ! remap cpmd arrays to pmflib arrays
    call remap_qmmm_coords_in(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,tmp_x)
    call remap_qmmm_coords_in(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,CVELP,tmp_v)

    ! do shake
    call pmf_core_vv_shake_SFR(tmp_x,tmp_v)

    ! remap arrays back
    call remap_qmmm_coords_back(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,tmp_x)
    call remap_qmmm_coords_back(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,CVELP,tmp_v)

    return

end subroutine pmf_cpmd_shake_qmmm

!===============================================================================
! subroutine pmf_cpmd_force_qmmm
!===============================================================================

subroutine pmf_cpmd_force_qmmm(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,FION)

    use pmf_sizes
    use pmf_core_vv
    use pmf_cpmd_dat
    use pmf_cpmd_common

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)            ! atom indexes
    integer        :: CPSP(NTAT)            ! atom types
    real(PMFDP)    :: TAUP(3,NAX,NSX)      ! coordinates
    real(PMFDP)    :: FION(3,NAX,NSX)      ! forces
    real(PMFDP)    :: ENE
    ! -----------------------------------------------
    real(PMFDP)    :: epmf
    ! --------------------------------------------------------------------------

    ! copy coordinates to tmp_x array
    call remap_qmmm_coords_in(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,TAUP,tmp_x)

    ! it is necessary to clean tmp_f because forces are added
    tmp_f(:,:) = 0.0d0
    call remap_qmmm_force_in_add(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,FION,tmp_f)

    epmf = 0.0d0

! FIXME
!    call pmf_core_vv_force_SRF(tmp_x,tmp_f,ENE,epmf)

    ! convert forces back ----------------------------------
    call remap_qmmm_force_back(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,FION,tmp_f)

    return

end subroutine pmf_cpmd_force_qmmm

!===============================================================================
! subroutine pmf_cpmd_rattle_qmmm
!===============================================================================

subroutine pmf_cpmd_rattle_qmmm(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,CVELP)

    use pmf_dat
    use pmf_sizes
    use pmf_cpmd_common
    use pmf_core_vv
    use pmf_cpmd_dat

    implicit none
    integer        :: NTAT                 ! total number of atoms
    integer        :: NSAT                 ! number of solute atoms
    integer        :: NDAT                 ! number of dummy atoms
    integer        :: NAX                  ! max number of atom of same type
    integer        :: NSX                  ! max number of atom types
    integer        :: CPAT(NTAT)           ! atom indexes
    integer        :: CPSP(NTAT)           ! atom types
    real(PMFDP)    :: CVELP(3,NAX,NSX)     ! velocities
    ! --------------------------------------------------------------------------

    if( .not. cst_enabled ) return

    ! remap cpmd arrays to pmflib arrays
    call remap_qmmm_coords_in(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,CVELP,tmp_v)

    ! do rattle
    call pmf_core_vv_rattle_SFR(tmp_v)

    ! remap arrays back
    call remap_qmmm_coords_back(NTAT,NSAT,NDAT,NAX,NSX,CPAT,CPSP,CVELP,tmp_v)

    return

end subroutine pmf_cpmd_rattle_qmmm

!===============================================================================




