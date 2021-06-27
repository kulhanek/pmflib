!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

module pmf_xdynbp

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_xdynbp_preinit(nsteps,stepsize,temp0)

    use xbp_topfile_dat
    use xbp_confile_dat
    use pmf_core
    use pmf_dat
    use pmf_timers
    use pmf_init
    use pmf_utils
    use pmf_pbc
    use pmf_mask

    implicit none
    integer             :: nsteps
    real(PMFDP)         :: stepsize
    real(PMFDP)         :: temp0
    ! --------------------------------------------
    integer             :: has_box
    real(PMFDP)         :: cbox(3)
    ! --------------------------------------------------------------------------

    ! setup conversion factors
    MassConv         = 1.0d0               ! g/mol -> g/mol
    LengthConv       = 1.0d0               ! A -> A
    AngleConv        = 1.0d0               ! rad -> rad
    TimeConv         = 1.0d0               ! fs -> fs
    VelocityConv     = 1.0d0               ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0               ! kcal/mol -> kcal/mol
    ForceConv        = -1.0d0              ! derivatives kcal/mol/A -> kcal/mol/A

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Library Pre-Initialization', '-')

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_LEAP_FROG,natom,0,nsteps,stepsize,0.0d0,temp0)

    ! init mask subsystem
    has_box = 0
    cbox(:) = 0.0d0
    call pmf_pbc_get_cbox(has_box,cbox)
    call pmf_mask_topo_init(natom,nres,has_box,cbox(1),cbox(2),cbox(3))

end subroutine pmf_xdynbp_preinit

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_xdynbp_init(ax,amass)

    use pmf_dat
    use pmf_core
    use pmf_utils
    use pmf_init
    use pmf_mask
    use pmf_control
    use xbp_topfile_dat

    use xbp_confile_dat

    implicit none
    real(PMFDP)     :: ax(:,:)
    real(PMFDP)     :: amass(:)
    ! --------------------------------------------
    integer         :: i
    ! ----------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Library Final Initialization', '-')

    ! init mask topology atom masses and positions
    do i=1,fnatoms
        call pmf_mask_set_topo_atom_mcrd(i,amass(i),ax(1,i),ax(2,i),ax(3,i))
    end do

    do i=1,nres
        call pmf_mask_set_topo_residue(i,res(i)%name,res(i)%start)
    end do

    ! finalize mask
    call pmf_mask_topo_finalize()

    ! print PMFLib header with system description
    call pmf_init_title('XDynBP')

    ! read pmf groups
    call pmf_control_read_pmflib_group(ControlPrmfile)

    ! read method CV setup
    call pmf_control_read_method_cvs_and_paths(ControlPrmfile)

    ! init all methods
    call pmf_init_all(amass,ax)

    ! init methods
    call pmf_init_pmf_methods

    if( cst_enabled .or. rst_enabled .or. mtd_enabled .or. abf_enabled .or. abp_enabled .or. mon_enabled ) then
       ! call pmf_init_profiling
    end if

    return

end subroutine pmf_xdynbp_init

!===============================================================================
! Subroutine: pmf_xdynbp_force
!===============================================================================

subroutine pmf_xdynbp_force(x,v,d,E)

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers
    use xbp_sizes
    use xbp_types

    implicit none
    real(DP)        :: x(:,:)
    real(DP)        :: v(:,:)
    real(DP)        :: d(:,:)
    type(ENERGIES)  :: E
    ! --------------------------------------------
    real(PMFDP)     :: epot
    real(PMFDP)     :: ekin
    real(PMFDP)     :: epmf
    real(PMFDP)     :: temp,bathtemp
    logical         :: remd_updated
    ! --------------------------------------------------------------------------

    epot = E%potential
    ekin = E%kinetic
    epmf = 0.0d0
    call pmf_timers_start_timer(PMFLIB_TIMER)
        call pmf_core_lf_update_xv(remd_updated,x,v,temp,bathtemp)
        call pmf_core_lf_force(x,v,d,epot,epmf)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    E%restraint%protein = E%restraint%protein + epmf
    E%restraint%total = E%restraint%total + epmf
    E%classical = E%classical + epmf
    E%potential = E%potential + epmf

    return

end subroutine pmf_xdynbp_force

!===============================================================================
! Subroutine: pmf_xdynbp_constraints
!===============================================================================

subroutine pmf_xdynbp_constraints(x)

    use pmf_dat
    use xbp_sizes
    use xbp_types
    use pmf_core_lf
    use pmf_timers

    implicit none
    real(DP)            :: x(:,:)  ! positions in t+dt
    ! -------------------------------------------------------------------

    if( .not. cst_enabled ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
        call pmf_core_lf_shake(x)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_xdynbp_constraints

!===============================================================================
! Subroutine: pmf_xdynbp_cst_checkatom
!===============================================================================

logical function pmf_xdynbp_cst_checkatom(atomid)

    use pmf_dat
    use cst_shake

    implicit none
    integer    :: atomid
    ! --------------------------------------------------------------------------

    pmf_xdynbp_cst_checkatom = .false.
    if( .not. cst_enabled ) return

    pmf_xdynbp_cst_checkatom = cst_shake_checkatom(atomid)

    return

end function pmf_xdynbp_cst_checkatom

!===============================================================================
! Subroutine:  pmf_xdynbp_num_of_cons
!===============================================================================

integer function pmf_xdynbp_num_of_cons()

    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    pmf_xdynbp_num_of_cons = NumOfCONs

end function pmf_xdynbp_num_of_cons

!===============================================================================
! Function:  pmf_xdynbp_allocate_shake
!===============================================================================

subroutine pmf_xdynbp_allocate_shake(num)

    use pmf_dat
    use cst_shake

    implicit none
    integer    :: num ! number of shake constraints
    ! --------------------------------------------------------------------------

    if( .not. cst_enabled ) return
    call cst_shake_allocate(num)

    return

end subroutine pmf_xdynbp_allocate_shake

!===============================================================================
! Function:  pmf_xdynbp_set_shake
!===============================================================================

subroutine pmf_xdynbp_set_shake(id,at1,at2,value)

    use pmf_dat
    use cst_shake

    implicit none
    integer        :: id       ! id of constraint
    integer        :: at1      ! id of first atom
    integer        :: at2      ! id of second atom
    real(PMFDP)    :: value    ! value of DS constraint
    ! --------------------------------------------------------------------------

    if( .not. cst_enabled ) return

    call cst_shake_set(id,at1,at2,value)

    return

end subroutine pmf_xdynbp_set_shake

!===============================================================================
! Subroutine: pmf_xdynbp_finalize
!===============================================================================

subroutine pmf_xdynbp_finalize

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

 return

end subroutine pmf_xdynbp_finalize

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_xdynbp
