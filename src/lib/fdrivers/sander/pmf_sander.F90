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

!sander.F90:   use pmf_sander
!sander.F90:    call pmf_sander_init_taskid_mpi(mytaskid)
!sander.F90:        call pmf_sander_init_preinit(mdin,natom,nres,ntb,ntc,nstlim,dt,temp0,a,b,c,alpha,beta,gamma)
!sander.F90:            call pmf_sander_set_residue(i,ih(m02-1+i),ix(i02-1+i))
!sander.F90:            call pmf_sander_set_atom(i,ih(m04-1+i),ih(m06-1+i))
!sander.F90:        call pmf_sander_finalize_preinit(natom,x(lmass),x(lcrd))
!sander.F90:        call pmf_sander_cst_init_collisions(ntc,nbonh,ix(iifstwt),ix(iibh),ix(ijbh),x(l50))
!sander.F90:        call pmf_sander_init(natom,x(lmass),x(lcrd))
!sander.F90:    call pmf_sander_bcast_dat_mpi(natom,numtasks,iparpt)
!sander.F90:        call pmf_sander_finalize()

! >> MD Loop
! a(t) <- dV/dr
! if( ntb > 0 ) call pmf_sander_update_box(a,b,c,alpha,beta,gamma): all at "t"
! call pmf_sander_force: r(t), v(t-dt/2), a(t), epot(t), ekin(t-dt)
! pressure coupling
! velocity update: v(t+dt/2), vold(t-dt/2)
! update position: r(t) -> r(t+dt)
! SHAKE + call pmf_sander_constraints: r(t+dt)
! fix velocities: from position changes: v(t+dt/2)
! thermostat + calculate KE: Ekin(t) <- [0.5( v(t+dt/2) + v(t-dt/2))]^2


module pmf_sander

use pmf_sizes

implicit none

! energy blob
type PMFSanderEnergy
    ! input
    real(PMFDP) :: epot         ! potential energy in t
    real(PMFDP) :: ekinvv       ! velocity-verlet kinetic energy in t-dt
    real(PMFDP) :: ekinlf       ! leapfrog kinetic energy in t-dt/2
    real(PMFDP) :: ekinv4       ! high-accuracy kinetic energy in t-2*dt
    ! output
    real(PMFDP) :: erst         ! PMF restrain energy
end type PMFSanderEnergy

contains

#ifdef MPI
!===============================================================================
! subroutine pmf_sander_init_taskid_mpi
!===============================================================================

subroutine pmf_sander_init_taskid_mpi(mytaskid)

    use pmf_init

    implicit none
    integer        :: mytaskid
    ! --------------------------------------------------------------------------

    call pmf_init_taskid_mpi(mytaskid,0)

end subroutine pmf_sander_init_taskid_mpi
!===============================================================================
#endif

!===============================================================================
! subroutine pmf_sander_init_preinit
!===============================================================================

subroutine pmf_sander_init_preinit(mdin,anatom,anres, &
                            antb,antc,ansteps,astepsize,atemp0, &
                            box_a,box_b,box_c,box_alpha,box_beta,box_gamma)

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_sander_dat
    use pmf_sander_control
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers

    implicit none
    character(*)   :: mdin
    integer        :: anatom                       ! number of atoms in AMBER topology
    integer        :: anres                        ! number of residues in AMBER topology
    integer        :: antb                         ! BOX type
    integer        :: antc                         ! shake mode
    integer        :: ansteps                      ! number of MD steps
    real(PMFDP)    :: astepsize                    ! step size
    real(PMFDP)    :: atemp0                       ! temperature
    real(PMFDP)    :: box_a,box_b,box_c            ! box dimmensions
    real(PMFDP)    :: box_alpha,box_beta,box_gamma
    ! -----------------------------------------------
    integer        :: has_box
    real(PMFDP)    :: cbox(3)
    ! --------------------------------------------------------------------------

    ! setup conversion factors
    MassConv         = 1.0d0        ! g/mol -> g/mol
    LengthConv       = 1.0d0        ! A -> A
    AngleConv        = PMF_D2R      ! deg -> rad
    TimeConv         = 1000.0d0     ! ps -> fs
    VelocityConv     = 1.0d0        ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0        ! kcal/mol -> kcal/mol
    ForceConv        = 1.0d0        ! kcal/mol/A -> kcal/mol/A

    ControlFileName  = mdin

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_LEAP_FROG,anatom,antb,ansteps,astepsize,0.0d0,atemp0)
    call pmf_pbc_set_box(box_a,box_b,box_c,box_alpha,box_beta,box_gamma)

    ! init mask subsystem
    call pmf_pbc_get_cbox(has_box,cbox)
    call pmf_mask_topo_init(anatom,anres,has_box,cbox(1),cbox(2),cbox(3))

    ! rewrite the setup
    fcanexmdloop     = .true.       ! the client is able to terminate md loop

    ! SHAKE
    fshake = .false.
    if( antc .ne. 1 ) then
        fshake = .true.
    end if

    return

end subroutine pmf_sander_init_preinit

!===============================================================================
! subroutine pmf_sander_set_residue
!===============================================================================

subroutine pmf_sander_set_residue(idx,name,first_atom)

    use pmf_mask

    implicit none
    integer            :: idx
    character(*)       :: name
    integer            :: first_atom
    ! --------------------------------------------------------------------------

    call pmf_mask_set_topo_residue(idx,name,first_atom)

end subroutine pmf_sander_set_residue

!===============================================================================
! subroutine pmf_sander_set_atom
!===============================================================================

subroutine pmf_sander_set_atom(idx,name,atype)

    use pmf_mask

    implicit none
    integer            :: idx
    character(*)       :: name
    character(*)       :: atype
    ! --------------------------------------------------------------------------

    call pmf_mask_set_topo_atom(idx,name,atype)

end subroutine pmf_sander_set_atom

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_finalize_preinit(anatom,amass,ax)

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_mask
    use pmf_core
    use pmf_init
    use pmf_sander_control
    use pmf_utils

    implicit none
    integer        :: anatom       ! number of atoms
    real(PMFDP)    :: amass(anatom)
    real(PMFDP)    :: ax(3,anatom)
    ! -----------------------------------------------
    integer        :: i, alloc_failed
    ! --------------------------------------------------------------------------

    ! init mask topology atom masses and positions
    allocate(frmass(fnatoms), stat= alloc_failed )
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate frmass!')
    end if

    do i=1,fnatoms
        frmass(i) = amass(i)
        call pmf_mask_set_topo_atom_mcrd(i,amass(i),ax(1,i),ax(2,i),ax(3,i))
    end do

    ! finalize mask
    call pmf_mask_topo_finalize()

    ! print PMFLib header with system description
    call pmf_init_title('SANDER')

    ! load method setups and CVs definition
    call pmf_sander_process_control

    return

end subroutine pmf_sander_finalize_preinit

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_init(anatom,amass,ax)

    use pmf_constants
    use pmf_sizes
    use pmf_init
    use pmf_dat
    use pmf_sander_dat

    implicit none
    integer        :: anatom       ! number of atoms
    real(PMFDP)    :: amass(anatom)
    real(PMFDP)    :: ax(3,anatom)
    ! --------------------------------------------------------------------------

    ! init all methods
    call pmf_init_all(amass,ax)

    ! init methods
    call pmf_init_pmf_methods

    return

end subroutine pmf_sander_init

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_shouldexit(exitcode)

    use pmf_constants
    use pmf_sizes
    use pmf_dat

    implicit none
    integer        :: exitcode       ! MD loop exit code
    ! --------------------------------------------------------------------------

    ! we cannot bcast the value here as pmf_sander_shouldexit is called on master
    ! protected sections several times per MD cycle
    ! thus fexit_mdloop is distributed in pmf_sander_force_mpi

    exitcode = fexit_mdloop

    return

end subroutine pmf_sander_shouldexit


#ifdef MPI
!===============================================================================
! subroutine pmf_sander_bcast_dat_mpi
!===============================================================================

subroutine pmf_sander_bcast_dat_mpi(anatom,anumtasks,aiparpt)

    use pmf_init
    use pmf_sander_dat
    use pmf_dat
    use pmf_utils

    integer     :: anatom
    integer     :: anumtasks
    integer     :: aiparpt(:)
    ! -----------------------------------------------
    integer     :: alloc_failed
    integer     :: i,j
    integer     :: istart,iend
    ! --------------------------------------------------------------------------

    call pmf_init_bcast_dat_mpi

    ! allocate temporary arrays
    if( fmaster ) then
        allocate(tmp_a(3,fnatoms), &
                 tmp_b(3,fnatoms), &
                 tmp_c(3,fnatoms),stat=alloc_failed)

         if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[SANDER] Unable to allocate memory for tmp_a/b/c arrays!')
         end if
    end if

    ! on all nodes !
    ! allocate owner map
    allocate(atm_owner_map(anatom),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory used in pmf_init_owner_map')
    end if

    ! init owner map
    do i=1,anumtasks
        istart = aiparpt(i)+1       ! aiparpt(1:anumtasks+1)
        iend = aiparpt(i+1)
        do j=istart,iend
            atm_owner_map(j) = i-1
        end do
    end do

end subroutine pmf_sander_bcast_dat_mpi

!===============================================================================
#endif

!===============================================================================
! subroutine pmf_sander_update_box
!===============================================================================

subroutine pmf_sander_update_box(a,b,c,alpha,beta,gamma)

    use pmf_dat
    use pmf_pbc
    use pmf_timers

    real(PMFDP)    :: a,b,c
    real(PMFDP)    :: alpha,beta,gamma
    ! --------------------------------------------------------------------------

    if( fsystype .eq. SYS_NT ) return  ! no box -> return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_pbc_set_box(a,b,c,alpha,beta,gamma)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_sander_update_box

!===============================================================================
! subroutine pmf_sander_force
!===============================================================================

subroutine pmf_sander_force(anatom,x,v,f,spmfene)

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers
    use pmf_dat

    implicit none
    integer                 :: anatom       ! number of atoms
    real(PMFDP)             :: x(3,anatom)  ! in
    real(PMFDP)             :: v(3,anatom)  ! in
    real(PMFDP)             :: f(3,anatom)  ! inout
    type(PMFSanderEnergy)   :: spmfene      ! inout
    ! --------------------------------------------------------------------------

    spmfene%erst = 0.0d0

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_update_step
    call pmf_core_lf_force(x,v,f,spmfene%epot,spmfene%ekinvv,spmfene%ekinlf,spmfene%ekinv4,spmfene%erst)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_sander_force

!===============================================================================
! Subroutine: pmf_sander_force_lng
!===============================================================================

subroutine pmf_sander_force_lng(anatom,flng)

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_timers

    implicit none
    integer        :: anatom            ! number of atoms
    real(PMFDP)    :: flng(3,anatom)    ! forces from Langevin
    ! --------------------------------------------------------------------------

    if( .not. lng_force_required ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_langevin_forces(flng)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_sander_force_lng

!===============================================================================
! subroutine pmf_sander_rstforce
!===============================================================================

subroutine pmf_sander_rstforce(anatom,x,f,epot,epmf)

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers

    implicit none
    integer        :: anatom        ! in - number of atoms
    real(PMFDP)    :: x(3,anatom)   ! in
    real(PMFDP)    :: f(3,anatom)   ! inout
    real(PMFDP)    :: epot          ! in
    real(PMFDP)    :: epmf          ! out
    ! --------------------------------------------------------------------------

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_rstforce(x,f,epot,epmf)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_sander_rstforce

!===============================================================================
! Subroutine: pmf_sander_constraints
!===============================================================================

subroutine pmf_sander_num_of_pmflib_cst(numofcst)

    use cst_dat

    implicit none
    real(PMFDP)    :: numofcst       ! number of CST constraints
    ! --------------------------------------------------------------------------

    numofcst = 0
    if( .not. cst_enabled ) return

    numofcst = NumOfCONs - NumOfSHAKECONs
    return

end subroutine pmf_sander_num_of_pmflib_cst

!===============================================================================
! Subroutine: pmf_sander_constraints
!===============================================================================

subroutine pmf_sander_constraints(leapfrog_mode,anatom,xbar,x,modified)

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_timers
    use pmf_utils

    implicit none
    integer        :: leapfrog_mode
    integer        :: anatom            ! number of atoms
    real(PMFDP)    :: xbar(3,anatom)    ! positions in t+dt - without shake for leapfrog_mode .eq. 1, otherwise position in t
    real(PMFDP)    :: x(3,anatom)       ! positions in t+dt
    logical        :: modified          ! was constraint applied?
    ! --------------------------------------------------------------------------

    modified = .false.
    if( .not. (cst_enabled .or. shake_force_required) ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_shake(x)
    if( shake_force_required .and. (fshake .or. cst_enabled) ) then
        if( leapfrog_mode .eq. 1 ) then
             call pmf_core_lf_shake_forces(xbar,x)
        else
            call pmf_utils_exit(PMF_OUT,1,'leapfrog_mode .eq. 1 is required in pmf_sander_constraints and shake_force_required')
        end if
    end if
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    modified = cst_enabled
    return

end subroutine pmf_sander_constraints

#ifdef MPI

!===============================================================================
! subroutine pmf_sander_force_mpi
!===============================================================================

subroutine pmf_sander_force_mpi(anatom,x,v,f,spmfene)

    use pmf_sizes
    use pmf_core_lf
    use pmf_sander_dat
    use pmf_dat
    use pmf_timers
    use pmf_utils
    use mpi

    implicit none
    integer                 :: anatom       ! number of atoms
    real(PMFDP)             :: x(3,anatom)  ! in
    real(PMFDP)             :: v(3,anatom)  ! in
    real(PMFDP)             :: f(3,anatom)  ! inout
    type(PMFSanderEnergy)   :: spmfene      ! inout
    ! ------------------------------------------------------
    integer         :: i,ierr
    ! --------------------------------------------------------------------------

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    ! gather data
    call pmf_sander_gather_array_mpi(tmp_a,x,atm_owner_map,1)
    call pmf_sander_gather_array_mpi(tmp_b,v,atm_owner_map,2)
    call pmf_sander_gather_array_mpi(tmp_c,f,atm_owner_map,3)

    if( fmaster ) then
        if( fdebug ) then
            write(PMF_DEBUG+fmytaskid,'(A)') '>> Input coordinates (force): '
            do i=1,NumOfLAtoms
                write(PMF_DEBUG+fmytaskid,'(3X,3F10.3)') tmp_a(:,RIndexes(i))
            end do
            write(PMF_DEBUG+fmytaskid,*)
        end if
        call pmf_core_lf_update_step
        call pmf_core_lf_force(tmp_a,tmp_b,tmp_c,spmfene%epot,spmfene%ekinvv,spmfene%ekinlf,spmfene%ekinv4,spmfene%erst)
    end if

    ! update data
    call pmf_sander_scatter_array_mpi(tmp_c,f,atm_owner_map,6)

    ! broadcast MD exit status
    call mpi_bcast(fexit_mdloop, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast fexit_mdloop variable!')
    end if

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

end subroutine pmf_sander_force_mpi

!===============================================================================
! Subroutine: pmf_sander_force_lng_mpi
!===============================================================================

subroutine pmf_sander_force_lng_mpi(anatom,flng)

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_timers

    implicit none
    integer        :: anatom            ! number of atoms
    real(PMFDP)    :: flng(3,anatom)    ! forces from Langevin
    ! --------------------------------------------------------------------------

    if( .not. lng_force_required ) return

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    call pmf_sander_gather_array_mpi(tmp_a,flng,atm_owner_map,1)
    call pmf_core_lf_langevin_forces(tmp_a)

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

    return

end subroutine pmf_sander_force_lng_mpi

!===============================================================================
! Subroutine: pmf_sander_constraints_mpi
!===============================================================================

! subroutine pmf_sander_constraints(leapfrog_mode,anatom,xbar,x,modified)

subroutine pmf_sander_constraints_mpi(leapfrog_mode,anatom,xbar,x,modified)

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_sander_dat
    use pmf_dat
    use pmf_timers
    use pmf_utils

    implicit none
    integer         :: leapfrog_mode
    integer         :: anatom
    real(PMFDP)     :: xbar(3,anatom)       ! positions in t+dt - without shake for leapfrog_mode .eq. 1, otherwise position in t
    real(PMFDP)     :: x(3,anatom)          ! positions in t+dt
    logical         :: modified             ! was constraint applied?
    ! ------------------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    modified = .false.
    if( .not. (cst_enabled .or. shake_force_required) ) return

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    ! gather data
    call pmf_sander_gather_array_mpi(tmp_a,x,atm_owner_map,1)

    if( fmaster ) then
        if( fdebug ) then
            write(PMF_DEBUG+fmytaskid,'(A)') '>> Input coordinates (shake): '
            do i=1,NumOfLAtoms
                write(PMF_DEBUG+fmytaskid,'(3X,3F10.3)') tmp_a(:,RIndexes(i))
            end do
            write(PMF_DEBUG+fmytaskid,*)
        end if
        call pmf_core_lf_shake(tmp_a)
    end if

    ! update data
    call pmf_sander_scatter_array_mpi(tmp_a,x,atm_owner_map,4)

    if( shake_force_required .and. (fshake .or. cst_enabled) ) then
        if( leapfrog_mode .eq. 1 ) then

            ! gather data
            call pmf_sander_gather_array_mpi(tmp_b,xbar,atm_owner_map,5)

            if( fmaster ) then
            end if
        else
            call pmf_utils_exit(PMF_OUT,1,'leapfrog_mode .eq. 1 is required in pmf_sander_constraints_mpi')
        end if
    end if

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

    modified = cst_enabled

end subroutine pmf_sander_constraints_mpi
!===============================================================================
#endif

!===============================================================================
! Subroutine: pmf_sander_cst_init_collisions
!===============================================================================

subroutine pmf_sander_cst_init_collisions(ntc,nbt,ifstwt,ib,jb,conp)

    use pmf_dat
    use cst_shake
    use pmf_sizes
    use pmf_core
    use pmf_utils

    implicit none
    integer        :: ntc          !
    integer        :: nbt          ! number of X-H bonds
    integer        :: ifstwt(*)    ! determine if bond is part of water
    integer        :: ib(*)        ! the first atom of bond
    integer        :: jb(*)        ! the second atom of bond
    real(PMFDP)    :: conp(*)
    ! -----------------------------------------------
    integer        :: ll, i, j, num
    ! -----------------------------------------------------------------------------

    if( .not. cst_enabled ) return
    if( .not. fmaster ) return ! only master can init shake constraint in collision
    if( ntc .eq. 1 ) return    ! no SHAKE

    if( ntc .ne. 2 ) then
        call pmf_utils_exit(PMF_OUT,1,'ntc has to be either one or two for PMFLib constrained dynamics!')
    end if

    ! determine number of SHAKE constraints in collision
    num = 0

    do ll = 1,nbt
        i = ib(ll)/3+1
        j  = jb(ll)/3+1
        if( ifstwt(ll) == 1 ) then
            if( cst_shake_checkatom(i) .or. cst_shake_checkatom(j) ) then
                call pmf_utils_exit(PMF_OUT,1,'Water atom cannot be a part of PMFLib constraint!')
            end if
        end if
        if( (cst_shake_checkatom(i) .or. cst_shake_checkatom(j)) .and. &
            (.not. (cst_shake_checkatom(i) .and. cst_shake_checkatom(j))) ) then
            num = num + 1
            cycle
        end if
    end do

    if( num .eq. 0 ) return    ! nothing in collision

    ! set SHAKE constraints in collisions
    call cst_shake_allocate(num)

    num = 1
    do ll = 1,nbt
        if (ifstwt(ll) == 1) cycle
        i  = ib(ll)/3+1
        j  = jb(ll)/3+1
        if( (cst_shake_checkatom(i) .or. cst_shake_checkatom(j)) .and. &
            (.not. (cst_shake_checkatom(i) .and. cst_shake_checkatom(j))) ) then
            call cst_shake_set(num,i,j,conp(ll))
            num = num + 1
            cycle
        end if
    end do

    return

end subroutine pmf_sander_cst_init_collisions

!===============================================================================
! Subroutine: pmf_sander_cst_checkatom
!===============================================================================

logical function pmf_sander_cst_checkatom(atomid)

    use pmf_dat
    use cst_shake

    implicit none
    integer    :: atomid
    ! --------------------------------------------------------------------------

    pmf_sander_cst_checkatom = .false.
    if( .not. cst_enabled ) return

    pmf_sander_cst_checkatom = cst_shake_checkatom(atomid)

    return

end function pmf_sander_cst_checkatom

!===============================================================================
! subroutine pmf_sander_finalize
!===============================================================================

subroutine pmf_sander_finalize

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

end subroutine pmf_sander_finalize

!===============================================================================

#ifdef MPI
!===============================================================================
! Subroutine:  pmf_sander_gather_array_mpi
!===============================================================================

subroutine pmf_sander_gather_array_mpi(lx,x,atm_owner_map,id_offset)

    use mpi
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)       :: lx(:,:)  ! local
    real(PMFDP)       :: x(:,:)   ! scattered data
    integer           :: atm_owner_map(:)
    integer           :: id_offset
    ! ------------------------------------------------------
    integer           :: i,at_id,ow_id,ierr,mpi_id
    integer           :: status(MPI_STATUS_SIZE)
    real(PMFDP)       :: tmp_x(3)
    ! --------------------------------------------------------------------------

    ierr = MPI_SUCCESS

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> gather_array'
    end if

    do i=1,NumOfLAtoms
        at_id = RIndexes(i)
        mpi_id = at_id + id_offset*fnatoms
        ow_id = atm_owner_map(at_id)
        if( ow_id .eq. fmytaskid ) then
            if( fmaster ) then
                ! get local copy
                lx(:,at_id) = x(:,at_id)
                if( fdebug ) then
                    write(PMF_DEBUG+fmytaskid,'(A,I5,1X,3F10.3)') '   local atom: ',at_id,lx(:,at_id)
                end if
            else
                ! send data
                tmp_x(:) = x(:,at_id)
                if( fdebug ) then
                    write(PMF_DEBUG+fmytaskid,'(A,I5,1X,3F10.3)') '   send atom:  ',at_id,tmp_x(:)
                end if
                call mpi_send(tmp_x(1),3,mpi_double_precision,0,mpi_id,mpi_comm_world,ierr)
            end if
        else
            if( fmaster ) then
                call mpi_recv(lx(1,at_id),3,mpi_double_precision,ow_id,mpi_id,mpi_comm_world,status,ierr)
                if( fdebug ) then
                    write(PMF_DEBUG+fmytaskid,'(A,I5,1X,3F10.3)') '   recv atom:  ',at_id,lx(:,at_id)
                end if
            end if
        end if
        if( ierr .ne. MPI_SUCCESS ) then
            call pmf_utils_exit(PMF_OUT,1,'Unable to transfer coords data via MPI (gather_array)!')
        end if
    end do

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*)
    end if

end subroutine pmf_sander_gather_array_mpi

!===============================================================================
! Subroutine:  pmf_sander_scatter_array_mpi
!===============================================================================

subroutine pmf_sander_scatter_array_mpi(lx,x,atm_owner_map,id_offset)

    use mpi
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)       :: lx(:,:)
    real(PMFDP)       :: x(:,:)
    integer           :: atm_owner_map(:)
    integer           :: id_offset
    ! ------------------------------------------------------
    integer           :: i,at_id,ow_id,ierr,mpi_id
    integer           :: status(MPI_STATUS_SIZE)
    real(PMFDP)       :: tmp_x(3)
    ! --------------------------------------------------------------------------

    ierr = MPI_SUCCESS

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> scatter_array'
    end if

    do i=1,NumOfLAtoms
        at_id = RIndexes(i)
        mpi_id = at_id + id_offset*fnatoms
        ow_id = atm_owner_map(at_id)
        if( ow_id .eq. fmytaskid ) then
            if( fmaster ) then
                if( fdebug ) then
                    write(PMF_DEBUG+fmytaskid,'(A,I5,1X,3F10.3)') '   local atom: ',at_id,lx(:,at_id)
                end if
                ! set local copy
                x(:,at_id) = lx(:,at_id)
            else
                ! recv data
                call mpi_recv(x(1,at_id),3,mpi_double_precision,0,mpi_id,mpi_comm_world,status,ierr)
                if( fdebug ) then
                    write(PMF_DEBUG+fmytaskid,'(A,I5,1X,3F10.3)') '   recv atom:  ',at_id,x(:,at_id)
                end if
            end if
        else
            if( fmaster ) then
                tmp_x(:) = lx(:,at_id)
                if( fdebug ) then
                    write(PMF_DEBUG+fmytaskid,'(A,I5,1X,3F10.3)') '   send atom:  ',at_id,tmp_x(1)
                end if
                call mpi_send(tmp_x(1),3,mpi_double_precision,ow_id,mpi_id,mpi_comm_world,ierr)
            end if
        end if
        if( ierr .ne. MPI_SUCCESS ) then
            call pmf_utils_exit(PMF_OUT,1,'Unable to transfer coords data via MPI (scatter_array)!')
        end if
    end do

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*)
    end if

end subroutine pmf_sander_scatter_array_mpi
!===============================================================================

#endif

end module pmf_sander




