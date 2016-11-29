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

module pmf_pmemd

implicit none
contains

#ifdef MPI
!===============================================================================
! subroutine pmf_pmemd_init_taskid_mpi
!===============================================================================

subroutine pmf_pmemd_init_taskid_mpi(mytaskid)

    use pmf_init

    implicit none
    integer        :: mytaskid
    ! --------------------------------------------------------------------------

    call pmf_init_taskid_mpi(mytaskid,0)

end subroutine pmf_pmemd_init_taskid_mpi
!===============================================================================
#endif

!===============================================================================
! subroutine pmf_pmemd_init_preinit
!===============================================================================

subroutine pmf_pmemd_init_preinit(mdin,anatom,anres, &
                            antb,ansteps,astepsize,atemp0, &
                            box_a,box_b,box_c,box_alpha,box_beta,box_gamma)

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_utils
    use pmf_pmemd_dat
    use pmf_pmemd_control
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers

    implicit none
    character(*)   :: mdin
    integer        :: anatom                       ! number of atoms in AMBER topology
    integer        :: anres                        ! number of residues in AMBER topology
    integer        :: antb                         ! BOX type
    integer        :: ansteps                      ! number of MD steps
    real(PMFDP)    :: astepsize                    ! step size
    real(PMFDP)    :: atemp0                       ! temperature
    real(PMFDP)    :: box_a,box_b,box_c            ! box dimmensions
    real(PMFDP)    :: box_alpha,box_beta,box_gamma
    ! ------------------------------------------------------
    integer        :: has_box
    real(PMFDP)    :: cbox(3)
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    ! setup conversion factors
    MassConv         = 1.0d0               ! g/mol -> g/mol
    LengthConv       = 1.0d0               ! A -> A
    AngleConv        = PMF_D2R             ! deg -> rad
    TimeConv         = 1000.0d0            ! ps -> fs
    VelocityConv     = 1.0d0               ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0               ! kcal/mol -> kcal/mol
    ForceConv        = 1.0d0               ! kcal/mol/A -> kcal/mol/A

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

end subroutine pmf_pmemd_init_preinit

!===============================================================================
! subroutine pmf_pmemd_set_residue
!===============================================================================

subroutine pmf_pmemd_set_residue(idx,name,first_atom)

    use pmf_mask
    use pmf_dat

    implicit none
    integer            :: idx
    character(*)       :: name
    integer            :: first_atom
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_mask_set_topo_residue(idx,name,first_atom)

end subroutine pmf_pmemd_set_residue

!===============================================================================
! subroutine pmf_pmemd_set_atom
!===============================================================================

subroutine pmf_pmemd_set_atom(idx,name,atype)

    use pmf_mask

    implicit none
    integer            :: idx
    character(*)       :: name
    character(*)       :: atype
    ! --------------------------------------------------------------------------

    call pmf_mask_set_topo_atom(idx,name,atype)

end subroutine pmf_pmemd_set_atom

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_finalize_preinit(amass,ax)

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_mask
    use pmf_core
    use pmf_init
    use pmf_pmemd_control
    use pmf_utils

    implicit none
    real(PMFDP)    :: amass(:)
    real(PMFDP)    :: ax(:,:)
    ! ------------------------------------------------------
    integer        :: i, alloc_failed
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

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
    call pmf_init_title('PMEMD')

    ! load method setups and CVs definition
    call pmf_pmemd_process_control

end subroutine pmf_pmemd_finalize_preinit

#ifdef MPI
!===============================================================================
! subroutine pmf_pmemd_bcast_constraints_mpi
!===============================================================================

subroutine pmf_pmemd_bcast_constraints_mpi

    use cst_init
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) 'PMFLib MPI: TaskID: ', fmytaskid
    end if

    ! this HAS TO BE called by all processes
    call cst_init_mpi_bcast_constraints

end subroutine pmf_pmemd_bcast_constraints_mpi
#endif

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_init(amass,ax)

    use pmf_constants
    use pmf_sizes
    use pmf_init
    use pmf_dat
    use cst_init
    use pmf_pmemd_control

    implicit none
    real(PMFDP)    :: amass(:)
    real(PMFDP)    :: ax(:,:)
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    ! init PMF subsystems
    call pmf_init_all(amass,ax)

    ! init methods
    call pmf_init_pmf_methods

    return

end subroutine pmf_pmemd_init

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_shouldexit(exitcode)

    use pmf_constants
    use pmf_sizes
    use pmf_dat

    implicit none
    integer        :: exitcode       ! MD loop exit code
    ! --------------------------------------------------------------------------

    ! we cannot bcast the value here as pmf_pmemd_shouldexit is called on master
    ! protected sections several times per MD cycle
    ! thus fexit_mdloop is distributed in pmf_pmemd_force_mpi

    exitcode = fexit_mdloop

    return

end subroutine pmf_pmemd_shouldexit


#ifdef MPI
!===============================================================================
! subroutine pmf_pmemd_bcast_dat_mpi
!===============================================================================

subroutine pmf_pmemd_bcast_dat_mpi

    use pmf_init
    use pmf_pmemd_dat
    use pmf_dat
    use pmf_utils

    integer :: alloc_failed
    ! --------------------------------------------------------------------------

    call pmf_init_bcast_dat_mpi

    ! allocate temporary arrays
    if( fmaster ) then
        allocate(tmp_a(3,fnatoms), &
                 tmp_b(3,fnatoms), &
                 tmp_c(3,fnatoms),stat=alloc_failed)

         if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[PMEMD] Unable to allocate memory for tmp_a/b/c arrays!')
         end if
    end if

end subroutine pmf_pmemd_bcast_dat_mpi

!===============================================================================
#endif

!===============================================================================
! subroutine pmf_pmemd_update_box
!===============================================================================

subroutine pmf_pmemd_update_box(a,b,c,alpha,beta,gamma)

    use pmf_dat
    use pmf_pbc

    real(PMFDP)    :: a,b,c
    real(PMFDP)    :: alpha,beta,gamma
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fsystype .eq. SYS_NT ) return  ! no box -> return

    call pmf_pbc_set_box(a,b,c,alpha,beta,gamma)

end subroutine pmf_pmemd_update_box

!===============================================================================
! subroutine pmf_pmemd_update_xv
!===============================================================================

subroutine pmf_pmemd_update_xv(updated,x,v,temp,bathtemp)

    use pmf_sizes
    use pmf_core_lf
    use pmf_dat
    use pmf_timers

    implicit none
    logical        :: updated
    real(PMFDP)    :: x(:,:)
    real(PMFDP)    :: v(:,:)
    real(PMFDP)    :: temp
    real(PMFDP)    :: bathtemp
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_update_xv(updated,x,v,temp,bathtemp)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_pmemd_update_xv

!===============================================================================
! subroutine pmf_pmemd_force
!===============================================================================

subroutine pmf_pmemd_force(x,v,f,epot,epmf)

    use pmf_sizes
    use pmf_core_lf
    use pmf_dat
    use pmf_timers

    implicit none
    real(PMFDP)    :: x(:,:)
    real(PMFDP)    :: v(:,:)
    real(PMFDP)    :: f(:,:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: epmf
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_force(x,v,f,epot,epmf)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_pmemd_force

!===============================================================================
! Subroutine: pmf_pmemd_constraints
!===============================================================================

subroutine pmf_pmemd_constraints(x,modified)

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_timers

    implicit none
    real(PMFDP)    :: x(:,:)     ! positions in t+dt
    logical        :: modified   ! was constraint applied?
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    modified = .false.
    if( .not. cst_enabled ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_shake(x)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    modified = .true.

end subroutine pmf_pmemd_constraints

#ifdef MPI
!===============================================================================
! subroutine pmf_pmemd_update_xv_mpi
!===============================================================================

subroutine pmf_pmemd_update_xv_mpi(updated,x,v,temp,bathtemp,atm_owner_map)

    use pmf_sizes
    use pmf_core_lf
    use pmf_pmemd_dat
    use pmf_utils
    use pmf_dat
    use pmf_timers
    use mpi

    implicit none
    logical        :: updated
    real(PMFDP)    :: x(:,:)
    real(PMFDP)    :: v(:,:)
    real(PMFDP)    :: temp
    real(PMFDP)    :: bathtemp
    integer        :: atm_owner_map(:) ! atom map among processes
    ! ------------------------------------------------------
    integer        :: ierr,i
    ! --------------------------------------------------------------------------

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    ! gather data
    call pmf_pmemd_gather_array_mpi(tmp_a,x,atm_owner_map,1)
    call pmf_pmemd_gather_array_mpi(tmp_b,v,atm_owner_map,2)

    ! do main staff
    if( fmaster ) then
        if( fdebug ) then
            write(PMF_DEBUG+fmytaskid,'(A)') '>> Input coordinates (update_xv): '
            do i=1,NumOfLAtoms
                write(PMF_DEBUG+fmytaskid,'(3X,3F10.3)') tmp_a(:,RIndexes(i))
            end do
            write(PMF_DEBUG+fmytaskid,*)
        end if
        call pmf_core_lf_update_xv(updated,tmp_a,tmp_b,temp,bathtemp)
    end if

    ! broadcast updated
    call mpi_bcast(updated, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to broadcast the value of updated!')
    end if

    ! update data if necessary
    if( updated ) then
        call pmf_pmemd_scatter_array_mpi(tmp_a,x,atm_owner_map,3)
        call pmf_pmemd_scatter_array_mpi(tmp_b,v,atm_owner_map,4)
    end if

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

end subroutine pmf_pmemd_update_xv_mpi

!===============================================================================
! subroutine pmf_pmemd_force_mpi
!===============================================================================

subroutine pmf_pmemd_force_mpi(x,v,f,epot,epmf,atm_owner_map)

    use pmf_sizes
    use pmf_core_lf
    use pmf_pmemd_dat
    use pmf_dat
    use pmf_timers
    use pmf_utils
    use mpi

    implicit none
    real(PMFDP)    :: x(:,:)
    real(PMFDP)    :: v(:,:)
    real(PMFDP)    :: f(:,:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: epmf
    integer        :: atm_owner_map(:) ! atom map among processes
    ! ------------------------------------------------------
    integer        :: i, ierr
    ! --------------------------------------------------------------------------

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    ! gather data
    call pmf_pmemd_gather_array_mpi(tmp_a,x,atm_owner_map,1)
    call pmf_pmemd_gather_array_mpi(tmp_b,v,atm_owner_map,2)
    call pmf_pmemd_gather_array_mpi(tmp_c,f,atm_owner_map,3)

    if( fmaster ) then
        if( fdebug ) then
            write(PMF_DEBUG+fmytaskid,'(A)') '>> Input coordinates (force): '
            do i=1,NumOfLAtoms
                write(PMF_DEBUG+fmytaskid,'(3X,3F10.3)') tmp_a(:,RIndexes(i))
            end do
            write(PMF_DEBUG+fmytaskid,*)
        end if

        call pmf_core_lf_force(tmp_a,tmp_b,tmp_c,epot,epmf)
    end if

    ! broadcast MD exit status
    call mpi_bcast(fexit_mdloop, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast fexit_mdloop variable!')
    end if

    ! update data
    call pmf_pmemd_scatter_array_mpi(tmp_a,x,atm_owner_map,4)
    call pmf_pmemd_scatter_array_mpi(tmp_b,v,atm_owner_map,5)
    call pmf_pmemd_scatter_array_mpi(tmp_c,f,atm_owner_map,6)

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

end subroutine pmf_pmemd_force_mpi

!===============================================================================
! Subroutine: pmf_pmemd_constraints_mpi
!===============================================================================

subroutine pmf_pmemd_constraints_mpi(x,modified,atm_owner_map)

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_pmemd_dat
    use pmf_dat
    use pmf_timers

    implicit none
    real(PMFDP)    :: x(:,:)               ! positions in t+dt
    logical        :: modified             ! was constraint applied?
    integer        :: atm_owner_map(:)     ! atom map among processes
    ! ------------------------------------------------------
    integer        :: i
    ! --------------------------------------------------------------------------

    modified = .false.
    if( .not. cst_enabled ) return

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    ! gather data
    call pmf_pmemd_gather_array_mpi(tmp_a,x,atm_owner_map,1)

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
    call pmf_pmemd_scatter_array_mpi(tmp_a,x,atm_owner_map,4)

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

    modified = .true.

end subroutine pmf_pmemd_constraints_mpi
!===============================================================================
#endif

!===============================================================================
! Subroutine: pmf_pmemd_cst_checkatom
!===============================================================================

logical function pmf_pmemd_cst_checkatom(atomid)

    use pmf_dat
    use cst_shake

    implicit none
    integer    :: atomid
    ! --------------------------------------------------------------------------

    pmf_pmemd_cst_checkatom = .false.
    if( .not. cst_enabled ) return

    pmf_pmemd_cst_checkatom = cst_shake_checkatom(atomid)

end function pmf_pmemd_cst_checkatom

!===============================================================================
! Subroutine: pmf_pmemd_cst_shake_allocate
!===============================================================================

subroutine pmf_pmemd_cst_shake_allocate(num)

    use pmf_dat
    use cst_shake

    implicit none
    integer    :: num ! number of shake constraints
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return
    if( .not. cst_enabled ) return

    call cst_shake_allocate(num)

end subroutine pmf_pmemd_cst_shake_allocate

!===============================================================================
! Function:  pmf_pmemd_cst_set_shake
!===============================================================================

subroutine pmf_pmemd_cst_set_shake(id,at1,at2,value)

    use pmf_sizes
    use pmf_dat
    use cst_shake

    implicit none
    integer        :: id       ! id of constraint
    integer        :: at1      ! id of first atom
    integer        :: at2      ! id of second atom
    real(PMFDP)    :: value    ! value of DS constraint
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return
    if( .not. cst_enabled ) return

    call cst_shake_set(id,at1,at2,value)

end subroutine pmf_pmemd_cst_set_shake

!===============================================================================
! subroutine pmf_pmemd_finalize
!===============================================================================

subroutine pmf_pmemd_finalize

    use pmf_constants
    use pmf_utils
    use pmf_finalize
    use pmf_timers
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Library Finalization', '-')

    call pmf_timers_stop_timer(PMFLIB_TOTAL_TIMER)
    call pmf_finalize_all(.true.)

end subroutine pmf_pmemd_finalize

#ifdef MPI
!===============================================================================
! Subroutine:  gather_array_mpi
!===============================================================================

subroutine pmf_pmemd_gather_array_mpi(lx,x,atm_owner_map,id_offset)

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

end subroutine pmf_pmemd_gather_array_mpi

!===============================================================================
! Subroutine:  scatter_array_mpi
!===============================================================================

subroutine pmf_pmemd_scatter_array_mpi(lx,x,atm_owner_map,id_offset)

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

end subroutine pmf_pmemd_scatter_array_mpi
!===============================================================================

#endif

!===============================================================================

end module pmf_pmemd




