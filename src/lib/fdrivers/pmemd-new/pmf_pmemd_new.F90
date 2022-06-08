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

module pmf_pmemd_new

implicit none
contains

#ifdef MPI
!===============================================================================
! subroutine pmf_pmemd_init_taskid_mpi
!===============================================================================

subroutine pmf_pmemd_init_taskid_mpi(mytaskid,numoftasks)

    use pmf_init

    implicit none
    integer        :: mytaskid
    integer        :: numoftasks
    ! --------------------------------------------------------------------------

    call pmf_init_taskid_mpi(mytaskid,numoftasks)

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
    use pmf_pmemd_new_dat
    use pmf_pmemd_new_control
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
    use pmf_pmemd_new_control
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
    call pmf_init_title('PMEMD-NEW')

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
    use pmf_pmemd_new_control

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

subroutine pmf_pmemd_bcast_dat_mpi()

    use pmf_init
    use pmf_pmemd_new_dat
    use pmf_dat
    use pmf_utils

    integer :: alloc_failed
    ! --------------------------------------------------------------------------

    call pmf_init_bcast_dat_mpi

    allocate(chunk_sizes(fnumoftasks), &
             chunk_offsets(fnumoftasks), &
             send_buffer(3,NumOfLAtoms), &
             recv_buffer(3,NumOfLAtoms), &
             stat=alloc_failed)

     if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,&
             '[PMEMD] Unable to allocate memory for chunck_sizes/chunck_offsets/send_buffer/recv_buffer arrays!')
     end if

    ! allocate temporary arrays
    if( fmaster ) then
        numofmpitransfers = 0
        fragmentation = 0
        numofchunkgaps = 0
        allocate(tmp_a(3,fnatoms), &
                 tmp_b(3,fnatoms), &
                 tmp_c(3,fnatoms), &
                 accu_chunk_sizes(fnumoftasks), &
                 stat=alloc_failed)

         if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[PMEMD] Unable to allocate memory for tmp_a/b/c arrays!')
         end if

        accu_chunk_sizes(:) = 0
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
! subroutine pmf_pmemd_force
!===============================================================================

subroutine pmf_pmemd_force(x,v,f,epot,ekin,epmf)

    use pmf_sizes
    use pmf_core_lf
    use pmf_dat
    use pmf_timers

    implicit none
    real(PMFDP)    :: x(:,:)            ! in
    real(PMFDP)    :: v(:,:)            ! in
    real(PMFDP)    :: f(:,:)            ! inout
    real(PMFDP)    :: epot              ! in
    real(PMFDP)    :: ekin              ! in
    real(PMFDP)    :: epmf              ! out
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_update_step
    call pmf_core_lf_force(x,v,f,epot,ekin,0.0d0,0.0d0,0.0d0,0.0d0,epmf)
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
! subroutine pmf_pmemd_force_mpi
!===============================================================================

subroutine pmf_pmemd_force_mpi(x,v,f,epot,ekin,epmf,atm_owner_map)

    use pmf_sizes
    use pmf_core_lf
    use pmf_pmemd_new_dat
    use pmf_dat
    use pmf_timers
    use pmf_utils

    implicit none

INCLUDE 'mpif.h'

    real(PMFDP)    :: x(:,:)            ! in
    real(PMFDP)    :: v(:,:)            ! in
    real(PMFDP)    :: f(:,:)            ! inout
    real(PMFDP)    :: epot              ! in
    real(PMFDP)    :: ekin              ! in
    real(PMFDP)    :: epmf              ! out
    integer        :: atm_owner_map(:)  ! in - atom map among processes
    ! ------------------------------------------------------
    integer        :: i, ierr
    ! --------------------------------------------------------------------------

    if(fmaster) then
        call pmf_timers_start_timer(PMFLIB_TIMER)
    end if

    ! gather data
    call pmf_pmemd_gather_array_mpi(tmp_a,x,atm_owner_map)
    call pmf_pmemd_gather_array_mpi(tmp_b,v,atm_owner_map)
    call pmf_pmemd_gather_array_mpi(tmp_c,f,atm_owner_map)

    if( fmaster ) then
        if( fdebug ) then
            write(PMF_DEBUG+fmytaskid,'(A)') '>> Input coordinates (force): '
            do i=1,NumOfLAtoms
                write(PMF_DEBUG+fmytaskid,'(3X,3F10.3)') tmp_a(:,RIndexes(i))
            end do
            write(PMF_DEBUG+fmytaskid,*)
        end if

        call pmf_core_lf_update_step
        call pmf_core_lf_force(tmp_a,tmp_b,tmp_c,epot,ekin,0.0d0,epmf)
    end if

    ! broadcast MD exit status
    call mpi_bcast(fexit_mdloop, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast fexit_mdloop variable!')
    end if

    ! update data
    call pmf_pmemd_scatter_array_mpi(tmp_c,f,atm_owner_map)

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
    use pmf_pmemd_new_dat
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
    call pmf_pmemd_gather_array_mpi(tmp_a,x,atm_owner_map)

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
    call pmf_pmemd_scatter_array_mpi(tmp_a,x,atm_owner_map)

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

#ifdef MPI
    call pmf_pmemd_mpistat
#endif

    call pmf_timers_stop_timer(PMFLIB_TOTAL_TIMER)
    call pmf_finalize_all(.true.)

end subroutine pmf_pmemd_finalize

#ifdef MPI

!===============================================================================
! Subroutine:  pmf_pmemd_mpistat
!===============================================================================

subroutine pmf_pmemd_mpistat

    use pmf_constants
    use pmf_dat
    use pmf_pmemd_new_dat

    implicit none
    integer           :: i
    real(PMFDP)       :: f,tot
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    write(PMF_OUT,'(A)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [MPI] <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

    if( numofmpitransfers .le. 0 ) return

    write(PMF_OUT,10) numofmpitransfers
    write(PMF_OUT,15) fnumoftasks
    f = real(fragmentation,PMFDP) / real(numofmpitransfers,PMFDP)
    write(PMF_OUT,20) f
    f = real(numofchunkgaps,PMFDP) / real(numofmpitransfers,PMFDP)
    write(PMF_OUT,23) f

    write(PMF_OUT,25)
    write(PMF_OUT,30)
    tot = 0.0
    do i=1,fnumoftasks
        f = real(accu_chunk_sizes(i),PMFDP) / real(3*numofmpitransfers,PMFDP)
        write(PMF_OUT,35) i, f
        tot = tot + f
    end do
    write(PMF_OUT,40)
    write(PMF_OUT,45) tot, NumOfLAtoms

 10 format('| MPI> Num of gather/scatter calls        = ',I12)
 15 format('| MPI> Num of CPUs                        = ',I12)
 20 format('| MPI> Average fragmentation per all CPUs = ',F17.4)
 23 format('| MPI> Average num of gaps between chunks = ',F17.4)
 25 format('| MPI> CPU <chunksize(atoms)>')
 30 format('| MPI> --- ------------------')
 35 format('| MPI> ',I3,1X,F18.4)
 40 format('| MPI> --- ------------------')
 45 format('| MPI> Total = ',F14.4,' (NumOfLAtoms =',I6,')')

end subroutine pmf_pmemd_mpistat

!===============================================================================
! Subroutine:  gather_array_mpi
!===============================================================================

subroutine pmf_pmemd_gather_array_mpi(lx,x,atm_owner_map)

    use pmf_dat
    use pmf_utils
    use pmf_pmemd_new_dat

    implicit none

INCLUDE 'mpif.h'

    real(PMFDP)       :: lx(:,:)  ! local
    real(PMFDP)       :: x(:,:)   ! scattered data
    integer           :: atm_owner_map(:)
    ! ------------------------------------------------------
    integer           :: i,at_id,ow_id,ierr,sat_id,slen
    logical           :: start
    real(PMFDP)       :: f
    ! --------------------------------------------------------------------------

    ierr = MPI_SUCCESS

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> gather_array'
    end if

    numofmpitransfers = numofmpitransfers + 1

    ! prepare input data
    sat_id = 1
    chunk_sizes(:) = 0
    do i=1,NumOfLAtoms
        at_id = RIndexes(i)
        ow_id = atm_owner_map(at_id)
        if( ow_id .eq. fmytaskid ) then
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,'(I8,A,I5,1X,3F10.3)') numofmpitransfers, '   > send atom: ',at_id,x(:,at_id)
            end if
            send_buffer(:,sat_id) = x(:,at_id)
            sat_id = sat_id + 1
        end if
        chunk_sizes(ow_id+1) = chunk_sizes(ow_id+1) + 3
    end do
    chunk_offsets(:) = 0
    slen = 0
    do i=2,fnumoftasks
        slen = slen + chunk_sizes(i-1)
        chunk_offsets(i) = slen
    end do

    ! MPI statistics
    if( fmaster ) then
        accu_chunk_sizes(:) = accu_chunk_sizes(:) + chunk_sizes(:)
        slen = 0
        start = .false.
        do i=1,fnumoftasks
            if( chunk_sizes(i) .ne. 0 ) then
                if( .not. start ) start = .true.
                fragmentation = fragmentation + 1
                slen = slen + chunk_sizes(i)
            end if
            if( slen .eq. 3*NumOfLAtoms ) start = .false.
            if( start .and. chunk_sizes(i) .eq. 0 ) then
                numofchunkgaps = numofchunkgaps + 1
            end if
        end do
    end if

    if( fdebug ) then
        do i=1,fnumoftasks
            write(PMF_DEBUG+fmytaskid,'(A,I3,A,I7,A,I7)') '    cpu id: ',i,', offset = ', &
                                                          chunk_offsets(i)/3,', size = ',chunk_sizes(i)/3
        end do
        write(PMF_DEBUG+fmytaskid,'(A,I3,A,I7)') '   task id: ', fmytaskid, &
                                                          ', number of atoms: ', chunk_sizes(fmytaskid+1)/3
    end if

    ! sanity check
    slen = 0
    do i=1,fnumoftasks
        if( chunk_sizes(i) .ne. 0 ) then
            slen = chunk_offsets(i) + chunk_sizes(i)
        end if
    end do
    if( slen .ne. 3*NumOfLAtoms ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Problem with data distribution in pmf_pmemd_gather_array_mpi!')
    end if

    ! write info about fragmentation
    if( fmaster .and. frepmpifrag ) then
        f = real(fragmentation,PMFDP) / real(numofmpitransfers,PMFDP)
        if( f .gt. 1.0d0 ) then
            write(PMF_DEBUG-1,20) numofmpitransfers, f
        end if
    end if

    ! gather data
    call mpi_gatherv(send_buffer,chunk_sizes(fmytaskid+1),mpi_double_precision, &
                     recv_buffer,chunk_sizes(1),chunk_offsets(1),mpi_double_precision, &
                     0,mpi_comm_world,ierr)

    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to mpi_gatherv!')
    end if

    ! recollect data - only master
    if( fmaster ) then
        chunk_sizes(:) = 0
        do i=1,NumOfLAtoms
            at_id = RIndexes(i)
            ow_id = atm_owner_map(at_id) + 1
            chunk_sizes(ow_id) = chunk_sizes(ow_id) + 1
            lx(:,at_id) = recv_buffer(:,chunk_sizes(ow_id)+chunk_offsets(ow_id)/3)
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,'(I8,A,I5,1X,3F10.3)') numofmpitransfers,'   < recv atom: ',at_id,lx(:,at_id)
            end if
        end do
    end if

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*)
    end if

 20 format('| MPI> ',I9,'Average fragmentation per all CPUs = ',F17.4)

end subroutine pmf_pmemd_gather_array_mpi

!===============================================================================
! Subroutine:  scatter_array_mpi
!===============================================================================

subroutine pmf_pmemd_scatter_array_mpi(lx,x,atm_owner_map)

    use pmf_dat
    use pmf_utils
    use pmf_pmemd_new_dat

    implicit none

INCLUDE 'mpif.h'

    real(PMFDP)       :: lx(:,:)
    real(PMFDP)       :: x(:,:)
    integer           :: atm_owner_map(:)
    ! ------------------------------------------------------
    integer           :: i,at_id,ow_id,ierr,sat_id,slen
    logical           :: start
    ! --------------------------------------------------------------------------

    ierr = MPI_SUCCESS

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,'(A)') '>> scatter_array'
    end if

    numofmpitransfers = numofmpitransfers + 1

    ! prepare input data
    sat_id = 1
    chunk_sizes(:) = 0
    do i=1,NumOfLAtoms
        at_id = RIndexes(i)
        ow_id = atm_owner_map(at_id) + 1
        chunk_sizes(ow_id) = chunk_sizes(ow_id) + 3
    end do
    chunk_offsets(:) = 0
    slen = 0
    do i=2,fnumoftasks
        slen = slen + chunk_sizes(i-1)
        chunk_offsets(i) = slen
    end do

    if( fmaster ) then
        chunk_sizes(:) = 0
        send_buffer(:,:) = 0
        do i=1,NumOfLAtoms
            at_id = RIndexes(i)
            ow_id = atm_owner_map(at_id) + 1
            chunk_sizes(ow_id) = chunk_sizes(ow_id) + 3
            send_buffer(:,(chunk_sizes(ow_id)+chunk_offsets(ow_id))/3) = lx(:,at_id)
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,'(I8,A,I5,1X,3F10.3)') numofmpitransfers, '   > send atom: ',at_id,lx(:,at_id)
            end if
        end do
        ! MPI stat
        accu_chunk_sizes(:) = accu_chunk_sizes(:) + chunk_sizes(:)
        slen = 0
        start = .false.
        do i=1,fnumoftasks
            if( chunk_sizes(i) .ne. 0 ) then
                if( .not. start ) start = .true.
                fragmentation = fragmentation + 1
                slen = slen + chunk_sizes(i)
            end if
            if( slen .eq. 3*NumOfLAtoms ) start = .false.
            if( start .and. chunk_sizes(i) .eq. 0 ) then
                numofchunkgaps = numofchunkgaps + 1
            end if
        end do
    end if

    if( fdebug ) then
        do i=1,fnumoftasks
            write(PMF_DEBUG+fmytaskid,'(A,I3,A,I7,A,I7)') '    cpu id: ',i,', offset = ', &
                                                          chunk_offsets(i)/3,', size = ',chunk_sizes(i)/3
        end do
        write(PMF_DEBUG+fmytaskid,'(A,I3,A,I7)') '   task id: ', fmytaskid, ', number of atoms: ', chunk_sizes(fmytaskid+1)/3
    end if

    ! sanity check
    slen = 0
    do i=1,fnumoftasks
        if( chunk_sizes(i) .ne. 0 ) then
            slen = chunk_offsets(i) + chunk_sizes(i)
        end if
    end do
    if( slen .ne. 3*NumOfLAtoms ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Problem with data distribution in pmf_pmemd_scatter_array_mpi!')
    end if

    ! scatter data
    call mpi_scatterv(send_buffer,chunk_sizes(1),chunk_offsets(1),mpi_double_precision, &
                     recv_buffer,chunk_sizes(fmytaskid+1),mpi_double_precision, &
                     0,mpi_comm_world,ierr)

    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to mpi_scatterv!')
    end if

    ! restore data
    sat_id = 1
    do i=1,NumOfLAtoms
        at_id = RIndexes(i)
        ow_id = atm_owner_map(at_id)
        if( ow_id .eq. fmytaskid ) then
            x(:,at_id) = recv_buffer(:,sat_id)
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,'(I8,A,I5,1X,3F10.3)') numofmpitransfers,'   < recv atom: ',at_id,x(:,at_id)
            end if
            sat_id = sat_id + 1
        end if
    end do

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*)
    end if

end subroutine pmf_pmemd_scatter_array_mpi
!===============================================================================

#endif

!===============================================================================

end module pmf_pmemd_new




