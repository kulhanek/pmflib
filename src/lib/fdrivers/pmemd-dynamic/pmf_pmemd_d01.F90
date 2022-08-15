!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_pmemd_d01

use pmf_sizes
use iso_c_binding
use pmf_pmemd_dat_d01

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_check_interface(rnum,inum,ekin_len,setup_len,str1,str1_len,str2,str2_len) &
           bind(c,name='int_pmf_pmemd_check_interface')

    use pmf_constants
    use pmf_utils

    implicit none
    real(CPMFDP)        :: rnum
    integer(CPMFINT)    :: inum
    integer(CPMFINT)    :: ekin_len
    integer(CPMFINT)    :: setup_len
    character(CPMFCHAR) :: str1(*)
    integer(CPMFINT)    :: str1_len
    character(CPMFCHAR) :: str2(*)
    integer(CPMFINT)    :: str2_len
    ! -------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    if( rnum .ne. PMFLIB_CHECK_R81 ) then
        call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for PMFLIB_CHECK_R81!')
    end if
    if( inum .ne. PMFLIB_CHECK_INT1 ) then
        call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for PMFLIB_CHECK_INT1!')
    end if

    if( ekin_len .ne. PMFLIB_EKIN_SIZE ) then
        call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for PMFLIB_EKIN_SIZE!')
    end if
    if( setup_len .ne. PMFLIB_SETUP_SIZE ) then
        call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for PMFLIB_SETUP_SIZE!')
    end if

    if( str1_len .ne. len(PMFLIB_CHECK_STR1) ) then
        call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for len(PMFLIB_CHECK_STR1)!')
    end if
    do i=1,min(str1_len,len(PMFLIB_CHECK_STR1))
        if( str1(i) .ne. PMFLIB_CHECK_STR1(i:i) ) then
            call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for PMFLIB_CHECK_STR1!')
        end if
    end do

    if( str2_len .ne. len(PMFLIB_CHECK_STR2) ) then
        call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for len(PMFLIB_CHECK_STR2)!')
    end if
    do i=1,min(str2_len,len(PMFLIB_CHECK_STR2))
        if( str2(i) .ne. PMFLIB_CHECK_STR2(i:i) ) then
            call pmf_utils_exit(PMF_OUT,1,'Driver interface compromised for PMFLIB_CHECK_STR2!')
        end if
    end do

end subroutine pmf_pmemd_check_interface

!===============================================================================
! subroutine pmf_pmemd_init_preinit
!===============================================================================

subroutine pmf_pmemd_init_preinit(mdin,mdin_len,anatom,anres,   &
                            antb,antc,ansteps,astepsize,atemp0, &
                            box_a,box_b,box_c,                  &
                            box_alpha,box_beta,box_gamma) bind(c,name='int_pmf_pmemd_init_preinit')

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_utils
    use pmf_pmemd_control_d01
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers

    implicit none
    character(CPMFCHAR) :: mdin(*)
    integer(CPMFINT)    :: mdin_len
    integer(CPMFINT)    :: anatom                       ! number of atoms in AMBER topology
    integer(CPMFINT)    :: anres                        ! number of residues in AMBER topology
    integer(CPMFINT)    :: antb                         ! BOX type
    integer(CPMFINT)    :: antc                         ! shake mode
    integer(CPMFINT)    :: ansteps                      ! number of MD steps
    real(CPMFDP)        :: astepsize                    ! step size
    real(CPMFDP)        :: atemp0                       ! temperature
    real(CPMFDP)        :: apress0                      ! pressure
    real(CPMFDP)        :: box_a,box_b,box_c            ! box dimensions
    real(CPMFDP)        :: box_alpha,box_beta,box_gamma
    ! ------------------------------------------------------
    integer             :: has_box, i
    real(PMFDP)         :: cbox(3)
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    ! setup conversion factors
    MassConv         = 1.0d0                ! g/mol -> g/mol
    LengthConv       = 1.0d0                ! A -> A
    AngleConv        = PMF_D2R              ! deg -> rad
    TimeConv         = 1000.0d0             ! ps -> fs
    VelocityConv     = 1.0d0                ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0                ! kcal/mol -> kcal/mol
    ForceConv        = 1.0d0                ! kcal/mol/A -> kcal/mol/A
    TemperatureConv  = 1.0d0                ! K
    PressureConv     = 1.0d5                ! Bar -> Pa

    ControlFileName = ''
    do i=1,min(mdin_len,len(ControlFileName))
        ControlFileName(i:i)  = mdin(i)
    end do

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_LEAP_FROG,anatom,antb,ansteps,astepsize,0.0d0,atemp0,apress0)
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

end subroutine pmf_pmemd_init_preinit

!===============================================================================
! subroutine pmf_pmemd_set_residue
!===============================================================================

subroutine pmf_pmemd_set_residue(idx,name,name_len,first_atom) bind(c,name='int_pmf_pmemd_set_residue')

    use pmf_mask
    use pmf_dat

    implicit none
    integer(CPMFINT)    :: idx
    character(CPMFCHAR) :: name(*)
    integer(CPMFINT)    :: name_len
    integer(CPMFINT)    :: first_atom
    ! -------------------------------------------
    character(len=PMF_MAX_MASK_NAME)    :: lname
    integer                             :: i
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    lname=''
    do i=1,min(name_len,PMF_MAX_MASK_NAME)
        lname(i:i) = name(i)
    end do

    call pmf_mask_set_topo_residue(idx,lname,first_atom)

end subroutine pmf_pmemd_set_residue

!===============================================================================
! subroutine pmf_pmemd_set_atom
!===============================================================================

subroutine pmf_pmemd_set_atom(idx,name,name_len,atype,atype_len) bind(c,name='int_pmf_pmemd_set_atom')

    use pmf_mask
    use pmf_dat

    implicit none
    integer(CPMFINT)    :: idx
    character(CPMFCHAR) :: name(*)
    integer(CPMFINT)    :: name_len
    character(CPMFCHAR) :: atype(*)
    integer(CPMFINT)    :: atype_len
    ! -------------------------------------------
    character(len=PMF_MAX_MASK_NAME)    :: lname
    character(len=PMF_MAX_MASK_NAME)    :: latype
    integer                             :: i
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    lname=''
    do i=1,min(name_len,PMF_MAX_MASK_NAME)
        lname(i:i) = name(i)
    end do

    latype=''
    do i=1,min(atype_len,PMF_MAX_MASK_NAME)
        latype(i:i) = atype(i)
    end do

    call pmf_mask_set_topo_atom(idx,lname,latype)

end subroutine pmf_pmemd_set_atom

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_finalize_preinit(amass,ax) bind(c,name='int_pmf_pmemd_finalize_preinit')

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_mask
    use pmf_core
    use pmf_init
    use pmf_pmemd_control_d01
    use pmf_utils

    implicit none
    real(CPMFDP)   :: amass(:)
    real(CPMFDP)   :: ax(:,:)
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

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_init(amass,ax) bind(c,name='int_pmf_pmemd_init')

    use pmf_constants
    use pmf_sizes
    use pmf_init
    use pmf_dat
    use cst_init
    use pmf_pmemd_control_d01

    implicit none
    real(CPMFDP)   :: amass(:)
    real(CPMFDP)   :: ax(:,:)
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

subroutine pmf_pmemd_get_setup(setup) bind(c,name='int_pmf_pmemd_get_setup')

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use abf_dat
    use pmf_cvs
    use pmf_utils

    implicit none
#ifdef MPI
INCLUDE 'mpif.h'
#endif

    integer(CPMFINT)    :: setup(:)
    ! -------------------------------------------
    integer             :: i
#ifdef MPI
    integer             :: ierr
#endif
    ! --------------------------------------------------------------------------

    ! reset setup
    setup(:) = 0

    if( fmaster ) then
        if( cst_enabled ) then
            call pmf_utils_exit(PMF_OUT,1,'[cst] method is not supported with this PMFLib driver!')
        end if

        ! [EPOT] requires ENE and FRC
        do i=1,NumOfCVs
            if( CVList(i)%cv%ctype .eq. 'EPOT' ) then
                setup(PMFLIB_SETUP_FORCE_NEED_ENE) = 1
                setup(PMFLIB_SETUP_FORCE_NEED_FRC) = 1
                exit
            end if
        end do

        if( fmode .eq. 2 ) then
            ! we need velocities for this ABF mode
            setup(PMFLIB_SETUP_FORCE_NEED_VEL) = 1
        end if
    end if

#ifdef MPI
    ! broadcast to other processes
    call mpi_bcast(setup, size(setup), mpi_integer, 0, mpi_comm_world, ierr)
    if( ierr .ne. MPI_SUCCESS ) then
        call pmf_utils_exit(PMF_OUT, 1,'[PMF] Unable to broadcast setup array in pmf_pmemd_get_setup!')
    end if
#endif

    return

end subroutine pmf_pmemd_get_setup

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_pmemd_shouldexit(exitcode) bind(c,name='int_pmf_pmemd_shouldexit')

    use pmf_constants
    use pmf_sizes
    use pmf_dat

    implicit none
    integer(CPMFINT)    :: exitcode       ! MD loop exit code
    ! --------------------------------------------------------------------------

    ! we cannot bcast the value here as pmf_pmemd_shouldexit is called on master
    ! protected sections several times per MD cycle
    ! thus fexit_mdloop is distributed in pmf_pmemd_force_mpi

    exitcode = fexit_mdloop

    return

end subroutine pmf_pmemd_shouldexit

!===============================================================================
! subroutine pmf_pmemd_finalize
!===============================================================================

subroutine pmf_pmemd_finalize() bind(c,name='int_pmf_pmemd_finalize')

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

!===============================================================================
! subroutine pmf_pmemd_update_box
!===============================================================================

subroutine pmf_pmemd_update_box(a,b,c,alpha,beta,gamma) bind(c,name='int_pmf_pmemd_update_box')

    use pmf_dat
    use pmf_pbc

    real(CPMFDP)    :: a,b,c
    real(CPMFDP)    :: alpha,beta,gamma
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fsystype .eq. SYS_NT ) return  ! no box -> return

    call pmf_pbc_set_box(a,b,c,alpha,beta,gamma)

end subroutine pmf_pmemd_update_box

!===============================================================================
! subroutine pmf_pmemd_force
!===============================================================================

subroutine pmf_pmemd_force(x,v,f,epot,epmf) bind(c,name='int_pmf_pmemd_force')

    use pmf_sizes
    use pmf_core_lf
    use pmf_dat
    use pmf_timers

    implicit none
    real(CPMFDP)    :: x(:,:)       ! in
    real(CPMFDP)    :: v(:,:)       ! in
    real(CPMFDP)    :: f(:,:)       ! inout
    real(CPMFDP)    :: epot         ! in
    real(CPMFDP)    :: epmf         ! out
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_update_step
    call pmf_core_lf_force(x,v,f,epot,epmf)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_pmemd_force

!===============================================================================
! subroutine pmf_pmemd_register_ekin
!===============================================================================

subroutine pmf_pmemd_register_ekin(ekin,valid) bind(c,name='int_pmf_pmemd_register_ekin')

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers
    use pmf_dat
    use pmf_utils

    implicit none
    real(CPMFDP)            :: ekin(:)      ! in
    integer(CPMFINT)        :: valid
    ! --------------------------------------------
    type(PMFKineticEnergy)  :: sekin
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_timers_start_timer(PMFLIB_TIMER)

    sekin%KinEneVV = ekin(PMFLIB_EKIN_VV)
    sekin%KinEneLF = ekin(PMFLIB_EKIN_LF)
    sekin%KinEneHA = ekin(PMFLIB_EKIN_HA)
    sekin%Valid = .false.
    if( valid .eq. 1 ) sekin%Valid = .true.

    call pmf_core_lf_register_ekin(sekin)

    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_pmemd_register_ekin

#ifdef MPI

!===============================================================================
! subroutine pmf_pmemd_init_taskid_mpi
!===============================================================================

subroutine pmf_pmemd_init_taskid_mpi(mytaskid,numoftasks) bind(c,name='int_pmf_pmemd_init_taskid_mpi')

    use pmf_init

    implicit none
    integer(CPMFINT)    :: mytaskid
    integer(CPMFINT)    :: numoftasks
    ! --------------------------------------------------------------------------

    call pmf_init_taskid_mpi(mytaskid,numoftasks)

end subroutine pmf_pmemd_init_taskid_mpi

!===============================================================================
! subroutine pmf_pmemd_bcast_dat_mpi
!===============================================================================

subroutine pmf_pmemd_bcast_dat_mpi() bind(c,name='int_pmf_pmemd_bcast_dat_mpi')

    use pmf_init
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
! subroutine pmf_pmemd_force_mpi
!===============================================================================

subroutine pmf_pmemd_force_mpi(x,v,f,epot,epmf,atm_owner_map) bind(c,name='int_pmf_pmemd_force_mpi')

    use pmf_sizes
    use pmf_core_lf
    use pmf_dat
    use pmf_timers
    use pmf_utils

    implicit none

INCLUDE 'mpif.h'

    real(CPMFDP)        :: x(:,:)           ! in
    real(CPMFDP)        :: v(:,:)           ! in
    real(CPMFDP)        :: f(:,:)           ! inout
    real(CPMFDP)        :: epot             ! in
    real(CPMFDP)        :: epmf             ! out
    integer(CPMFINT)    :: atm_owner_map(:) ! in - atom map among processes
    ! ------------------------------------------------------
    integer             :: i, ierr
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
        call pmf_core_lf_force(tmp_a,tmp_b,tmp_c,epot,epmf)
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
! Subroutine:  pmf_pmemd_mpistat
!===============================================================================

subroutine pmf_pmemd_mpistat

    use pmf_constants
    use pmf_dat

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

end module pmf_pmemd_d01




