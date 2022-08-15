!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022 Petr Kulhanek, kulhanek@chemi.muni.cz
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

! ==============================================================================
! this a client side of the SANDER driver
! ==============================================================================

module PMFLibSANDER_d01

use iso_c_binding

implicit none

#if defined PMFLIB

! ==============================================================================
! type sizes
integer, parameter                          :: CPMFDP        = c_double
integer, parameter                          :: CPMFINT       = c_int
integer, parameter                          :: CPMFCHAR      = c_char
! interface binding check
integer(CPMFINT), parameter                 :: PMFLIB_CHECK_INT1 = 1089523658
real(CPMFDP), parameter                     :: PMFLIB_CHECK_R81  = 1.78493547
character(kind=CPMFCHAR,len=10), parameter  :: PMFLIB_CHECK_STR1 = 'PMFLib v06'
character(kind=CPMFCHAR,len=10), parameter  :: PMFLIB_CHECK_STR2 = 'DRVABI d01'
! ==============================================================================

! constants
integer(c_int), parameter   :: RTLD_NOW                     = 2         ! extracted from libc
integer(c_int), parameter   :: RTLD_DI_ORIGIN               = 6         ! extracted from libc
integer, parameter          :: MAX_PATH                     = 4096      ! Linux

! energy array
integer, parameter          :: PMFLIB_EKIN_VV               = 1
integer, parameter          :: PMFLIB_EKIN_LF               = 2
integer, parameter          :: PMFLIB_EKIN_HA               = 3
integer, parameter          :: PMFLIB_EKIN_SIZE             = PMFLIB_EKIN_HA

! setup array
integer, parameter          :: PMFLIB_SETUP_ISCHEME         = 1
integer, parameter          :: PMFLIB_SETUP_SIZE            = PMFLIB_SETUP_ISCHEME

! ==============================================================================
! interface to linux API
interface
    function dlopen(filename,mode) bind(c,name="dlopen")
        ! void *dlopen(const char *filename, int mode);
        use iso_c_binding
        implicit none
        type(c_ptr)                     :: dlopen
        character(c_char), intent(in)   :: filename(*)
        integer(c_int), value           :: mode
    end function
    ! -------------------------------------------------------------------------
    function dlinfo(handle,request,info) bind(c,name="dlinfo")
        ! int dlinfo(void *handle, int request, void *info);
        use iso_c_binding
        implicit none
        integer(c_int)                  :: dlinfo
        type(c_ptr), value              :: handle
        integer(c_int), value           :: request
        character(c_char)               :: info(*)
    end function
    ! -------------------------------------------------------------------------
    function dlerror() bind(c,name="dlerror")
        ! char *dlerror(void);
        use iso_c_binding
        implicit none
        type(C_PTR)                     :: dlerror
    end function
    ! -------------------------------------------------------------------------
    function dlsym(handle,name) bind(c,name="dlsym")
        ! void *dlsym(void *handle, const char *name);
        use iso_c_binding
        implicit none
        type(c_funptr)                  :: dlsym
        type(c_ptr), value              :: handle
        character(c_char), intent(in)   :: name(*)
    end function
    ! -------------------------------------------------------------------------
    function dlclose(handle) bind(c,name="dlclose")
        ! int dlclose(void *handle);
        use iso_c_binding
        implicit none
        integer(c_int)      :: dlclose
        type(c_ptr), value  :: handle
    end function
    ! -------------------------------------------------------------------------
    function strlen(string) bind(c,name='strlen')
        ! int strlen(char *string)
        use iso_c_binding
        implicit none
        integer(c_int)      :: strlen
        type(c_ptr), value  :: string
    end function
end interface

! ==============================================================================
! interface to PMFLib driver
abstract interface

! SETUP ==========================================
    subroutine int_pmf_sander_check_interface(rnum,inum,ene_len,setup_len,str1,str1_len,str2,str2_len) bind(c)
        import
        implicit none
        real(CPMFDP)        :: rnum
        integer(CPMFINT)    :: inum
        integer(CPMFINT)    :: ene_len
        integer(CPMFINT)    :: setup_len
        character(CPMFCHAR) :: str1(*)
        integer(CPMFINT)    :: str1_len
        character(CPMFCHAR) :: str2(*)
        integer(CPMFINT)    :: str2_len
    end subroutine int_pmf_sander_check_interface
    ! -------------------------------------------------------------------------
    subroutine int_pmf_sander_init_preinit(mdin,mdin_len,anatom,anres,  &
                            antb,antc,ansteps,astepsize,atemp0,apress0, &
                            box_a,box_b,box_c,                          &
                            box_alpha,box_beta,box_gamma) bind(c)
        import
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
    end subroutine int_pmf_sander_init_preinit
    ! -------------------------------------------------------------------------
    subroutine int_pmf_sander_set_residue(idx,name,name_len,first_atom) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: idx
        character(CPMFCHAR) :: name(*)
        integer(CPMFINT)    :: name_len
        integer(CPMFINT)    :: first_atom
    end subroutine int_pmf_sander_set_residue
    ! -------------------------------------------------------------------------
    subroutine int_pmf_sander_set_atom(idx,name,name_len,atype,atype_len) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: idx
        character(CPMFCHAR) :: name(*)
        integer(CPMFINT)    :: name_len
        character(CPMFCHAR) :: atype(*)
        integer(CPMFINT)    :: atype_len
    end subroutine int_pmf_sander_set_atom
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_finalize_preinit(anatom,amass,ax) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        real(CPMFDP)        :: amass(anatom)
        real(CPMFDP)        :: ax(3,anatom)
    end subroutine int_pmf_sander_finalize_preinit
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_init(anatom,amass,ax) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        real(CPMFDP)        :: amass(anatom)
        real(CPMFDP)        :: ax(3,anatom)
    end subroutine int_pmf_sander_init
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_get_setup(setup,setup_len) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: setup(:)
        integer(CPMFINT)    :: setup_len
    end subroutine int_pmf_sander_get_setup
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_shouldexit(exitcode) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: exitcode
    end subroutine int_pmf_sander_shouldexit
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_finalize() bind(c)
        implicit none
    end subroutine int_pmf_sander_finalize

! FORCES =========================================
    subroutine int_pmf_sander_update_box(a,b,c,alpha,beta,gamma) bind(c)
        import
        implicit none
        real(CPMFDP)        :: a,b,c
        real(CPMFDP)        :: alpha,beta,gamma
    end subroutine int_pmf_sander_update_box
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_force(anatom,x,v,f,epot,epmf) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        real(CPMFDP)        :: x(*)         ! in
        real(CPMFDP)        :: v(*)         ! in
        real(CPMFDP)        :: f(*)         ! inout
        real(CPMFDP)        :: epot         ! in
        real(CPMFDP)        :: epmf         ! out
    end subroutine int_pmf_sander_force
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_register_ekin(ekin) bind(c)
        import
        implicit none
        real(CPMFDP)        :: ekin(:)      ! in
    end subroutine int_pmf_sander_register_ekin

! CST ============================================
    subroutine int_pmf_sander_cst_init_collisions(ntc,nbt,ifstwt,ib,jb,conp) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: ntc
        integer(CPMFINT)    :: nbt
        integer(CPMFINT)    :: ifstwt(*)
        integer(CPMFINT)    :: ib(*)
        integer(CPMFINT)    :: jb(*)
        real(CPMFDP)        :: conp(*)
    end subroutine int_pmf_sander_cst_init_collisions
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_num_of_pmflib_cst(numofcst) bind(c)
        import
        implicit none
        real(CPMFDP)    :: numofcst
    end subroutine int_pmf_sander_num_of_pmflib_cst
    ! --------------------------------------------------------------------------
    function int_pmf_sander_cst_checkatom(atomid) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: int_pmf_sander_cst_checkatom
        integer(CPMFINT)    :: atomid
    end function int_pmf_sander_cst_checkatom
    ! --------------------------------------------------------------------------
    subroutine int_pmf_sander_shake(anatom,x,modified) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        real(CPMFDP)        :: x(*)
        integer(CPMFINT)    :: modified
    end subroutine int_pmf_sander_shake

#ifdef MPI
! MPI ============================================
    subroutine int_pmf_sander_init_taskid_mpi(mytaskid,numoftasks) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: mytaskid
        integer(CPMFINT)    :: numoftasks
    end subroutine int_pmf_sander_init_taskid_mpi
    ! -------------------------------------------------------------------------
    subroutine int_pmf_sander_bcast_dat_mpi(anatom,anumtasks,aiparpt) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        integer(CPMFINT)    :: anumtasks
        integer(CPMFINT)    :: aiparpt(:)
    end subroutine int_pmf_sander_bcast_dat_mpi
    ! -------------------------------------------------------------------------
    subroutine int_pmf_sander_force_mpi(anatom,x,v,f,epot,epmf) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        real(CPMFDP)        :: x(*)
        real(CPMFDP)        :: v(*)
        real(CPMFDP)        :: f(*)
        real(CPMFDP)        :: epot         ! in
        real(CPMFDP)        :: epmf         ! out
    end subroutine int_pmf_sander_force_mpi
    ! -------------------------------------------------------------------------
    subroutine int_pmf_sander_shake_mpi(anatom,x,modified) bind(c)
        import
        implicit none
        integer(CPMFINT)    :: anatom
        real(CPMFDP)        :: x(*)
        integer(CPMFINT)    :: modified
    end subroutine int_pmf_sander_shake_mpi
    ! -------------------------------------------------------------------------
#endif

end interface

! ==============================================================================
! driver setup
character(len=100)          :: pmf_sander_driver_name    = 'libfpmfdrv_sander_d01.so'
character(len=MAX_PATH)     :: pmf_sander_driver_path    = ''
character(len=MAX_PATH)     :: pmf_sander_driver_error   = ''
type(c_ptr)                 :: pmf_sander_driver_handle  = C_NULL_PTR

! driver symbol
procedure(int_pmf_sander_check_interface), bind(c), pointer         :: pmf_sander_check_interface

procedure(int_pmf_sander_init_preinit), bind(c), pointer            :: pmf_sander_init_preinit
procedure(int_pmf_sander_set_residue), bind(c), pointer             :: pmf_sander_set_residue
procedure(int_pmf_sander_set_atom), bind(c), pointer                :: pmf_sander_set_atom
procedure(int_pmf_sander_finalize_preinit), bind(c), pointer        :: pmf_sander_finalize_preinit
procedure(int_pmf_sander_init), bind(c), pointer                    :: pmf_sander_init
procedure(int_pmf_sander_get_setup), bind(c), pointer               :: pmf_sander_get_setup

procedure(int_pmf_sander_shouldexit), bind(c), pointer              :: pmf_sander_shouldexit
procedure(int_pmf_sander_finalize), bind(c), pointer                :: pmf_sander_finalize

procedure(int_pmf_sander_update_box), bind(c), pointer              :: pmf_sander_update_box
procedure(int_pmf_sander_force), bind(c), pointer                   :: pmf_sander_force
procedure(int_pmf_sander_register_ekin), bind(c), pointer           :: pmf_sander_register_ekin

procedure(int_pmf_sander_cst_init_collisions), bind(c), pointer     :: pmf_sander_cst_init_collisions
procedure(int_pmf_sander_num_of_pmflib_cst), bind(c), pointer       :: pmf_sander_num_of_pmflib_cst
procedure(int_pmf_sander_cst_checkatom), bind(c), pointer           :: pmf_sander_cst_checkatom
procedure(int_pmf_sander_shake), bind(c), pointer                   :: pmf_sander_shake

#ifdef MPI
procedure(int_pmf_sander_init_taskid_mpi), bind(c), pointer         :: pmf_sander_init_taskid_mpi
procedure(int_pmf_sander_bcast_dat_mpi), bind(c), pointer           :: pmf_sander_bcast_dat_mpi
procedure(int_pmf_sander_force_mpi), bind(c), pointer               :: pmf_sander_force_mpi
procedure(int_pmf_sander_shake_mpi), bind(c), pointer               :: pmf_sander_shake_mpi
#endif

! ==============================================================================
! run-time variables
logical             :: use_pmflib               = .false.
logical             :: pmflib_stop_on_failure   = .true.
integer             :: pmflib_cst_modified      = 0
integer             :: pmflib_exit              = 0
real(CPMFDP)        :: pmflib_ncst              = 0

real(CPMFDP)        :: pmflib_epot              = 0.0d0
real(CPMFDP)        :: pmflib_erst              = 0.0d0

real(CPMFDP)        :: pmflib_ekin(PMFLIB_EKIN_SIZE)
integer(CPMFINT)    :: pmflib_setup(PMFLIB_SETUP_SIZE)

! ==============================================================================

contains

!===============================================================================
! subroutine pmf_sander_bind_to_driver
!===============================================================================

subroutine pmf_sander_bind_to_driver(master)

    implicit none
    logical         :: master
    ! -------------------------------------------
    type(c_funptr)  :: proc_addr
    ! --------------------------------------------------------------------------

    use_pmflib = .false.

! print header
    if( master ) then
        write(6,*)
        write(6,10)
        write(6,20)
        write(6,25)
        write(6,30) trim(pmf_sander_driver_name)
    end if

! load the driver
    pmf_sander_driver_handle = dlopen(trim(pmf_sander_driver_name)//c_null_char, RTLD_NOW)
    if( .not. c_associated(pmf_sander_driver_handle) ) then
        if( master ) then
            call pmf_sander_get_dlerror
            write(6,40) trim(pmf_sander_driver_error)
            write(6,10)
            write(6,*)
        end if
        if( pmflib_stop_on_failure ) then
            stop 'Unable to load the PMFLib driver!'
        end if
        return
    end if

    pmf_sander_driver_path(:) = ' '
    if( master ) then
        if( dlinfo(pmf_sander_driver_handle,RTLD_DI_ORIGIN,pmf_sander_driver_path) .ne. 0 ) then
            pmf_sander_driver_path = '-not specified-'
        end if

        write(6,60) trim(pmf_sander_driver_path(1:max(0,index(pmf_sander_driver_path,c_null_char)-1)))
        write(6,50)
        write(6,70)
    end if

! bind procedures
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_check_interface"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_check_interface'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_check_interface)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_init_preinit"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_init_preinit'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_init_preinit)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_set_residue"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_set_residue'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_set_residue)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_set_atom"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_set_atom'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_set_atom)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_finalize_preinit"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure pmf_sander_finalize_preinit'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_finalize_preinit)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_init"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_init'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_init)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_get_setup"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_get_setup'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_get_setup)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_shouldexit"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_shouldexit'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_shouldexit)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_finalize"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_finalize'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_finalize)

! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_update_box"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_update_box'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_update_box)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_force"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_force'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_force)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_register_ekin"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_register_ekin'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_register_ekin)

! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_cst_init_collisions"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_cst_init_collisions'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_cst_init_collisions)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_num_of_pmflib_cst"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_num_of_pmflib_cst'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_num_of_pmflib_cst)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_cst_checkatom"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_cst_checkatom'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_cst_checkatom)
        ! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_shake"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_shake'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_shake)

#ifdef MPI
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_init_taskid_mpi"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_init_taskid_mpi'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_init_taskid_mpi)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_bcast_dat_mpi"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_bcast_dat_mpi'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_bcast_dat_mpi)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_force_mpi"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_force_mpi'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_force_mpi)
! ------------------
    proc_addr=dlsym(pmf_sander_driver_handle, "int_pmf_sander_shake_mpi"//c_null_char)
    if (.not. c_associated(proc_addr))then
        stop 'Unable to load the procedure int_pmf_sander_shake_mpi'
    end if
    call c_f_procpointer(proc_addr,pmf_sander_shake_mpi)
! ------------------
#endif

! print footer
    if( master ) then
        write(6,80)
        write(6,90)
        call pmf_sander_check_interface(PMFLIB_CHECK_R81,PMFLIB_CHECK_INT1,PMFLIB_EKIN_SIZE,PMFLIB_SETUP_SIZE, &
                                       PMFLIB_CHECK_STR1,len(PMFLIB_CHECK_STR1), &
                                       PMFLIB_CHECK_STR2,len(PMFLIB_CHECK_STR2))
        write(6,100)
        write(6,10)
        write(6,*)
    else
        call pmf_sander_check_interface(PMFLIB_CHECK_R81,PMFLIB_CHECK_INT1,PMFLIB_EKIN_SIZE,PMFLIB_SETUP_SIZE, &
                                       PMFLIB_CHECK_STR1,len(PMFLIB_CHECK_STR1), &
                                       PMFLIB_CHECK_STR2,len(PMFLIB_CHECK_STR2))
    end if

    use_pmflib = .true.
    return

 10 format('# -----------------------------------------------------------------------------')
 20 format('# PMFLib dynamic binding')
 25 format("# > Loading driver ...")
 30 format('#   Driver: ', A)
 60 format('#   Path:   ', A)
 40 format("#   Driver not loaded - disabling PMFLib functionality (",A,")!")
 50 format("#   Driver loaded successfully!")

 70 format("# > Loading symbols ...")
 80 format("#   All symbols loaded successfully!")

 90 format("# > Checking interface integrity ...")
100 format("#   Everything seems to be OK!")

end subroutine pmf_sander_bind_to_driver

!===============================================================================
! subroutine pmf_sander_update_setup
!===============================================================================

subroutine pmf_sander_update_setup(master,ischeme)

    implicit none
    logical         :: master
    integer         :: ischeme
    ! --------------------------------------------------------------------------

    if( .not. master ) return

    ! populate setup
    pmflib_setup(PMFLIB_SETUP_ISCHEME) = ischeme

    ! send the setup
    call pmf_sander_get_setup(pmflib_setup,PMFLIB_SETUP_SIZE)

end subroutine pmf_sander_update_setup

!===============================================================================
! subroutine pmf_sander_release_driver
!===============================================================================

subroutine pmf_sander_release_driver(master)

    implicit none
    logical         :: master
    ! --------------------------------------------------------------------------

    if( .not. use_pmflib ) return

    if( master ) then
        write(6,*)
        write(6,10)
        write(6,20)
        write(6,30)
    end if

    if( dlclose(pmf_sander_driver_handle) .ne. 0 ) then
        if( master ) then
            call pmf_sander_get_dlerror
            write(6,40) trim(pmf_sander_driver_error)
            write(6,10)
            write(6,*)
        end if
        return
    else
        if( master ) then
            write(6,50)
            write(6,10)
            write(6,*)
        end if
    end if

    use_pmflib = .false.
    return

 10 format('# -----------------------------------------------------------------------------')
 20 format('# PMFLib dynamic binding')
 30 format('# > Unloading driver ...')
 40 format('#   Unable to unload driver (',A,')!')
 50 format('#   Driver unloaded successfully!')

end subroutine pmf_sander_release_driver

!===============================================================================
! subroutine pmf_sander_get_dlerror
!===============================================================================

subroutine pmf_sander_get_dlerror

    implicit none
    type(c_ptr)         :: errstr
    character, pointer  :: strptr(:)
    integer             :: i
    ! --------------------------------------------------------------------------

    pmf_sander_driver_error = ''
    errstr = dlerror()
    call c_f_pointer(errstr,strptr,shape=[strlen(errstr)])

    do i=1,strlen(errstr)
        pmf_sander_driver_error(i:i) = strptr(i)
    end do

end subroutine pmf_sander_get_dlerror

!===============================================================================

#endif

end module PMFLibSANDER_d01

