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

module pmf_sander

use pmf_sizes
use pmf_sander_dat_d01
use iso_c_binding

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_check_interface(rnum,inum,ekin_len,setup_len,str1,str1_len,str2,str2_len) &
           bind(c,name='int_pmf_sander_check_interface')

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

end subroutine pmf_sander_check_interface

!===============================================================================
! subroutine pmf_sander_init_preinit
!===============================================================================

subroutine pmf_sander_init_preinit(mdin,mdin_len,ischeme,anatom,anres, &
                            antb,antc,ansteps,astepsize,atemp0,apress0, &
                            box_a,box_b,box_c,box_alpha,box_beta,box_gamma) &
                            bind(c,name='int_pmf_sander_init_preinit')

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_sander_dat_d01
    use pmf_sander_control_d01
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers

    implicit none
    character(CPMFCHAR) :: mdin(*)
    integer(CPMFINT)    :: mdin_len
    integer(CPMFINT)    :: ischeme
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
    ! -----------------------------------------------
    integer             :: has_box, i
    real(PMFDP)         :: cbox(3)
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_init_preinit'
    end if

    ! setup conversion factors
    MassConv         = 1.0d0        ! g/mol -> g/mol
    LengthConv       = 1.0d0        ! A -> A
    AngleConv        = PMF_D2R      ! deg -> rad
    TimeConv         = 1000.0d0     ! ps -> fs
    VelocityConv     = 1.0d0        ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0        ! kcal/mol -> kcal/mol
    ForceConv        = 1.0d0        ! kcal/mol/A -> kcal/mol/A
    TemperatureConv  = 1.0d0        ! K
    PressureConv     = 1.0d5        ! Bar -> Pa

    ControlFileName = ''
    do i=1,min(mdin_len,len(ControlFileName))
        ControlFileName(i:i)  = mdin(i)
    end do

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    if( ischeme .ne. 0 ) then
        call pmf_init_variables(IA_LF_MIDDLE,anatom,antb,ansteps,astepsize,0.0d0,atemp0,apress0)
    else
        call pmf_init_variables(IA_LEAP_FROG,anatom,antb,ansteps,astepsize,0.0d0,atemp0,apress0)
    end if

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

subroutine pmf_sander_set_residue(idx,name,name_len,first_atom) bind(c,name='int_pmf_sander_set_residue')

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

end subroutine pmf_sander_set_residue

!===============================================================================
! subroutine pmf_sander_set_atom
!===============================================================================

subroutine pmf_sander_set_atom(idx,name,name_len,atype,atype_len) bind(c,name='int_pmf_sander_set_atom')

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

end subroutine pmf_sander_set_atom

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_finalize_preinit(anatom,amass,ax) bind(c,name='int_pmf_sander_finalize_preinit')

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_mask
    use pmf_core
    use pmf_init
    use pmf_sander_control_d01
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: anatom       ! number of atoms
    real(CPMFDP)        :: amass(anatom)
    real(CPMFDP)        :: ax(3,anatom)
    ! -----------------------------------------------
    integer             :: i, alloc_failed
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_finalize_preinit'
    end if

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
    call pmf_init_title('SANDER-DYNAMIC')

    ! load method setups and CVs definition
    call pmf_sander_process_control

    return

end subroutine pmf_sander_finalize_preinit

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_init(anatom,amass,ax) bind(c,name='int_pmf_sander_init')

    use pmf_constants
    use pmf_sizes
    use pmf_init
    use pmf_dat
    use pmf_utils
    use pmf_sander_dat_d01

    implicit none
    integer(CPMFINT)    :: anatom       ! number of atoms
    real(CPMFDP)        :: amass(anatom)
    real(CPMFDP)        :: ax(3,anatom)
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_init'
    end if

    ! init all methods
    call pmf_init_all(amass,ax)

    ! init methods
    call pmf_init_pmf_methods

    return

end subroutine pmf_sander_init

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sender_get_setup(setup,setup_len) bind(c,name='int_pmf_sander_get_setup')

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: setup(:)
    integer(CPMFINT)    :: setup_len
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sender_get_setup'
    end if

    if( setup_len .ne. PMFLIB_SETUP_SIZE ) then
        call pmf_utils_exit(PMF_OUT,1,'Incompatible PMFLIB_SETUP_SIZE - driver interface compromised in pmf_sender_get_setup!')
    end if

    ! currently nothing here
    ignored_arg__ = size(setup) .ne. 0

end subroutine pmf_sender_get_setup

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_sander_shouldexit(exitcode) bind(c,name='int_pmf_sander_shouldexit')

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: exitcode       ! MD loop exit code
    ! --------------------------------------------------------------------------

    ! we cannot bcast the value here as pmf_sander_shouldexit is called on master
    ! protected sections several times per MD cycle
    ! thus fexit_mdloop is distributed in pmf_sander_force_mpi

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_shouldexit'
    end if

    exitcode = fexit_mdloop

    return

end subroutine pmf_sander_shouldexit

!===============================================================================
! subroutine pmf_sander_update_box
!===============================================================================

subroutine pmf_sander_update_box(a,b,c,alpha,beta,gamma) bind(c,name='int_pmf_sander_update_box')

    use pmf_dat
    use pmf_pbc
    use pmf_timers
    use pmf_utils

    real(CPMFDP)    :: a,b,c
    real(CPMFDP)    :: alpha,beta,gamma
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return
    if( fsystype .eq. SYS_NT ) return  ! no box -> return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_update_box'
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_pbc_set_box(a,b,c,alpha,beta,gamma)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_sander_update_box

!===============================================================================
! subroutine pmf_sander_force
!===============================================================================

subroutine pmf_sander_force(anatom,x,v,f,epot,epmf) bind(c,name='int_pmf_sander_force')

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers
    use pmf_dat
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: anatom       ! number of atoms
    real(CPMFDP)        :: x(3,anatom)  ! in
    real(CPMFDP)        :: v(3,anatom)  ! in
    real(CPMFDP)        :: f(3,anatom)  ! inout
    real(CPMFDP)        :: epot         ! in
    real(CPMFDP)        :: epmf         ! out
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_force'
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)

    call pmf_core_lf_update_step
    call pmf_core_lf_force(x,v,f,epot,epmf)

    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_sander_force

!===============================================================================
! subroutine pmf_sander_register_ekin
!===============================================================================

subroutine pmf_sander_register_ekin(ekin) bind(c,name='int_pmf_sander_register_ekin')

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers
    use pmf_dat
    use pmf_utils

    implicit none
    real(CPMFDP)            :: ekin(:)      ! in
    ! --------------------------------------------
    type(PMFKineticEnergy)  :: sekin
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_register_ekin'
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)

    sekin%KinEneVV = ekin(PMFLIB_EKIN_VV)
    sekin%KinEneLF = ekin(PMFLIB_EKIN_LF)
    sekin%KinEneHA = ekin(PMFLIB_EKIN_HA)
    sekin%Valid = .true.

    call pmf_core_lf_register_ekin(sekin)

    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_sander_register_ekin

!===============================================================================
! subroutine pmf_sander_rstforce
! this is used in geometry optimization
!===============================================================================

subroutine pmf_sander_rstforce(anatom,x,f,epot,epmf) bind(c,name='int_pmf_sander_rstforce')

    use pmf_sizes
    use pmf_core_lf
    use pmf_timers
    use pmf_dat
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: anatom        ! in - number of atoms
    real(CPMFDP)        :: x(3,anatom)   ! in
    real(CPMFDP)        :: f(3,anatom)   ! inout
    real(CPMFDP)        :: epot          ! in
    real(CPMFDP)        :: epmf          ! out
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_rstforce'
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)
    call pmf_core_lf_rstforce(x,f,epot,epmf)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

    return

end subroutine pmf_sander_rstforce

!===============================================================================
! Subroutine: pmf_sander_cst_init_collisions
!===============================================================================

subroutine pmf_sander_cst_init_collisions(ntc,nbt,ifstwt,ib,jb,conp) bind(c,name='int_pmf_sander_cst_init_collisions')

    use pmf_dat
    use cst_mdcon_cvs
    use abf_constraints
    use pmf_sizes
    use pmf_core
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: ntc          !
    integer(CPMFINT)    :: nbt          ! number of X-H bonds
    integer(CPMFINT)    :: ifstwt(*)    ! determine if bond is part of water
    integer(CPMFINT)    :: ib(*)        ! the first atom of bond
    integer(CPMFINT)    :: jb(*)        ! the second atom of bond
    real(CPMFDP)        :: conp(*)
    ! -----------------------------------------------
    integer             :: ll, i, j, num, num2
    ! -----------------------------------------------------------------------------

    if( .not. fmaster ) return ! only master can init shake constraint in collision
    if( ntc .eq. 1 ) return    ! no SHAKE

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_cst_init_collisions'
    end if

    if( cst_enabled .or. abf_enabled ) then
        write(PMF_OUT,10)
        write(PMF_OUT,15)
    end if

    if( cst_enabled .and. abf_enabled ) then
        write(PMF_OUT,20)
    end if

! CST --------------------------------------------
    if( cst_enabled ) then
        if( (ntc .ne. 1) .and. (ntc .ne. 2) ) then
            call pmf_utils_exit(PMF_OUT,1,'ntc has to be either one or two for PMFLib constrained dynamics!')
        end if

        ! determine number of SHAKE constraints in collision
        num = 0
        num2 = 0

        do ll = 1,nbt
            i = ib(ll)/3+1
            j  = jb(ll)/3+1
            if( ifstwt(ll) == 1 ) then
                if( cst_mdcon_cvs_checkatom(i,1) .or. cst_mdcon_cvs_checkatom(j,1) ) then
                    call pmf_utils_exit(PMF_OUT,1,'A fast water atom cannot be a part of CST CV')
                end if
            end if
            if( cst_mdcon_cvs_checkatom(i,1) .or. cst_mdcon_cvs_checkatom(j,1) ) then
                num = num + 1
                write(PMF_OUT,100) i, j
            else if( cst_mdcon_cvs_checkatom(i,0) .or. cst_mdcon_cvs_checkatom(j,0) ) then
                num2 = num2 + 1
                write(PMF_OUT,110) i, j
            end if
        end do

        if( num2 .gt. 0 ) then
            write(PMF_OUT,120) num2
            call cst_mdcon_cvs_setnumofexcluded(num2)
        end if

        if( num .gt. 0 ) then
            write(PMF_OUT,130) num

            ! set SHAKE constraints in collisions
            call cst_mdcon_cvs_allocate(num)

            num = 1
            do ll = 1,nbt
                if (ifstwt(ll) == 1) cycle
                i  = ib(ll)/3+1
                j  = jb(ll)/3+1
                if( cst_mdcon_cvs_checkatom(i,1) .or. cst_mdcon_cvs_checkatom(j,1) ) then
                    call cst_mdcon_cvs_set(num,i,j,conp(ll))
                    num = num + 1
                    cycle
                end if
            end do
        else
            write(PMF_OUT,140)
        end if

        write(PMF_OUT,10)

        return
    end if

! ABF --------------------------------------------
    if( abf_enabled ) then
        if( ntc .ne. 2 ) then
            call pmf_utils_exit(PMF_OUT,1,'ntc has to be either one or two for PMFLib ABF dynamics!')
        end if

        ! determine number of SHAKE constraints in collision
        num = 0
        num2 = 0    ! exc but not added to PMFLib

        do ll = 1,nbt
            i = ib(ll)/3+1
            j  = jb(ll)/3+1
            if( ifstwt(ll) == 1 ) then
                if( abf_constraints_checkatom(i,1) .or. abf_constraints_checkatom(j,1) ) then
                    call pmf_utils_exit(PMF_OUT,1,'A fast water atom cannot be a part of ABF CV!')
                end if
            end if
            if( abf_constraints_checkatom(i,1) .or. abf_constraints_checkatom(j,1) ) then
                num = num + 1
                write(PMF_OUT,200) i, j
            else if( abf_constraints_checkatom(i,0) .or. abf_constraints_checkatom(j,0) ) then
                num2 = num2 + 1
                write(PMF_OUT,210) i, j
            end if
        end do

        if( num2 .gt. 0 ) then
            write(PMF_OUT,220) num2
        end if

        if( num .gt. 0 ) then
            write(PMF_OUT,230) num

            ! set SHAKE constraints in collisions
            call abf_constraints_allocate(num)

            num = 1
            do ll = 1,nbt
                if (ifstwt(ll) == 1) cycle
                i  = ib(ll)/3+1
                j  = jb(ll)/3+1
                if( abf_constraints_checkatom(i,1) .or. abf_constraints_checkatom(j,1) ) then
                    call abf_constraints_set(num,i,j,conp(ll))
                    num = num + 1
                    cycle
                end if
            end do
        else
            write(PMF_OUT,240)
        end if

        write(PMF_OUT,10)

        return
    end if

    return

 10 format('# ==============================================================================')
 15 format('# ******************* PMFLib: SHAKE <> CV collisions ***************************')
 20 format('# WARNING: BOTH CST and ABF enabled - CST has priority over ABF!')

100 format("# CST: SHAKEn bond (",I6,"-",I6,") excluded from MD engine and added into PMFLib CST CV list")
110 format("# CST: SHAKEn bond (",I6,"-",I6,") excluded from MD engine")
120 format('# CST: Number of MD constraints excluded: ', I6)
130 format('# CST: Number of MD constraints in collision: ', I6)
140 format('# CST: No MD constraints in collision!')

200 format("# ABF: SHAKEn bond (",I6,"-",I6,") added into PMFLib CV list")
210 format("# ABF: SHAKEn bond (",I6,"-",I6,") excluded from MD engine")
220 format('# ABF: Number of SHAKE constraints removed: ', I6)
230 format('# ABF: Number of SHAKE constraints in collision: ', I6)
240 format('# ABF: No SHAKE in collision added to PMFLib!')

end subroutine pmf_sander_cst_init_collisions

!===============================================================================
! Subroutine: pmf_sander_num_of_pmflib_cst
!===============================================================================

subroutine pmf_sander_num_of_pmflib_cst(numofcst) bind(c,name='int_pmf_sander_num_of_pmflib_cst')

    use cst_dat
    use pmf_utils

    implicit none
    real(CPMFDP)    :: numofcst       ! number of CST constraints
    ! --------------------------------------------------------------------------

    numofcst = 0
    if( .not. cst_enabled ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_num_of_pmflib_cst'
    end if

    numofcst = NumOfCONs - NumOfExcMDCONs ! we must return excluded MD constraints to MD
    return

end subroutine pmf_sander_num_of_pmflib_cst

!===============================================================================
! Subroutine: pmf_sander_cst_checkatom
! this instruct MD engine to remove SHAKE constraint from its internal lists
!===============================================================================

function pmf_sander_cst_checkatom(atomid) bind(c,name='int_pmf_sander_cst_checkatom')

    use pmf_dat
    use cst_mdcon_cvs
    use abf_constraints

    implicit none
    integer(CPMFINT)    :: atomid
    integer(CPMFINT)    :: pmf_sander_cst_checkatom
    ! --------------------------------------------------------------------------

    pmf_sander_cst_checkatom = 0

    if( cst_enabled ) then
        if( cst_mdcon_cvs_checkatom(atomid,0) ) pmf_sander_cst_checkatom = 1
    end if

    if( abf_enabled ) then
        if( abf_constraints_checkatom(atomid,0) ) pmf_sander_cst_checkatom = 1
    end if

    return

end function pmf_sander_cst_checkatom

!===============================================================================
! Subroutine: pmf_sander_shake
!===============================================================================

subroutine pmf_sander_shake(anatom,xp,modified) bind(c,name='int_pmf_sander_shake')

    use pmf_sizes
    use pmf_dat
    use pmf_timers
    use pmf_utils
    use pmf_core_lf

    implicit none
    integer(CPMFINT)    :: anatom            ! number of atoms
    real(CPMFDP)        :: xp(3,anatom)      ! positions in t+dt
    integer(CPMFINT)    :: modified          ! was constraint applied?
    ! --------------------------------------------------------------------------

    modified = 0
    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_shake'
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)
        call pmf_core_lf_shake(xp,modified)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_sander_shake

!===============================================================================
! Subroutine: pmf_sander_rattlev
!===============================================================================

subroutine pmf_sander_rattlev(anatom,xp,vp,cid,modified) bind(c,name='int_pmf_sander_rattlev')

    use pmf_sizes
    use pmf_dat
    use pmf_timers
    use pmf_utils
    use pmf_core_lf

    implicit none
    integer(CPMFINT)    :: anatom            ! number of atoms
    real(CPMFDP)        :: xp(3,anatom)      ! positions at FIXME
    real(CPMFDP)        :: vp(3,anatom)      ! velocities at FIXME
    integer(CPMFINT)    :: cid               ! callid
    integer(CPMFINT)    :: modified          ! was constraint applied?
    ! --------------------------------------------------------------------------

    modified = 0
    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_rattlev'
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)
        call pmf_core_lf_rattlev(xp,vp,cid,modified)
    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_sander_rattlev

!===============================================================================
! subroutine pmf_sander_finalize
!===============================================================================

subroutine pmf_sander_finalize() bind(c,name='int_pmf_sander_finalize')

    use pmf_constants
    use pmf_utils
    use pmf_dat
    use pmf_finalize
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fdebug ) then
        write(PMF_DEBUG+fmytaskid,*) '>>TR: pmf_sander_finalize'
    end if

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'PMF Library Finalization', '-')

    call pmf_timers_stop_timer(PMFLIB_TOTAL_TIMER)
    call pmf_finalize_all(.true.)

end subroutine pmf_sander_finalize

!===============================================================================

#ifdef MPI

!===============================================================================
! subroutine pmf_sander_force_mpi
!===============================================================================

subroutine pmf_sander_force_mpi(anatom,x,v,f,epot,epmf) bind(c,name='int_pmf_sander_force_mpi')

    use pmf_sizes
    use pmf_core_lf
    use pmf_sander_dat_d01
    use pmf_dat
    use pmf_timers
    use pmf_utils
    use mpi

    implicit none
    integer(CPMFINT)    :: anatom       ! number of atoms
    real(CPMFDP)        :: x(3,anatom)  ! in
    real(CPMFDP)        :: v(3,anatom)  ! in
    real(CPMFDP)        :: f(3,anatom)  ! inout
    real(CPMFDP)        :: epot         ! in
    real(CPMFDP)        :: epmf         ! out
    ! ------------------------------------------------------
    integer             :: i,ierr
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
        call pmf_core_lf_force(tmp_a,tmp_b,tmp_c,epot,epmf)
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
! Subroutine: pmf_sander_shake_mpi
!===============================================================================

subroutine pmf_sander_shake_mpi(anatom,x,modified) bind(c,name='int_pmf_sander_shake_mpi')

    use pmf_sizes
    use pmf_dat
    use pmf_core_lf
    use pmf_sander_dat_d01
    use pmf_dat
    use pmf_timers
    use pmf_utils

    implicit none
    integer(CPMFINT)    :: anatom
    real(CPMFDP)        :: x(3,anatom)          ! positions in t+dt
    integer(CPMFINT)    :: modified             ! was constraint applied?
    ! ------------------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    modified = 0
    if( .not. cst_enabled ) return

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

    if(fmaster) then
        call pmf_timers_stop_timer(PMFLIB_TIMER)
    end if

    modified = 1

end subroutine pmf_sander_shake_mpi

!===============================================================================
! subroutine pmf_sander_bcast_dat_mpi
!===============================================================================

subroutine pmf_sander_bcast_dat_mpi(anatom,anumtasks,aiparpt) bind(c,name='int_pmf_sander_bcast_dat_mpi')

    use pmf_init
    use pmf_sander_dat_d01
    use pmf_dat
    use pmf_utils

    integer(CPMFINT)    :: anatom
    integer(CPMFINT)    :: anumtasks
    integer(CPMFINT)    :: aiparpt(:)
    ! -----------------------------------------------
    integer             :: alloc_failed
    integer             :: i,j
    integer             :: istart,iend
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
! subroutine pmf_sander_init_taskid_mpi
!===============================================================================

subroutine pmf_sander_init_taskid_mpi(mytaskid,numoftasks) bind(c,name='int_pmf_sander_init_taskid_mpi')

    use pmf_init

    implicit none
    integer(CPMFINT)    :: mytaskid
    integer(CPMFINT)    :: numoftasks
    ! --------------------------------------------------------------------------

    call pmf_init_taskid_mpi(mytaskid,numoftasks)

end subroutine pmf_sander_init_taskid_mpi

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




