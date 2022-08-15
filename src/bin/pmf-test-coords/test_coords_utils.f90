! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
!
!     This program is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License along
!     with this program; if not, write to the Free Software Foundation, Inc.,
!     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ==============================================================================

module test_coords_utils

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine: initialization_pass_I
!===============================================================================

subroutine initialization_pass_I

    use pmf_utils
    use pmf_constants
    use pmf_init
    use pmf_pbc
    use pmf_mask
    use pmf_core
    use test_coords_dat

    implicit none
    integer        :: has_box
    real(PMFDP)    :: cbox(3)
    integer        :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init subsystems -------------------------------

! setup conversion factors
    MassConv         = 1.0d0               ! g/mol -> g/mol
    LengthConv       = 1.0d0               ! A -> A
    AngleConv        = PMF_D2R             ! deg -> rad
    TimeConv         = 1000.0d0            ! ps -> fs
    VelocityConv     = 1.0d0               ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0               ! kcal/mol -> kcal/mol
    ForceConv        = 1.0d0               ! kcal/mol/A -> kcal/mol/A

! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_LEAP_FROG,max_atoms,0,1,0.001d0,0.0d0,300.0d0,100000.0d0)

! init mask subsystem
    call pmf_pbc_get_cbox(has_box,cbox)
    call pmf_mask_topo_init(max_atoms,0,has_box,cbox(1),cbox(2),cbox(3))

! finalize mask
    call pmf_mask_topo_finalize()

! print PMFLib header with system description
    call pmf_init_title('TEST')

! init fake frmass
    allocate(frmass(fnatoms),   &
          stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for arrays - initialization_pass_I!')
    end if

    ! this is required by tests in cvs, which test number of atoms by cv_get_group_rmass
        frmass(:) = 1.0d0


 return

end subroutine initialization_pass_I

!===============================================================================
! Subroutine: initialization_pass_II
!===============================================================================

subroutine initialization_pass_II

    use pmf_utils
    use pmf_constants
    use pmf_init
    use pmf_pbc
    use pmf_mask
    use pmf_core
    use test_coords_dat

    implicit none
    integer        :: alloc_failed
    ! --------------------------------------------------------------------------

! init internal arrays
    call pmf_init_pmf

! general arrays --------------------------------
    allocate( symbols(fnatoms),                           &
            rx(3,fnatoms),                                &
            lx(3,NumOfLAtoms),                            &
            loc_x(3,NumOfLAtoms),                         &
            fdx_ana(3,NumOfLAtoms),                       &
            fdx_num(3,NumOfLAtoms),                       &
            tmp_context%CVsValues(NumOfCVs),              &
            tmp_context%CVsDrvs(3,NumOfLAtoms,NumOfCVs),  &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for arrays - initialization_pass_II!')
    end if

    if( UseExternalSnapshosts ) then
        call open_xyz_stream
    end if

end subroutine initialization_pass_II

!===============================================================================
! Subroutine:  test_coordinates
!===============================================================================

subroutine test_coordinates

    use pmf_dat
    use pmf_core
    use pmf_cvs
    use pmf_utils
    use pmf_control
    use test_coords_dat
    use pmf_alloc_cv

    implicit none
    integer         :: i,j
    ! --------------------------------------------------------------------------

    write(PMF_OUT,*)
    call pmf_utils_heading(PMF_OUT,'TESTING', ':')

    if( NumOfCVs .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'No coordinate specified!')
    end if

    write(PMF_OUT,110) NumOfCVs

    select case(test_type)
        case(TEST_DERIV)
            do i=1, NumOfCVs
                write(PMF_OUT,*)
                call pmf_utils_heading(PMF_OUT,'Testing', '-')
                write(PMF_OUT,*)
                write(PMF_OUT,130) i,trim(CVList(i)%cv%ctype)
                write(PMF_OUT,*)

                ! run tests -------------------------------------
                do j=1,num_of_tests
                   call run_coord_test_deriv(i,j)
                end do
            end do
        case(TEST_COMBINED)
            ! run tests -------------------------------------
            do j=1,num_of_tests
               call run_coord_test_deriv_combined(j)
            end do
        case(TEST_IDENTITY)
            if( NumOfCVs .ne. 2 ) then
                call pmf_utils_exit(PMF_OUT,1,'Identity test requires only two coordinate definitions!')
            end if
            ! run tests -------------------------------------
            do i=1,num_of_tests
                call run_coord_test_identity(i)
            end do
        case(TEST_VALUE)
        case default
            call pmf_utils_exit(PMF_OUT,1,'Not implemented type of test!')
    end select

    return

110 format('Number of coordinates          : ',I2)
130 format('=== Coordinate #',I2.2,' of type "',A,'"')

end subroutine test_coordinates

!===============================================================================
! Subroutine: run_coord_test_deriv(test_id)
!===============================================================================

subroutine run_coord_test_deriv(cvid,test_id)

    use test_coords_dat
    use pmf_cvs
    use pmf_utils

    implicit none
    integer            :: cvid,test_id
    ! -----------------------------------------------
    real(PMFDP)        :: value_a,value_n
    logical            :: passed
    logical            :: no_more_snapshots
    ! --------------------------------------------------------------------------

    total_num_of_tests = total_num_of_tests + 1

! reset gradients and cvs --------------------
    lx(:,:) = 0.0d0
    loc_x(:,:) = 0.0d0
    mass(:) = 0.0d0
    frmass(:) = 0.0d0
    fdx_ana(:,:) = 0.0d0
    fdx_num(:,:) = 0.0d0
    tmp_context%CVsValues(:) = 0.0d0
    tmp_context%CVsDrvs(:,:,:) = 0.0d0

! generate new cvs ---------------------------
    call new_xyz(no_more_snapshots)
    if( no_more_snapshots .eqv. .true. ) return

! analytical gradient ---------------------------
    call CVList(cvid)%cv%calculate_cv(lx,tmp_context)
    value_a = tmp_context%CVsValues(cvid)
    fdx_ana(:,:) = tmp_context%CVsDrvs(:,:,cvid)

! numerical gradient ----------------------------
    call numerical_derivatives(CVList(cvid)%cv,lx,value_n,fdx_num)

    passed = compare_gradients(fdx_ana,fdx_num)

    if( passed ) then
        write(PMF_OUT,10) test_id,value_a
    else
        write(PMF_OUT,20) test_id,value_a
        num_of_failed_tests = num_of_failed_tests + 1
        if( verbose_print ) then
            write(PMF_OUT,*)
            write(PMF_OUT,30)
            call print_xyz(cvid)
            write(PMF_OUT,*)
            write(PMF_OUT,40)
            call print_grad(cvid,fdx_ana)
            write(PMF_OUT,*)
            write(PMF_OUT,50)
            call print_grad(cvid,fdx_num)
            write(PMF_OUT,*)
            write(PMF_OUT,60)
            call print_grad(cvid,fdx_ana-fdx_num)
            write(PMF_OUT,*)
        end if
        if( stop_when_error ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested stop in the case of error!')
        end if
    end if

    return

10 format('Test ID = ',I10,1X,'Value = ',E20.12, ' PASSED')
20 format('Test ID = ',I10,1X,'Value = ',E20.12, ' FAILED')

30 format('Coordinates:')
40 format('Analytical gradients:')
50 format('Numerical gradients:')
60 format('Difference between gradients (ana-num):')

end subroutine run_coord_test_deriv

!===============================================================================
! Subroutine: run_coord_test_deriv(test_id)
!===============================================================================

subroutine run_coord_test_deriv_combined(test_id)

    use test_coords_dat
    use pmf_cvs
    use pmf_utils

    implicit none
    integer            :: test_id
    ! -----------------------------------------------
    integer            :: i
    real(PMFDP)        :: value_a,value_n
    logical            :: passed
    logical            :: no_more_snapshots
    ! --------------------------------------------------------------------------

    total_num_of_tests = total_num_of_tests + 1

! reset gradients and cvs --------------------
    lx(:,:) = 0.0d0
    loc_x(:,:) = 0.0d0
    mass(:) = 0.0d0
    frmass(:) = 0.0d0
    fdx_ana(:,:) = 0.0d0
    fdx_num(:,:) = 0.0d0
    tmp_context%CVsValues(:) = 0.0d0
    tmp_context%CVsDrvs(:,:,:) = 0.0d0

! generate new cvs ---------------------------
    call new_xyz(no_more_snapshots)
    if( no_more_snapshots .eqv. .true. ) return

! analytical gradient ---------------------------
    do i=1,NumOfCVs
    ! analytical gradient ---------------------------
        call CVList(i)%cv%calculate_cv(lx,tmp_context)
        value_a = tmp_context%CVsValues(i)
        fdx_ana(:,:) = tmp_context%CVsDrvs(:,:,i)

    ! numerical gradient ----------------------------
        call numerical_derivatives_fromall(i,lx,value_n,fdx_num)

        passed = compare_gradients(fdx_ana,fdx_num)

        if( passed ) then
            write(PMF_OUT,10) test_id,value_a,i,trim(CVList(i)%cv%ctype)
        else
            write(PMF_OUT,20) test_id,value_a,i,trim(CVList(i)%cv%ctype)
            num_of_failed_tests = num_of_failed_tests + 1
            if( verbose_print ) then
                write(PMF_OUT,*)
                write(PMF_OUT,30)
                call print_xyz(i)
                write(PMF_OUT,*)
                write(PMF_OUT,40)
                call print_grad(i,fdx_ana)
                write(PMF_OUT,*)
                write(PMF_OUT,50)
                call print_grad(i,fdx_num)
                write(PMF_OUT,*)
            end if
            if( stop_when_error ) then
                call pmf_utils_exit(PMF_OUT,1,'Requested stop in the case of error!')
            end if
        end if
    end do

    return

10 format('Test ID = ',I10,1X,'Value = ',E20.12,1X,'Coordinate #',I2.2,' of type "',A,'"',' PASSED')
20 format('Test ID = ',I10,1X,'Value = ',E20.12,1X,'Coordinate #',I2.2,' of type "',A,'"',' FAILED')

 30 format('Coordinates:')
 40 format('Analytical gradients:')
 50 format('Numerical gradients:')

end subroutine run_coord_test_deriv_combined

!===============================================================================
! Subroutine: run_coord_test_identity(test_id)
!===============================================================================

subroutine run_coord_test_identity(test_id)

    use test_coords_dat
    use pmf_cvs
    use pmf_utils

    implicit none
    integer            :: test_id
    ! -----------------------------------------------
    real(PMFDP)        :: value_1,value_2
    logical            :: no_more_snapshots
    ! --------------------------------------------------------------------------

    total_num_of_tests = total_num_of_tests + 1

! reset gradients and cvs --------------------
    lx(:,:) = 0.0d0
    loc_x(:,:) = 0.0d0
    mass(:) = 0.0d0
    frmass(:) = 0.0d0
    fdx_ana(:,:) = 0.0d0
    fdx_num(:,:) = 0.0d0
    tmp_context%CVsValues(:) = 0.0d0
    tmp_context%CVsDrvs(:,:,:) = 0.0d0

! generate new cvs ---------------------------
    call new_xyz(no_more_snapshots)
    if( no_more_snapshots .eqv. .true. ) return

! analytical gradient ---------------------------
    call CVList(1)%cv%calculate_cv(lx,tmp_context)
    value_1 = tmp_context%CVsValues(1)
    call CVList(2)%cv%calculate_cv(lx,tmp_context)
    value_2 = tmp_context%CVsValues(2)

    if( abs(value_1-value_2) .lt. alarm_treshold ) then
        write(PMF_OUT,10) test_id,value_1,value_2
    else
        write(PMF_OUT,20) test_id,value_1,value_2
        num_of_failed_tests = num_of_failed_tests + 1
        if( verbose_print ) then
            write(PMF_OUT,*)
            write(PMF_OUT,30)
            call print_xyz(1)
            write(PMF_OUT,*)
            write(PMF_OUT,40)
            call print_grad(1,tmp_context%CVsDrvs(:,:,1))
            write(PMF_OUT,*)
            write(PMF_OUT,30)
            call print_xyz(2)
            write(PMF_OUT,*)
            write(PMF_OUT,40)
            call print_grad(1,tmp_context%CVsDrvs(:,:,2))
        end if
        if( stop_when_error ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested stop in the case of error!')
        end if
    end if

    return

10 format('Test ID = ',I10,1X,'Value1 = ',E14.6,1X,'Value2 = ',E14.6, ' PASSED')
20 format('Test ID = ',I10,1X,'Value1 = ',E14.6,1X,'Value2 = ',E14.6, ' FAILED')

30 format('Coordinates:')
40 format('Analytical gradients:')

end subroutine run_coord_test_identity

!===============================================================================
! Subroutine: numerical_derivatives
!===============================================================================

subroutine numerical_derivatives(cst_item,x,value,fd)

    use test_coords_dat
    use pmf_cvs

    implicit none
    class(CVType)      :: cst_item
    real(PMFDP)        :: x(:,:)
    real(PMFDP)        :: value
    real(PMFDP)        :: fd(:,:)
    ! --------------------------------------------------------------------------
    integer            :: j,k
    real(PMFDP)        :: mv1,mv2
    logical            :: processed
    ! --------------------------------------------------------------------------

! calculate value
    tmp_context%CVsDrvs(:,:,:) = 0.0d0
    tmp_context%CVsValues(:) = 0.0d0
    call cst_item%calculate_cv(x,tmp_context)
    value = tmp_context%CVsValues(cst_item%idx)

! calculate derivatives
    fd(:,:) = 0.0d0
    do j=1,cst_item%natoms
        ! was atom already processed? this is necessary because groups can overlap
        processed = .false.
        do k=1,j-1
            if( cst_item%lindexes(j) .eq. cst_item%lindexes(k) ) then
                processed = .true.
                exit
            end if
        end do
        if( .not. processed ) then
            do k=1,3
                ! right point ----------------------------------------
                loc_x = x
                loc_x(k,cst_item%lindexes(j)) = loc_x(k,cst_item%lindexes(j)) + num_diff

                tmp_context%CVsDrvs(:,:,:) = 0.0d0
                tmp_context%CVsValues(:) = 0.0d0
                call cst_item%calculate_cv(loc_x,tmp_context)
                mv1 = tmp_context%CVsValues(cst_item%idx)

                ! left point ----------------------------------------
                loc_x = lx
                loc_x(k,cst_item%lindexes(j)) = loc_x(k,cst_item%lindexes(j)) - num_diff

                tmp_context%CVsDrvs(:,:,:) = 0.0d0
                tmp_context%CVsValues(:) = 0.0d0
                call cst_item%calculate_cv(loc_x,tmp_context)
                mv2 = tmp_context%CVsValues(cst_item%idx)

                fd(k,cst_item%lindexes(j)) = fd(k,cst_item%lindexes(j)) + (mv1 - mv2)/(2.0*num_diff)

            end do
        end if
    end do

end subroutine numerical_derivatives

!===============================================================================
! Subroutine: numerical_derivatives_fromall
!===============================================================================

subroutine numerical_derivatives_fromall(cvidx,x,value,fd)

    use test_coords_dat
    use pmf_cvs

    implicit none
    integer            :: cvidx
    real(PMFDP)        :: x(:,:)
    real(PMFDP)        :: value
    real(PMFDP)        :: fd(:,:)
    ! --------------------------------------------------------------------------
    integer            :: i,j,k
    real(PMFDP)        :: mv1,mv2
    ! --------------------------------------------------------------------------

! calculate value
    tmp_context%CVsDrvs(:,:,:) = 0.0d0
    do i=1,cvidx
        call CVList(i)%cv%calculate_cv(x,tmp_context)
        value = tmp_context%CVsValues(i)
    end do

! calculate derivatives
    fd(:,:) = 0.0d0
    do j=1,NumOfLAtoms
        do k=1,3
            ! right point ----------------------------------------
            loc_x = x
            loc_x(k,j) = loc_x(k,j) + num_diff

            tmp_context%CVsDrvs(:,:,:) = 0.0d0
            do i=1,cvidx
                call CVList(i)%cv%calculate_cv(loc_x,tmp_context)
                mv1 = tmp_context%CVsValues(i)
            end do

            ! left point ----------------------------------------
            loc_x = lx
            loc_x(k,j) = loc_x(k,j) - num_diff

            tmp_context%CVsDrvs(:,:,:) = 0.0d0
            do i=1,cvidx
                call CVList(i)%cv%calculate_cv(loc_x,tmp_context)
                mv2 = tmp_context%CVsValues(i)
            end do

            fd(k,j) = fd(k,j) + (mv1 - mv2)/(2.0*num_diff)

        end do
    end do

end subroutine numerical_derivatives_fromall

!===============================================================================
! Subroutine: compare_gradients
!===============================================================================

logical function compare_gradients(fd1,fd2)

    use test_coords_dat

    implicit none
    real(PMFDP)        :: fd1(:,:)
    real(PMFDP)        :: fd2(:,:)
    ! --------------------------------------------------------------------------
    integer            :: i,k
    ! --------------------------------------------------------------------------

    compare_gradients = .true.

    do i=1,NumOfLAtoms
        do k=1,3
            if( abs(fd1(k,i)-fd2(k,i)) .gt. alarm_treshold ) then
                compare_gradients = .false.
                return
            end if
        end do
    end do

    return

end function compare_gradients

!===============================================================================
! Subroutine: print_xyz
!===============================================================================

subroutine print_all_xyz

    use test_coords_dat

    implicit none
    integer            :: i
    ! --------------------------------------------------------------------------

    ! print cvs ----------------------------------
    write(PMF_OUT,*)
    write(PMF_OUT,10)
    write(PMF_OUT,20)
    do i=1,NumOfLAtoms
        write(PMF_OUT,30) i,mass(i),lx(1,i),lx(2,i),lx(3,i)
    end do

    return

10 format(' Atom    Mass       X         Y          Z     ')
20 format('------ -------- ---------- ---------- ----------')
30 format(I6,1X,F8.3,1X,F10.4,1X,F10.4,1X,F10.4)

end subroutine print_all_xyz

!===============================================================================
! Subroutine: print_xyz
!===============================================================================

subroutine print_xyz(cvid)

    use test_coords_dat

    implicit none
    integer            :: cvid
    integer            :: i, li
    ! --------------------------------------------------------------------------

! print cvs ----------------------------------
    write(PMF_OUT,*)
    write(PMF_OUT,10)
    write(PMF_OUT,20)
    do i=1,CVList(cvid)%cv%natoms
        li = CVList(cvid)%cv%lindexes(i)
        write(PMF_OUT,30) li,mass(li), &
                            lx(1,li),  &
                            lx(2,li),  &
                            lx(3,li)
    end do

     return

10 format(' Atom    Mass       X         Y          Z     ')
20 format('------ -------- ---------- ---------- ----------')
30 format(I6,1X,F8.3,1X,F10.4,1X,F10.4,1X,F10.4)

end subroutine print_xyz

!===============================================================================
! Subroutine: print_grad
!===============================================================================

subroutine print_grad(cvid,fd)

    use test_coords_dat

    implicit none
    integer            :: cvid
    real(PMFDP)        :: fd(:,:)
    ! -----------------------------------------------
    integer            :: i, li
    ! --------------------------------------------------------------------------

! print grad ----------------------------------
    write(PMF_OUT,*)
    write(PMF_OUT,10)
    write(PMF_OUT,20)

    do i=1,CVList(cvid)%cv%natoms
        li = CVList(cvid)%cv%lindexes(i)
        write(PMF_OUT,30) li,mass(li), &
                            fd(1,li),  &
                            fd(2,li),  &
                            fd(3,li)
    end do

    return

10 format(' Atom    Mass      dX        dY         dZ     ')
20 format('------ -------- ---------- ---------- ----------')
30 format(I6,1X,F8.3,1X,F10.4,1X,F10.4,1X,F10.4)

end subroutine print_grad

!===============================================================================
! Subroutine: new_xyz
!===============================================================================

subroutine new_xyz(no_more_snapshots)

    use test_coords_dat

    implicit none
    logical            :: no_more_snapshots
    ! -------------------------
    integer            :: i
    ! --------------------------------------------------------------------------

    no_more_snapshots = .false.

    if( UseExternalSnapshosts ) then
        call read_xyz_stream(no_more_snapshots)
    else
        ! generate random coordinates
        do i=1,NumOfLAtoms
            call rand_xyz(max_radius,lx(1,i),lx(2,i),lx(3,i))
        end do
    end if

! generate random masses
    do i=1,NumOfLAtoms
        if( uniform_mass .eq. 0 ) then
            call rand_mass(max_mass,mass(i))
            frmass(i) = mass(i)
        else
            mass(i) = uniform_mass
            frmass(i) = uniform_mass
        end if
    end do

end subroutine new_xyz

!===============================================================================
! Subroutine: rand_xyz
!===============================================================================

subroutine rand_xyz(radius,x,y,z)

    use test_coords_dat

    implicit none
    real(PMFDP)       :: radius
    real(PMFDP)       :: x,y,z
    ! -----------------------------------------------
    real(PMFDP)       :: r2
    ! --------------------------------------------------------------------------

    r2 = 0.0d0

    do while( (r2 .gt. 1.0d0) .or. (r2 .eq. 0d0) )
        ! choose x,y,z in uniform cube (-1,-1,-1) to (+1,+1,+1)
        x = -1 + 2 * rand_uniform()
        y = -1 + 2 * rand_uniform()
        z = -1 + 2 * rand_uniform()

        ! see if it is in the unit circle
        r2 = x**2 + y**2 + z**2
    end do

    ! scale
    x = x * radius
    y = y * radius
    z = z * radius

end subroutine rand_xyz

!===============================================================================
! Subroutine: rand_mass
!===============================================================================

subroutine rand_mass(mass_limit,mass)

    implicit none
    real(PMFDP)       :: mass_limit
    real(PMFDP)       :: mass
    ! --------------------------------------------------------------------------

    mass = 1.0d0 + mass_limit * rand_uniform()

end subroutine rand_mass

!===============================================================================
! Function: rand_uniform
!===============================================================================

real(PMFDP) function rand_uniform()

    implicit none
    real(PMFSP)   vector(1)
    ! --------------------------------------------------------------------------

    call RANLUX(vector,1)
    rand_uniform = vector(1)

end function rand_uniform

!===============================================================================
! Subroutine: print_usage
!===============================================================================

subroutine print_usage

    implicit none
    ! --------------------------------------------------------------------------

    write(PMF_OUT,'(/,a,/)') '=== [usage] ===================================================================='

    write(PMF_OUT,10)

    return

10 format('    pmf-test-coords <control_file> <cvs_file>')

end subroutine print_usage

!===============================================================================
! Subroutine: open_xyz_stream
!===============================================================================

subroutine open_xyz_stream

    use pmf_constants
    use test_coords_dat
    use pmf_utils

    implicit none
    !---------------------------------------------------------------------------

! try to open XYZ file
    open(unit=PMF_XYZ_S, file=SnapshotFile, status='unknown', form='formatted', err=100)
    return

! error handlers ---------------------------------
    100  call pmf_utils_exit(PMF_OUT,1,&
            'Opening XYZ stream file failed!')

end subroutine open_xyz_stream

!===============================================================================
! Subroutine: read_xyz_stream
!===============================================================================

subroutine read_xyz_stream(error)

    use pmf_constants
    use test_coords_dat

    implicit none
    logical                            :: error
    ! -----------------------------------------------
    integer                            :: i,nats
    character(len=80)                  :: xyz_comment
    !---------------------------------------------------------------------------

    error = .false.

! read header
    read(PMF_XYZ_S,fmt=*,err=110,end=110) nats
    read(PMF_XYZ_S,fmt=*,err=120,end=110) xyz_comment

    if( nats .ne. fnatoms ) then
        error = .true.
        return
    end if

! read data
    do i=1,nats
        read(PMF_XYZ_S,fmt=*,err=130,end=110) symbols(i), rx(1,i),rx(2,i),rx(3,i)
    end do

! map to local atoms
    do i=1,NumOfLAtoms
        lx(:,i) = rx(:,Rindexes(i))
    end do

    return

! error handlers ---------------------------------
110  error = .false.
120  error = .false.
130  error = .false.

    return

end subroutine read_xyz_stream

!===============================================================================
! Subroutine: close_xyz_stream
!===============================================================================

subroutine close_xyz_stream

    use pmf_constants

    implicit none
    !---------------------------------------------------------------------------

    ! close file
    close(PMF_XYZ_S)

    return

end subroutine close_xyz_stream

! ==============================================================================

end module test_coords_utils
