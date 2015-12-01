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
! Subroutine: initialization
!===============================================================================

subroutine initialization

 use pmf_utils
 use pmf_constants
 use pmf_init
 use pmf_pbc
 use pmf_mask
 use pmf_core
 use test_coords_dat

 implicit none
 integer        :: alloc_failed
 integer        :: has_box
 real(PMFDP)    :: cbox(3)
 ! -----------------------------------------------------------------------------

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
 call pmf_init_variables(IA_LEAP_FROG,max_atoms,0,1,0.001d0,0.0d0,300.0d0)

 ! init mask subsystem
 call pmf_pbc_get_cbox(has_box,cbox)
 call pmf_mask_topo_init(max_atoms,0,has_box,cbox(1),cbox(2),cbox(3))

 ! finalize mask
 call pmf_mask_topo_finalize()

 ! print PMFLib header with system description
 call pmf_init_title('PMEMD')

 ! general arrays --------------------------------
 allocate(lx(3,fnatoms),                     &
          loc_x(3,fnatoms),                  &
          mass(fnatoms),                     &
          fdx_ana(3,fnatoms),                &
          fdx_num(3,fnatoms),                &
          tmp_context%CVsValues(1),         &
          tmp_context%CVsDrvs(3,fnatoms,1),  &
          symbols(fnatoms),                  &
          stat= alloc_failed )

 if( alloc_failed .ne. 0 ) then
    call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for arrays!')
 end if

 if( UseExternalSnapshosts ) then
    call open_xyz_stream
 end if

 return

end subroutine initialization

!===============================================================================
! Subroutine: test_coord_deriv
!===============================================================================

subroutine test_coord_deriv

 use test_coords_dat

 implicit none
 integer            :: i
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,*)

 ! init cv indexes ----------------------------
 do i=1,rc1%natoms
    rc1%lindexes(i) = rc1%rindexes(i)
 end do

 ! run tests -------------------------------------
 do i=1,num_of_tests
    call run_coord_test_deriv(i)
 end do

end subroutine test_coord_deriv

!===============================================================================
! Subroutine: test_coord_identity
!===============================================================================

subroutine test_coord_identity

 use test_coords_dat

 implicit none
 integer            :: i
 ! -----------------------------------------------------------------------------

 write(PMF_OUT,*)

 ! init cv indexes ----------------------------
 do i=1,rc1%natoms
    rc1%lindexes(i) = rc1%rindexes(i)
 end do

 do i=1,rc2%natoms
    rc2%lindexes(i) = rc2%rindexes(i)
 end do

 ! run tests -------------------------------------
 do i=1,num_of_tests
    call run_coord_test_identity(i)
 end do

end subroutine test_coord_identity

!===============================================================================
! Subroutine: run_coord_test_deriv(test_id)
!===============================================================================

subroutine run_coord_test_deriv(test_id)

 use test_coords_dat
 use pmf_cvs
 use pmf_utils

 implicit none
 integer            :: test_id
 ! -----------------------------------------------
 real(PMFDP)        :: value_a,value_n
 logical            :: passed
 logical            :: no_more_snapshots
 ! -----------------------------------------------------------------------------

 total_num_of_tests = total_num_of_tests + 1

 ! reset gradients and cvs --------------------
 lx(:,:) = 0.0d0
 loc_x(:,:) = 0.0d0
 mass(:) = 0.0d0
 fdx_ana(:,:) = 0.0d0
 fdx_num(:,:) = 0.0d0
 tmp_context%CVsValues(:) = 0.0d0
 tmp_context%CVsDrvs(:,:,:) = 0.0d0

 ! generate new cvs ---------------------------
 call new_xyz(no_more_snapshots)
 if( no_more_snapshots .eqv. .true. ) return

 ! analytical gradient ---------------------------
 call rc1%calculate_cv(lx,tmp_context)
 value_a = tmp_context%CVsValues(1)
 fdx_ana(:,:) = tmp_context%CVsDrvs(:,:,1)

 ! numerical gradient ----------------------------
 call numerical_derivatives(rc1,lx,value_n,fdx_num)

 passed = compare_gradients(fdx_ana,fdx_num)

 if( passed ) then
    write(PMF_OUT,10) test_id,value_a
 else
    write(PMF_OUT,20) test_id,value_a
    num_of_failed_tests = num_of_failed_tests + 1
    if( verbose_print ) then
        write(PMF_OUT,*)
        write(PMF_OUT,30)
        call print_xyz
        write(PMF_OUT,*)
        write(PMF_OUT,40)
        call print_grad(fdx_ana)
        write(PMF_OUT,*)
        write(PMF_OUT,50)
        call print_grad(fdx_num)
        write(PMF_OUT,*)
    end if
    if( stop_when_error ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> INFO: Requested stop in the case of error.')
    end if
 end if

 return

10 format('Test ID = ',I10,1X,'Value = ',E20.12, ' PASSED')
20 format('Test ID = ',I10,1X,'Value = ',E20.12, ' FAILED')

30 format('Coordinates:')
40 format('Analytical gradients:')
50 format('Numerical gradients:')

end subroutine run_coord_test_deriv

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
 ! -----------------------------------------------------------------------------

 total_num_of_tests = total_num_of_tests + 1

 ! reset gradients and cvs --------------------
 lx(:,:) = 0.0d0
 loc_x(:,:) = 0.0d0
 mass(:) = 0.0d0
 fdx_ana(:,:) = 0.0d0
 fdx_num(:,:) = 0.0d0
 tmp_context%CVsValues(:) = 0.0d0
 tmp_context%CVsDrvs(:,:,:) = 0.0d0

 ! generate new cvs ---------------------------
 call new_xyz(no_more_snapshots)
 if( no_more_snapshots .eqv. .true. ) return

 ! analytical gradient ---------------------------
 call rc1%calculate_cv(lx,tmp_context)
 value_1 = tmp_context%CVsValues(1)
 call rc2%calculate_cv(lx,tmp_context)
 value_2 = tmp_context%CVsValues(1)

 if( abs(value_1-value_2) .lt. alarm_treshold ) then
    write(PMF_OUT,10) test_id,value_1,value_2
 else
    write(PMF_OUT,20) test_id,value_1,value_2
    num_of_failed_tests = num_of_failed_tests + 1
    if( verbose_print ) then
        write(PMF_OUT,*)
        write(PMF_OUT,30)
        call print_xyz
        write(PMF_OUT,*)
    end if
    if( stop_when_error ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> INFO: Requested stop in the case of error.')
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
 ! -----------------------------------------------------------------------------
 integer            :: j,k
 real(PMFDP)        :: mv1,mv2
 logical            :: processed
 ! -----------------------------------------------------------------------------

 ! calculate value
 call cst_item%calculate_cv(x,tmp_context)
 value = tmp_context%CVsValues(1)

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

            call cst_item%calculate_cv(loc_x,tmp_context)
            mv1 = tmp_context%CVsValues(1)

            ! left point ----------------------------------------
            loc_x = lx
            loc_x(k,cst_item%lindexes(j)) = loc_x(k,cst_item%lindexes(j)) - num_diff

            call cst_item%calculate_cv(loc_x,tmp_context)
            mv2 = tmp_context%CVsValues(1)

            fd(k,cst_item%lindexes(j)) = fd(k,cst_item%lindexes(j)) + (mv1 - mv2)/(2.0*num_diff)

        end do
    end if
 end do

end subroutine numerical_derivatives

!===============================================================================
! Subroutine: compare_gradients
!===============================================================================

logical function compare_gradients(fd1,fd2)

 use test_coords_dat

 implicit none
 real(PMFDP)        :: fd1(:,:)
 real(PMFDP)        :: fd2(:,:)
 ! -----------------------------------------------------------------------------
 integer            :: i,k
 ! -----------------------------------------------------------------------------

 compare_gradients = .true.

 do i=1,fnatoms
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
 ! -----------------------------------------------------------------------------

 ! print cvs ----------------------------------
 write(PMF_OUT,*)
 write(PMF_OUT,10)
 write(PMF_OUT,20)
 do i=1,fnatoms
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

subroutine print_xyz

 use test_coords_dat

 implicit none
 integer            :: i
 ! -----------------------------------------------------------------------------

 ! print cvs ----------------------------------
 write(PMF_OUT,*)
 write(PMF_OUT,10)
 write(PMF_OUT,20)
 do i=1,rc1%natoms
    write(PMF_OUT,30) rc1%lindexes(i),mass(rc1%lindexes(i)), &
                        lx(1,rc1%lindexes(i)),  &
                        lx(2,rc1%lindexes(i)),  &
                        lx(3,rc1%lindexes(i))
 end do

 return

10 format(' Atom    Mass       X         Y          Z     ')
20 format('------ -------- ---------- ---------- ----------')
30 format(I6,1X,F8.3,1X,F10.4,1X,F10.4,1X,F10.4)

end subroutine print_xyz

!===============================================================================
! Subroutine: print_grad
!===============================================================================

subroutine print_grad(fd)

 use test_coords_dat

 implicit none
 real(PMFDP)        :: fd(:,:)
 ! -----------------------------------------------
 integer            :: i
 ! -----------------------------------------------------------------------------

 ! print cvs ----------------------------------
 write(PMF_OUT,*)
 write(PMF_OUT,10)
 write(PMF_OUT,20)

 do i=1,rc1%natoms
    write(PMF_OUT,30) rc1%lindexes(i),mass(rc1%lindexes(i)), &
                        fd(1,rc1%lindexes(i)),  &
                        fd(2,rc1%lindexes(i)),  &
                        fd(3,rc1%lindexes(i))
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
 ! -----------------------------------------------------------------------------

 no_more_snapshots = .false.

 if( UseExternalSnapshosts ) then
    call read_xyz_stream(no_more_snapshots)
 else
    ! generate random coordinates
    do i=1,fnatoms
        call rand_xyz(max_radius,lx(1,i),lx(2,i),lx(3,i))
    end do
 end if

 ! generate random masses
 do i=1,fnatoms
    if( uniform_mass .eq. 0 ) then
        call rand_mass(max_mass,mass(i))
    else
        mass(i) = uniform_mass
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
 ! -----------------------------------------------------------------------------

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
 ! -----------------------------------------------------------------------------

 mass = 1.0d0 + mass_limit * rand_uniform()

end subroutine rand_mass

!===============================================================================
! Function: rand_uniform
!===============================================================================

real(PMFDP) function rand_uniform()

 implicit none
 real(PMFSP)   vector(1)
 ! -----------------------------------------------------------------------------

 call RANLUX(vector,1)
 rand_uniform = vector(1)

end function rand_uniform

!===============================================================================
! Subroutine: print_usage
!===============================================================================

subroutine print_usage

 implicit none
 ! -----------------------------------------------------------------------------

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
 !------------------------------------------------------------------------------

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
 integer                            :: i,lnatom
 character(len=80)                  :: xyz_comment
 !------------------------------------------------------------------------------

 error = .false.

 ! read header
 read(PMF_XYZ_S,fmt=*,err=110,end=110) lnatom
 read(PMF_XYZ_S,fmt=*,err=120,end=110) xyz_comment

 if( lnatom .gt. fnatoms ) then
    error = .true.
    return
 end if

 ! read data
 do i=1,lnatom
    read(PMF_XYZ_S,fmt=*,err=130,end=110) symbols(i), lx(1,i),lx(2,i),lx(3,i)
 end do

 return

! error handlers ---------------------------------
100  error = .false.
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
 !------------------------------------------------------------------------------

 ! close file
 close(PMF_XYZ_S)

 return

end subroutine close_xyz_stream

! ==============================================================================

end module test_coords_utils
