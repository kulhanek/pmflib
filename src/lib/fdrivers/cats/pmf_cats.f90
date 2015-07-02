!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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
! subroutine pmf_cats_begin_init
!===============================================================================

subroutine pmf_cats_begin_init(mdin,anatom,anres, &
                            antb,box_a,box_b,box_c,box_alpha,box_beta,box_gamma)

    use pmf_utils
    use pmf_dat
    use pmf_init
    use pmf_utils
    use pmf_pbc
    use pmf_core
    use pmf_mask
    use pmf_timers
    use pmf_cats_dat
    use pmf_cats_control

    implicit none
    character(*)   :: mdin                         ! control file name
    integer        :: anatom                       ! number of atoms in AMBER topology
    integer        :: anres                        ! number of residues in AMBER topology
    integer        :: antb                         ! BOX type
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
    TimeConv         = 1.0d0               ! fs -> fs
    VelocityConv     = 1.0d0               ! pmflib velocity -> pmflib velocity
    EnergyConv       = 1.0d0               ! kcal/mol -> kcal/mol
    ForceConv        = 1.0d0               ! kcal/mol/A -> kcal/mol/A

    ControlFileName  = mdin

    ! init timers
    call pmf_timers_init_top()
    call pmf_timers_start_timer(PMFLIB_TOTAL_TIMER)

    ! init basic PMF setup
    call pmf_init_dat()
    call pmf_init_variables(IA_LEAP_FROG,anatom,antb,1,1.0d0,0.0d0,300.0d0)
    call pmf_pbc_set_box(box_a,box_b,box_c,box_alpha,box_beta,box_gamma)

    ! init mask subsystem
    call pmf_pbc_get_cbox(has_box,cbox)
    call pmf_mask_topo_init(anatom,anres,has_box,cbox(1),cbox(2),cbox(3))

    ! rewrite the setup
    fcanexmdloop     = .true.       ! the client is able to terminate md loop

end subroutine pmf_cats_begin_init

!===============================================================================
! subroutine pmf_cats_set_residue
!===============================================================================

subroutine pmf_cats_set_residue(idx,name,first_atom)

    use pmf_mask
    use pmf_dat

    implicit none
    integer            :: idx
    character(*)       :: name
    integer            :: first_atom
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    call pmf_mask_set_topo_residue(idx,name,first_atom)

end subroutine pmf_cats_set_residue

!===============================================================================
! subroutine pmf_cats_set_atom
!===============================================================================

subroutine pmf_cats_set_atom(idx,name,atype)

    use pmf_mask

    implicit none
    integer            :: idx
    character(*)       :: name
    character(*)       :: atype
    ! --------------------------------------------------------------------------

    call pmf_mask_set_topo_atom(idx,name,atype)

end subroutine pmf_cats_set_atom

!===============================================================================
! subroutine pmf_cats_end_init
!===============================================================================

subroutine pmf_cats_end_init(anatom,amass,ax)

    use pmf_constants
    use pmf_sizes
    use pmf_dat
    use pmf_mask
    use pmf_core
    use pmf_init
    use pmf_cats_control

    implicit none
    integer        :: anatom                       ! number of atoms in AMBER topology
    real(PMFDP)    :: amass(anatom)
    real(PMFDP)    :: ax(3,anatom)
    ! ------------------------------------------------------
    integer        :: i
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    ! init mask topology atom masses and positions
    do i=1,fnatoms
        call pmf_mask_set_topo_atom_mcrd(i,amass(i),ax(1,i),ax(2,i),ax(3,i))
    end do

    ! finalize mask
    call pmf_mask_topo_finalize()

    ! print PMFLib header with system description
    call pmf_init_title('CATs')

    ! load method setups and CVs definition
    call pmf_cats_process_control

    ! init PMF subsystems
    call pmf_init_all(amass,ax)

    ! flush output streams
    flush(PMF_OUT)

end subroutine pmf_cats_end_init

!===============================================================================
! subroutine pmf_cats_update_box
!===============================================================================

subroutine pmf_cats_update_box(a,b,c,alpha,beta,gamma)

    use pmf_dat
    use pmf_pbc

    real(PMFDP)    :: a,b,c
    real(PMFDP)    :: alpha,beta,gamma
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( fsystype .eq. SYS_NT ) return  ! no box -> return

    call pmf_pbc_set_box(a,b,c,alpha,beta,gamma)

end subroutine pmf_cats_update_box

!===============================================================================
! subroutine pmf_cats_update_x
!===============================================================================

subroutine pmf_cats_update_x(natoms,x)

    use pmf_sizes
    use pmf_core_lf
    use pmf_dat
    use pmf_timers
    use pmf_utils

    implicit none
    integer         :: natoms
    real(PMFDP)     :: x(3,natoms)
    integer         :: i, ridx
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( natoms .ne. fnatoms ) then
        call pmf_utils_exit(PMF_OUT,1,'Inconsistent number of atoms in pmf_cats_update_x!')
    end if

    call pmf_timers_start_timer(PMFLIB_TIMER)
    do i=1,NumOfLAtoms
        ridx = RIndexes(i)
        Crd(:,i) = x(:,ridx)*LengthConv
        Vel(:,i) = 0.0d0
        Frc(:,i) = 0.0d0
    end do

    call pmf_timers_start_timer(PMFLIB_CVS_TIMER)
        ! clear previous data
        CVContext%CVsValues(:) = 0.0d0
        CVContext%CVsDrvs(:,:,:) = 0.0d0

        ! update CVs
        do i=1,NumOfCVs
            call CVList(i)%cv%calculate_cv(Crd,CVContext)
        end do
    call pmf_timers_stop_timer(PMFLIB_CVS_TIMER)

    call pmf_timers_stop_timer(PMFLIB_TIMER)

end subroutine pmf_cats_update_x

!===============================================================================
! subroutine pmf_cats_get_num_of_cvs
!===============================================================================

subroutine pmf_cats_get_num_of_cvs(numofcvsout)

    use pmf_dat
    use pmf_pbc

    integer    :: numofcvsout
    ! --------------------------------------------------------------------------

    numofcvsout = NumOfCVs

end subroutine pmf_cats_get_num_of_cvs

!===============================================================================
! subroutine pmf_cats_get_value
!===============================================================================

subroutine pmf_cats_get_value(value,name)

    use pmf_dat
    use cv_common

    real(PMFDP)     :: value
    character(*)    :: name
    integer         :: indx
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    indx = cv_common_find_cv(name)
    value = CVContext%CVsValues(indx)
    value = CVList(indx)%cv%get_rvalue(value)

end subroutine pmf_cats_get_value

!===============================================================================
! subroutine pmf_cats_get_value_by_indx
!===============================================================================

subroutine pmf_cats_get_value_by_indx(value,indx)

    use pmf_dat
    use cv_common
    use pmf_utils

    real(PMFDP)     :: value
    integer         :: indx
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( (indx .le. 0) .or. (indx .gt. NumOfCVs) ) then
        call pmf_utils_exit(PMF_OUT,1,'CV index is out of legal range!')
    end if

    value = CVContext%CVsValues(indx)
    value = CVList(indx)%cv%get_rvalue(value)

end subroutine pmf_cats_get_value_by_indx

!===============================================================================
! subroutine pmf_cats_get_name
!===============================================================================

subroutine pmf_cats_get_name(name,indx)

    use pmf_dat
    use pmf_utils

    character(*)    :: name
    integer         :: indx
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( (indx .le. 0) .or. (indx .gt. NumOfCVs) ) then
        call pmf_utils_exit(PMF_OUT,1,'CV index is out of legal range!')
    end if

    name = CVList(indx)%cv%name

end subroutine pmf_cats_get_name

!===============================================================================
! subroutine pmf_cats_get_type
!===============================================================================

subroutine pmf_cats_get_type(ctype,name)

    use pmf_dat
    use cv_common
    use pmf_utils

    character(*)    :: ctype
    character(*)    :: name
    integer         :: indx
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    indx = cv_common_find_cv(name)
    ctype = CVList(indx)%cv%ctype

end subroutine pmf_cats_get_type

!===============================================================================
! subroutine pmf_cats_get_type_by_indx
!===============================================================================

subroutine pmf_cats_get_type_by_indx(ctype,indx)

    use pmf_dat
    use pmf_utils

    character(*)    :: ctype
    integer         :: indx
    ! --------------------------------------------------------------------------

    if( .not. fmaster ) return

    if( (indx .le. 0) .or. (indx .gt. NumOfCVs) ) then
        call pmf_utils_exit(PMF_OUT,1,'CV index is out of legal range!')
    end if

    ctype = CVList(indx)%cv%ctype

end subroutine pmf_cats_get_type_by_indx

!===============================================================================
! subroutine pmf_cats_finalize
!===============================================================================

subroutine pmf_cats_finalize

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

end subroutine pmf_cats_finalize

!===============================================================================




