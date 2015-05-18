! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module explore

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine exp_startup

    use pmf_system_control
    use pmf_system
    use pmf_system_dat
    use smf_utils
    use smf_profiling
    use pmf_core
    use pmf_utils
    use pmf_paths

    implicit none
    integer     :: i, alloc_failed
    ! --------------------------------------------------------------------------

    call start_timer(INITIALIZATION_TIMER)

    ! read control data --------
    call process_control

    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Initialising explorer', '-')

    path_item => PathList(1)%path
    nbeads = path_item%nbeads

    write(PMF_OUT,10) NumOfPaths
    write(PMF_OUT,20) ncvs
    write(PMF_OUT,30) nbeads

    if( NumOfPaths .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'Only one PATH can be specified!')
    end if

    if( ncvs .ne. path_item%ncvs ) then
        call pmf_utils_exit(PMF_OUT,1,'Numbers of CVs in the PATH and FES differ!')
    end if

    do i=1,ncvs
        if( trim(cvs(i)%name) .ne. trim(CVList(i)%cv%name) ) then
            call pmf_utils_exit(PMF_OUT,1,'The sequence of CVs in the PATH and FES is not the same!')
        end if
        if( trim(cvs(i)%ctype) .ne. trim(CVList(i)%cv%ctype) ) then
            call pmf_utils_exit(PMF_OUT,1,'The sequence of CVs in the PATH and FES is not the same!')
        end if
    end do

    allocate( values(nbeads), &
              deriv(ncvs,nbeads), &
              uposition(ncvs,nbeads), &
              sposition(ncvs,nbeads), &
              rposition(ncvs,nbeads), &
              stat=alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for STM data!')
    end if

    ! initialize other data ----
    call start_timer(DERIVATIVES_TIMER)
        call calc_fes_deriv
    call stop_timer(DERIVATIVES_TIMER)

    call open_trajectory_file     ! open trajectory file if necessary

    call stop_timer(INITIALIZATION_TIMER)

    return

10 format('Number of paths : ',I7)
20 format('Number of CVs   : ',I7)
30 format('Number of beads : ',I7)

end subroutine exp_startup

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine exp_shutdown

    use pmf_system
    use pmf_system_dat
    use smf_profiling
    use smf_utils

    implicit none
    ! --------------------------------------------------------------------------

    call start_timer(FINALIZATION_TIMER)

    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Finalization', '-')

    ! write final path and results files
    call write_results

    ! close trajectory file
    call close_trajectory_file

    call stop_timer(FINALIZATION_TIMER)

    return

end subroutine exp_shutdown

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine exp_run

    use pmf_system_dat
    use pmf_system
    use smf_utils
    use smf_profiling
    use pmf_paths

    implicit none
    integer     :: b,i
    real(PMFDP) :: new_path_len, mov
    ! --------------------------------------------------------------------------

    call start_timer(CORE_TIMER)

    write(PMF_OUT,*)
    call centered_heading(PMF_OUT,'Path exploration', '-')
    call write_stmout_header

    ! optimize path
    call pmf_paths_optimize_alphas(path_item)
    old_path_len = pmf_paths_get_path_length(path_item)

    ! get values and derivatives
    call start_timer(FORCES_TIMER)
        do b=1,nbeads
            call fes_value_and_deriv(path_item%points(b,:),values(b),deriv(:,b))
        end do
    call stop_timer(FORCES_TIMER)

    call write_trajectory_snap

    do istep = 1, nsteps

        call update_path
        call smooth_path
        call reparametrize_path
        call check_boundaries

        MaxMovement = 0
        AveMovement = 0
        MaxMovementBead = 1
        do b=1,nbeads
            mov = 0
            do i=1,ncvs
                mov = mov + (path_item%points(b,i) - rposition(i,b))**2
            end do
            AveMovement =  AveMovement + mov ! // add mov square
            mov = sqrt(mov)
            if( mov > MaxMovement ) then
                MaxMovement = mov
                MaxMovementBead = b
            end if
        end do
        AveMovement = sqrt(AveMovement / nbeads)

        ! copy new path points to path
        do b=1,nbeads
            path_item%points(b,:) = rposition(:,b)
        end do

        call pmf_paths_optimize_alphas(path_item)
        new_path_len = pmf_paths_get_path_length(path_item)
        path_len_change = new_path_len - old_path_len

        ! get values and derivatives
        call start_timer(FORCES_TIMER)
            do b=1,nbeads
                call fes_value_and_deriv(path_item%points(b,:),values(b),deriv(:,b))
            end do
        call stop_timer(FORCES_TIMER)

        call write_stmout_results
        call write_trajectory_snap

        ! do we print intermediate results
        if( (output_freq .gt. 0) .and. (mod(istep,output_freq) .eq. 0) ) then
            call write_results
        end if

        old_path_len = new_path_len
    end do

    call stop_timer(CORE_TIMER)

    return

end subroutine exp_run

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine update_path

    use pmf_system_dat
    use pmf_core
    use smf_profiling
    use pmf_paths

    implicit none
    integer                    :: b, i
    REAL(PMFDP)                :: ntau,df,mov
    REAL(PMFDP)                :: tau(ncvs)
    REAL(PMFDP)                :: tdir(ncvs)
    ! --------------------------------------------------------------------------

    call start_timer(UPDATE_TIMER)

    ! update boundary points
    if( path_item%fixed(1) .neqv. .true. ) then ! avoid fixed points
        do i=1,ncvs
            mov = deriv(i,1)*stepsize
            if( (path_item%maxmoves(i) .gt. 0) .and. (abs(mov) .gt. path_item%maxmoves(i)) ) then
                mov = sign(path_item%maxmoves(i),mov)
            end if
            uposition(i,1) = path_item%points(1,i) - mov
        end do
    else
        uposition(:,1) = path_item%points(1,:)
    end if

    if( path_item%fixed(nbeads) .neqv. .true. ) then ! avoid fixed points
        do i=1,ncvs
            mov = deriv(i,nbeads)*stepsize
            if( (path_item%maxmoves(i) .gt. 0) .and. (abs(mov) .gt. path_item%maxmoves(i)) ) then
                mov = sign(path_item%maxmoves(i),mov)
            end if
            uposition(i,nbeads) = path_item%points(nbeads,i) - mov
        end do
    else
        uposition(:,nbeads) = path_item%points(nbeads,:)
    end if

    ! update middle points
    do b=2,nbeads-1
        call pmf_paths_get_intpoint_der(path_item,path_item%alphas(b),tau)
        df = dot_product(tau,deriv(:,b))
        ntau = sqrt(sum(tau(:)**2))
        if( ntau .gt. 0 ) then
            tau(:) = tau(:) / ntau
        end if
        tdir(:) = -deriv(:,b) + df*tau(:)
        do i=1,ncvs
            mov = tdir(i)*stepsize
            if( (path_item%maxmoves(i) .gt. 0) .and. (abs(mov) .gt. path_item%maxmoves(i)) ) then
                mov = sign(path_item%maxmoves(i),mov)
            end if
            uposition(i,b) = path_item%points(b,i) + mov
        end do
    end do

    call stop_timer(UPDATE_TIMER)

end subroutine update_path

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine smooth_path

    use pmf_system_dat
    use pmf_core
    use smf_profiling

    implicit none
    integer     :: b
    ! --------------------------------------------------------------------------

    call start_timer(SMOOTH_TIMER)

    if( .not. ((smooth_freq .gt. 0) .and. (mod(istep,smooth_freq) .eq. 0)) ) then
        sposition(:,:) = uposition(:,:)
        call stop_timer(SMOOTH_TIMER)
        return
    end if

    ! smooth the path
    sposition(:,1) = uposition(:,1)
    sposition(:,nbeads) = uposition(:,nbeads)
    do b=2,nbeads-1
        sposition(:,b) = (1.0d0 - smoothing)*uposition(:,b) + 0.5d0*smoothing*(uposition(:,b-1) + uposition(:,b+1))
    end do

    call stop_timer(SMOOTH_TIMER)

end subroutine smooth_path

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine reparametrize_path

    use pmf_system_dat
    use pmf_core
    use pmf_paths
    use smf_profiling

    implicit none
    integer     :: b,i
    real(PMFDP) :: new_alpha
    ! --------------------------------------------------------------------------

    call start_timer(REPARAM_TIMER)

    if( .not. ((reparam_freq .gt. 0) .and. (mod(istep,reparam_freq) .eq. 0)) ) then
        rposition(:,:) = sposition(:,:)
        call stop_timer(REPARAM_TIMER)
        return
    end if

    ! copy new path points to path
    do b=1,nbeads
        path_item%points(b,:) = sposition(:,b)
    end do

    ! setup new interpolation
    call pmf_paths_optimize_alphas(path_item)

    ! construct final path
    do b=1,path_item%nbeads
        new_alpha = real(b-1) / real(nbeads-1)
        do i=1,path_item%ncvs
            rposition(i,b) = pmf_paths_get_intcv(path_item,i,new_alpha)
        end do
    end do

    call stop_timer(REPARAM_TIMER)

end subroutine reparametrize_path

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine check_boundaries

    use pmf_system_dat
    use pmf_paths

    implicit none
    integer     :: b,i
    ! --------------------------------------------------------------------------

    do b=1,path_item%nbeads
        do i=1,path_item%ncvs
            if( rposition(i,b) .lt. path_item%minvalues(i) ) then
                rposition(i,b) = path_item%minvalues(i)
            end if
            if( rposition(i,b) .gt. path_item%maxvalues(i) ) then
                rposition(i,b) = path_item%maxvalues(i)
            end if
        end do
    end do

end subroutine check_boundaries

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module explore




