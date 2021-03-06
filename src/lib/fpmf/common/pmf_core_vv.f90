!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_core_vv

! SFR = ! Shake(t) -> Force(t+dt) -> Rattle (t+dt)

implicit none
contains

!===============================================================================
! Subroutine:  pmf_core_vv_shake_SFR
!===============================================================================

subroutine pmf_core_vv_shake_SFR(xp,vp)

    use pmf_dat
    use pmf_sizes
    use pmf_core
    use cst_core
    use pmf_timers

    implicit none
    real(PMFDP),intent(inout)  :: xp(:,:) ! r_u(t+dt)->R(t+dt)
    real(PMFDP),intent(inout)  :: vp(:,:) ! v_u(t')->??
    ! --------------------------------------------------------------------------

    if( .not. cst_enabled ) return

    call pmf_timers_start_timer(PMFLIB_METHODS_TIMER)
        call pmf_timers_start_timer(PMFLIB_CON_TIMER)

        ! update local data
        call pmf_core_in_data_xpvp(xp,vp)

        call cst_core_main_vv_shake

        ! update global data
        call pmf_core_out_data_xpvp(xp,vp)

        call pmf_timers_stop_timer(PMFLIB_CON_TIMER)
    call pmf_timers_stop_timer(PMFLIB_METHODS_TIMER)

    return

end subroutine pmf_core_vv_shake_SFR

!===============================================================================
! Subroutine:  pmf_core_vv_force_SRF
!===============================================================================

subroutine pmf_core_vv_force_SRF(x,f,epot,epmf)

    use pmf_dat
    use pmf_cvs
    use mon_output
    use rst_core
    use rst_dat
    use abf_core
    use abp_core
    use mtd_core
    use stm_core
    use mtd_dat
    use stm_dat
    use pmf_timers
    use pmf_core
    use pmf_paths
    use pdrv_core

    implicit none
    real(PMFDP),intent(in)     :: x(:,:)       ! position in t
    real(PMFDP),intent(inout)  :: f(:,:)       ! forces in t
    real(PMFDP),intent(in)     :: epot         ! potential energy of system in t
    real(PMFDP),intent(out)    :: epmf         ! energy from PMFLib
    ! -----------------------------------------------
    integer        :: i
    ! --------------------------------------------------------------------------

    epmf = 0.0d0

    if( .not. pmf_enabled ) return

    ! update time increments
    ftime = ftime + fdt
    fstep = fstep + 1

    ! convert potential energy
    PotEne = epot *EnergyConv

    ! update local data
    call pmf_core_in_data_xf(x,f)

!-------------------------------------------------
    call pmf_timers_start_timer(PMFLIB_CVS_TIMER)
        ! clear previous data
        CVContext%CVsValues(:) = 0.0d0
        CVContext%CVsDrvs(:,:,:) = 0.0d0

        ! update CVs
        do i=1,NumOfCVs
            call CVList(i)%cv%calculate_cv(Crd,CVContext)
        end do
    call pmf_timers_stop_timer(PMFLIB_CVS_TIMER)

!-------------------------------------------------
! this is extremely slow - it should be used only upon request
    if( fmonitor_paths ) then
        call pmf_timers_start_timer(PMFLIB_PATH_TIMER)
            ! update PATHs
            do i=1,NumOfPATHs
                call pmf_paths_get_path_current_alpha(PathList(i)%path,CVContext)
            end do
        call pmf_timers_stop_timer(PMFLIB_PATH_TIMER)
    end if

!-------------------------------------------------
    call pmf_timers_start_timer(PMFLIB_METHODS_TIMER)

        if( pdrv_enabled ) then
            call pmf_timers_start_timer(PMFLIB_PDRV_TIMER)
                call pdrv_core_main
            call pmf_timers_stop_timer(PMFLIB_PDRV_TIMER)
        end if

        if( mtd_enabled ) then
            call pmf_timers_start_timer(PMFLIB_MTD_TIMER)
                call mtd_core_main
                epmf = epmf + TotalMTDEnergy
            call pmf_timers_stop_timer(PMFLIB_MTD_TIMER)
        end if

        if( rst_enabled ) then
            call pmf_timers_start_timer(PMFLIB_RST_TIMER)
                call rst_core_main
                epmf = epmf + TotalRstEnergy
            call pmf_timers_stop_timer(PMFLIB_RST_TIMER)
        end if

        if( stm_enabled ) then
            call pmf_timers_start_timer(PMFLIB_STM_TIMER)
                call stm_core_main
                epmf = epmf + TotalSTMEnergy
            call pmf_timers_stop_timer(PMFLIB_STM_TIMER)
        end if

        if( mon_enabled ) then
            call pmf_timers_start_timer(PMFLIB_MON_TIMER)
                call mon_output_write_output
            call pmf_timers_stop_timer(PMFLIB_MON_TIMER)
        end if

        ! ABF has to be here because it could be influenced by restraints (wall restraints, etc.)
        if( abf_enabled ) then
            call pmf_timers_start_timer(PMFLIB_ABF_TIMER)
                call abf_core_main
            call pmf_timers_stop_timer(PMFLIB_ABF_TIMER)
        end if

        ! ABP has to be here because it could be influenced by restraints (wall restraints, etc.)
        if( abp_enabled ) then
            call pmf_timers_start_timer(PMFLIB_ABP_TIMER)
                call abp_core_main
            call pmf_timers_stop_timer(PMFLIB_ABP_TIMER)
        end if

    call pmf_timers_stop_timer(PMFLIB_METHODS_TIMER)

    epmf = epmf / EnergyConv

    ! update forces
    call pmf_core_out_data_f(f)

    ! debug
    ! write(RST_OUT,*) f(:,:)

end subroutine pmf_core_vv_force_SRF

!===============================================================================
! Subroutine:  pmf_core_vv_rattle_SFR
!===============================================================================

subroutine pmf_core_vv_rattle_SFR(vp)

    use pmf_sizes
    use pmf_dat
    use pmf_core
    use cst_core
    use pmf_timers

    implicit none
    real(PMFDP),intent(inout)  :: vp(:,:)
    ! --------------------------------------------------------------------------

    if( .not. cst_enabled ) return

    call pmf_timers_start_timer(PMFLIB_METHODS_TIMER)
        call pmf_timers_start_timer(PMFLIB_CON_TIMER)

        ! update local data
        call pmf_core_in_data_vp(vp)

        ! correct velocities
        call cst_core_main_vv_rattle()

        ! update global data
        call pmf_core_out_data_vp(vp)

        call pmf_timers_stop_timer(PMFLIB_CON_TIMER)
    call pmf_timers_stop_timer(PMFLIB_METHODS_TIMER)

    return

end subroutine pmf_core_vv_rattle_SFR

!===============================================================================

end module pmf_core_vv


