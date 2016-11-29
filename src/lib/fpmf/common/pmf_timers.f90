! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
! ==============================================================================

module pmf_timers

use pmf_sizes
use pmf_constants

implicit none

integer     :: PMFLIB_TOTAL_TIMER                       = -10
    integer     :: PMFLIB_TIMER                         = -20
        integer     :: PMFLIB_CVS_TIMER                 = -30
        integer     :: PMFLIB_PATH_TIMER                = -35
        integer     :: PMFLIB_METHODS_TIMER             = -40
            integer     :: PMFLIB_ABF_TIMER             = -50
                integer     :: PMFLIB_ABF_MWA_TIMER     = -55
            integer     :: PMFLIB_ABP_TIMER             = -56
                integer     :: PMFLIB_ABP_MWA_TIMER     = -57
            integer     :: PMFLIB_MTD_TIMER             = -60
                integer     :: PMFLIB_MTD_MWA_TIMER     = -65
            integer     :: PMFLIB_CON_TIMER             = -70
            integer     :: PMFLIB_RST_TIMER             = -80
            integer     :: PMFLIB_STM_TIMER             = -82
                integer     :: PMFLIB_STM_NET_TIMER     = -84
        integer     :: PMFLIB_EXTENSIONS_TIMER          = -90
            integer     :: PMFLIB_MON_TIMER             = -120
            integer     :: PMFLIB_PDRV_TIMER            = -130
            integer     :: PMFLIB_GAP_TIMER             = -140

contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_timers_init_top

    use smf_profiling
    use pmf_dat

    implicit none
    ! -------------------------------------------------------------------------

    call init_profiling(PMF_OUT,50)
    PMFLIB_TOTAL_TIMER = TOTAL_TIMER

end subroutine pmf_timers_init_top

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_timers_init

<<<<<<< HEAD
 use smf_profiling_dat
 use smf_profiling

 implicit none
 ! -----------------------------------------------------------------------------

 ! add standard timers --------------------------------
 PMFLIB_TIMER           = add_timer(PMFLIB_TOTAL_TIMER,'PMFLib')
    PMFLIB_CVS_TIMER              = add_timer(PMFLIB_TIMER,'Collective Variables')
    PMFLIB_PATH_TIMER              = add_timer(PMFLIB_TIMER,'Paths')
    PMFLIB_METHODS_TIMER          = add_timer(PMFLIB_TIMER,'Methods')
        PMFLIB_ABF_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Adaptive Biasing Force')
            PMFLIB_ABF_MWA_TIMER        = add_timer(PMFLIB_ABF_TIMER,'Multiple Walkers Approach')
        PMFLIB_ABP_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Adaptive Biasing Potential')
            PMFLIB_ABP_MWA_TIMER        = add_timer(PMFLIB_ABP_TIMER,'Multiple Walkers Approach')
        PMFLIB_MTD_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Metadynamics')
            PMFLIB_MTD_MWA_TIMER        = add_timer(PMFLIB_MTD_TIMER,'Multiple Walkers Approach')
        PMFLIB_STM_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'String Method')
            PMFLIB_STM_NET_TIMER        = add_timer(PMFLIB_STM_TIMER,'Network Communication')
        PMFLIB_CON_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Constrained Dynamics')
        PMFLIB_RST_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Restrained Dynamics')
    PMFLIB_EXTENSIONS_TIMER          = add_timer(PMFLIB_TIMER,'Extensions')
        PMFLIB_MON_TIMER        = add_timer(PMFLIB_EXTENSIONS_TIMER,'Monitoring')
        PMFLIB_PDRV_TIMER       = add_timer(PMFLIB_EXTENSIONS_TIMER,'Path Driving')
        PMFLIB_GAP_TIMER        = add_timer(PMFLIB_EXTENSIONS_TIMER,'Gaussian Approximation Potential')
=======
    use smf_profiling_dat
    use smf_profiling

    implicit none
    ! -------------------------------------------------------------------------

    ! add standard timers --------------------------------
    PMFLIB_TIMER           = add_timer(PMFLIB_TOTAL_TIMER,'PMFLib')
        PMFLIB_CVS_TIMER              = add_timer(PMFLIB_TIMER,'Collective Variables')
        PMFLIB_PATH_TIMER              = add_timer(PMFLIB_TIMER,'Paths')
        PMFLIB_METHODS_TIMER          = add_timer(PMFLIB_TIMER,'Methods')
            PMFLIB_ABF_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Adaptive Biasing Force')
                PMFLIB_ABF_MWA_TIMER        = add_timer(PMFLIB_ABF_TIMER,'Multiple Walkers Approach')
            PMFLIB_ABP_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Adaptive Biasing Potential')
                PMFLIB_ABP_MWA_TIMER        = add_timer(PMFLIB_ABP_TIMER,'Multiple Walkers Approach')
            PMFLIB_MTD_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Metadynamics')
                PMFLIB_MTD_MWA_TIMER        = add_timer(PMFLIB_MTD_TIMER,'Multiple Walkers Approach')
            PMFLIB_STM_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'String Method')
                PMFLIB_STM_NET_TIMER        = add_timer(PMFLIB_STM_TIMER,'Network Communication')
            PMFLIB_CON_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Constrained Dynamics')
            PMFLIB_RST_TIMER            = add_timer(PMFLIB_METHODS_TIMER,'Restrained Dynamics')
        PMFLIB_EXTENSIONS_TIMER          = add_timer(PMFLIB_TIMER,'Extensions')
            PMFLIB_MON_TIMER        = add_timer(PMFLIB_EXTENSIONS_TIMER,'Monitoring')
            PMFLIB_PDRV_TIMER       = add_timer(PMFLIB_EXTENSIONS_TIMER,'Path Driving')
            PMFLIB_REMD_TIMER       = add_timer(PMFLIB_EXTENSIONS_TIMER,'Replica Exchange Molecular Dynamics')
            PMFLIB_GAP_TIMER        = add_timer(PMFLIB_EXTENSIONS_TIMER,'Gaussian Approximation Potential')
>>>>>>> 48167413c2c057fa7683287f9ba6c96a630a57a4

end subroutine pmf_timers_init

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_timers_start_timer(id)

    use smf_profiling

    implicit none
    integer        :: id
    ! -------------------------------------------------------------------------

    call start_timer(id)

end subroutine pmf_timers_start_timer

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_timers_stop_timer(id)

    use smf_profiling

    implicit none
    integer        :: id
    ! -------------------------------------------------------------------------

    call stop_timer(id)

end subroutine pmf_timers_stop_timer

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine pmf_timers_finalize(do_profiling)

    use smf_profiling

    implicit none
    logical :: do_profiling
    ! -------------------------------------------------------------------------

    if( do_profiling ) then
        call write_timing
    end if
    call finalize_profiling

end subroutine pmf_timers_finalize

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module pmf_timers
