!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more detajls.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module abf_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine abf_core_main

    use abf_trajectory
    use abf_restart
    use abf_output
    use abf_client
    use abf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    call abf_core_update_history

    select case(fmode)
        ! standard algorithms
        case(1)
            call abf_core_force_3pV1
        case(2)
            call abf_core_force_3pV2
        ! standard algorithms
        case(3)
            call abf_core_force_5pV1
        case(4)
            call abf_core_force_5pV2
        case(5)
            call abf_core_force_3pV3
!        ! testing algorithms
!        case(10)
!            call abf_core_force_3pF
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    call abf_client_exchange_data(.false.)

end subroutine abf_core_main

!===============================================================================
! Subroutine:  abf_core_shake
! correct for forces from SHAKE
!===============================================================================

subroutine abf_core_shake

    use abf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------
    select case(fmode)
        case(10)
            ! get forces from SHAKE - abf_core_force_3pF
            shist(:,:,hist_len) = SHAKEFrc(:,:)
        case(1,2)
            ! ignored
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_shake!')
    end select

end subroutine abf_core_shake

!===============================================================================
! Subroutine:  abf_core_update_history
! apply ABF force and update history buffers
!===============================================================================

subroutine abf_core_update_history()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m,ci,ki
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        ekinvvhist(i)       = ekinvvhist(i+1)
        ekinlfhist(i)       = ekinlfhist(i+1)
        vhist(:,:,i)        = vhist(:,:,i+1)
        zdhist(:,:,:,i)     = zdhist(:,:,:,i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)

! shift ene
    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    select case(ftds_ekin_src)
    case(1)
        ekinhist(hist_len-1)    = KinEneVV - fekinaverage    ! shifted by -dt
    case(2)
        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
    case(3)
        ekinhist(hist_len-2)    = KinEneV4 - fekinaverage
    case(4)
        ekinhist(hist_len-3)    = KinEneV6 - fekinaverage
    case(5)
        ekinhist(hist_len-1)    = (KinEneVV - fekinaverage)*(300.0d0/267.63d0)    ! shifted by -dt
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_force_3pV!')
    end select

    !write(74895,*) fstep-2, KinEneV4, fstep-1, KinEneVV

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zdhist(m,j,i,hist_len) = v
            end do
        end do
    end do

! apply force filters
    la(:) = 0.0d0
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvhist(:,hist_len),la)
            case(1)
                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
            case(2)
                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
            case(3)
                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode in abf_core_force_3pV!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if
    micfhist(:,hist_len) = la(:)

    if( frecord ) then
        ! record time progress of data
        call abf_accu_add_data_record_lf(cvhist(:,hist_len),micfhist(:,hist_len), &
                                         epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
    end if

    return

end subroutine abf_core_update_history

!===============================================================================
! Subroutine:  abf_core_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV1()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1,v2
    ! --------------------------------------------------------------------------

    if( fstep .le. hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * (vhist(m,j,hist_len-2) - vhist(m,j,hist_len-3))
                ! velocity part
                v1 = v1 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-3)) * vhist(m,j,hist_len-2)
                v2 = v2 + (zdhist(m,j,i,hist_len-3)-zdhist(m,j,i,hist_len-4)) * vhist(m,j,hist_len-3)
            end do
        end do
        pxif(i) = f1*ifdtx
        pxiv(i) = 0.5d0*(v1+v2)*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV1

!===============================================================================
! Subroutine:  abf_core_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV3()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1,v2,v3,v4
    ! --------------------------------------------------------------------------

    if( fstep .le. hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        v4 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * (vhist(m,j,hist_len-2) - vhist(m,j,hist_len-3))
                ! velocity part
                v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                v3 = v3 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-3)) * vhist(m,j,hist_len-2)
                v4 = v4 + (zdhist(m,j,i,hist_len-3)-zdhist(m,j,i,hist_len-4)) * vhist(m,j,hist_len-3)
            end do
        end do
        pxif(i) = f1*ifdtx
        pxiv(i) = (1.0d0/16.0d0)*(-v1 + 9.0d0*v2 + 9.0d0*v3 -v4)*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV3

!===============================================================================
! Subroutine:  abf_core_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV2()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1
    ! --------------------------------------------------------------------------

    if( fstep .le. hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * (vhist(m,j,hist_len-2) - vhist(m,j,hist_len-3))
                ! velocity part
                v1 = v1 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-4)) * &
                          (vhist(m,j,hist_len-2) + vhist(m,j,hist_len-3))
            end do
        end do
        pxif(i) = f1*ifdtx
        pxiv(i) = 0.25d0*v1*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV2

!===============================================================================
! Subroutine:  abf_core_force_5pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_5pV1()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1,v2,v3,v4
    ! --------------------------------------------------------------------------

    if( fstep .le. hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        v4 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-2) * (- vhist(m,j,hist_len-0) + 27.0d0*vhist(m,j,hist_len-1) &
                                                      - 27.0d0*vhist(m,j,hist_len-2) + vhist(m,j,hist_len-3))
                ! velocity part
                v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                v3 = v3 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-3)) * vhist(m,j,hist_len-2)
                v4 = v4 + (zdhist(m,j,i,hist_len-3)-zdhist(m,j,i,hist_len-4)) * vhist(m,j,hist_len-3)
            end do
        end do
        pxif(i) = (1.0d0/24.0d0)*f1*ifdtx
        pxiv(i) = (1.0d0/16.0d0)*(-v1 + 9.0d0*v2 + 9.0d0*v3 -v4)*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_5pV1

!===============================================================================
! Subroutine:  abf_core_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_5pV2()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1
    ! --------------------------------------------------------------------------

    if( fstep .le. hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-2) * (- vhist(m,j,hist_len-0) + 27.0d0*vhist(m,j,hist_len-1) &
                                                      - 27.0d0*vhist(m,j,hist_len-2) + vhist(m,j,hist_len-3))
                ! velocity part
                v1 = v1 + (     - zdhist(m,j,i,hist_len-0) + 8.0d0*zdhist(m,j,i,hist_len-1) &
                           -8.0d0*zdhist(m,j,i,hist_len-3)       + zdhist(m,j,i,hist_len-4) ) * &
                           (       - vhist(m,j,hist_len-0) + 9.0d0*vhist(m,j,hist_len-1) &
                             + 9.0d0*vhist(m,j,hist_len-2)       - vhist(m,j,hist_len-3))
            end do
        end do
        pxif(i) = (1.0d0/24.0d0)*f1*ifdtx
        pxiv(i) = (1.0d0/192.0d0)*v1*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_5pV2

!!===============================================================================
!! Subroutine:  abf_core_force_3pV4
!! this is leap-frog ABF version, simplified algorithm
!! ICF from velocities + decomposition
!!===============================================================================
!
!subroutine abf_core_force_3pV4()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!
!    implicit none
!    integer                :: i,j,k,m
!    integer                :: ci,ki
!    real(PMFDP)            :: v1,s1,l1,f1
!    ! --------------------------------------------------------------------------
!
!! shift accuvalue history
!    do i=1,hist_len-1
!        cvhist(:,i)         = cvhist(:,i+1)
!        epothist(i)         = epothist(i+1)
!        ersthist(i)         = ersthist(i+1)
!        ekinhist(i)         = ekinhist(i+1)
!        ekinvvhist(i)         = ekinvvhist(i+1)
!        ekinlfhist(i)         = ekinlfhist(i+1)
!        vhist(:,:,i)        = vhist(:,:,i+1)
!        zdhist(:,:,:,i)     = zdhist(:,:,:,i+1)
!        micfhist(:,i)       = micfhist(:,i+1)
!    end do
!
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len) = CVContext%CVsValues(ci)
!    end do
!    vhist(:,:,hist_len)    = Vel(:,:)
!
!! shift ene
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!select case(ftds_ekin_src)
!    case(1)
!        ekinhist(hist_len-1)    = KinEneVV - fekinaverage    ! shifted by -dt
!    case(2)
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
!        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
!    case(3)
!        ! kinetic part
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
!        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
!        ! correction
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) &
!                                + 1.0d0/8.0d0*(epothist(hist_len)+ersthist(hist_len) + epothist(hist_len-2)+ersthist(hist_len-2)) &
!                                - 2.0d0/8.0d0*(epothist(hist_len-1)+ersthist(hist_len-1))
!    case(4)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-3/2dt and t-5/2dt
!                                  ! t-dt, t-2dt, t-3dt
!        ekinhist(hist_len-2)    = 0.5d0*(ekinlfhist(hist_len-2)+ekinlfhist(hist_len-3))  &
!                                - 1.0d0/8.0d0*(ekinvvhist(hist_len-1)-2.0d0*ekinvvhist(hist_len-2)+ekinvvhist(hist_len-3))
!    case(5)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-1/2dt t-3/2dt and t-5/2dt t-7/2dt
!        ekinhist(hist_len-2)    = 1.0d0/10.0d0*(-ekinlfhist(hist_len-1)+6.0d0*ekinlfhist(hist_len-2)&
!                                                +6.0d0*ekinlfhist(hist_len-3)-ekinlfhist(hist_len-4))
!    case(6)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-1/2dt t-3/2dt and t-5/2dt t-7/2dt
!        ekinhist(hist_len-2)    = 1.0d0/10.0d0*(-ekinlfhist(hist_len-1)+6.0d0*ekinlfhist(hist_len-2)&
!                                                +6.0d0*ekinlfhist(hist_len-3)-ekinlfhist(hist_len-4))
!        ekinhist(hist_len-2)    = 0.5d0*(ekinhist(hist_len-2) + ekinvvhist(hist_len-2))
!    case(7)
!        ekinhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!    case(8)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-3/2dt and t-5/2dt
!                                  ! t-dt, t-2dt, t-3dt
!        ekinhist(hist_len-2)    = (4.0d0/6.0d0)*(ekinlfhist(hist_len-2)+ekinlfhist(hist_len-3))  &
!                                - (1.0d0/6.0d0)*(ekinvvhist(hist_len-1)+ekinvvhist(hist_len-3))
!    case(9)
!        ekinhist(hist_len-2)    = KinEneV4 - fekinaverage
!    case default
!        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_force_3pV!')
!    end select
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!
!    do i=1,NumOfABFCVs
!        do j=1,NumOfLAtoms
!            do m=1,3
!                v1 = 0.0d0
!                do k=1,NumOfABFCVs
!                    ki = ABFCVList(k)%cvindx
!                    v1 = v1 + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
!                end do
!                zdhist(m,j,i,hist_len) = v1
!            end do
!        end do
!    end do
!
!! apply force filters
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode in abf_core_force_3pV!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len) = la(:)
!
!    if( frecord ) then
!        ! record time progress of data
!        call abf_accu_add_data_record_lf(cvhist(:,hist_len),micfhist(:,hist_len), &
!                                         epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!     end if
!
!! ABF part
!    if( fstep .ge. hist_len ) then
!        do i=1,NumOfABFCVs
!            f1 = 0.0d0
!            s1 = 0.0d0 ! zero - partly in f1
!            l1 = 0.0d0 ! zero
!            v1 = 0.0d0
!            do j=1,NumOfLAtoms
!                do m=1,3
!                    ! force part
!                    f1 = f1 + zdhist(m,j,i,hist_len-2) * (- vhist(m,j,hist_len-0) + 27.0d0*vhist(m,j,hist_len-1) &
!                                                          - 27.0d0*vhist(m,j,hist_len-2) + vhist(m,j,hist_len-3))
!                    ! velocity part
!                    v1 = v1 + (     - zdhist(m,j,i,hist_len-0) + 8.0d0*zdhist(m,j,i,hist_len-1) &
!                               -8.0d0*zdhist(m,j,i,hist_len-3)       + zdhist(m,j,i,hist_len-4) ) * &
!                               (       - vhist(m,j,hist_len-0) + 7.0d0*vhist(m,j,hist_len-1) &
!                                 + 7.0d0*vhist(m,j,hist_len-2)       - vhist(m,j,hist_len-3))
!                end do
!            end do
!            pxi0(i) = (1.0d0/24.0d0)*f1*ifdtx
!            pxi1(i) = s1
!            pxi2(i) = (1.0d0/144.0d0)*v1*ifdtx
!            pxi3(i) = l1
!        end do
!
!        ! write(7894,*) fstep, ekinhist(hist_len-2)
!
!        ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
!        call abf_core_register_rawdata(cvhist(:,hist_len-2),pxi0,pxi1,pxi2,pxi3,micfhist(:,hist_len-2), &
!                               epothist(hist_len-2),ersthist(hist_len-2),ekinhist(hist_len-2))
!    end if
!
!    return
!
!end subroutine abf_core_force_3pV4
!
!!===============================================================================
!! Subroutine:  abf_core_force_3pV6
!! this is leap-frog ABF version, simplified algorithm
!! ICF from velocities + decomposition
!!===============================================================================
!
!subroutine abf_core_force_3pV6()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!
!    implicit none
!    integer                :: i,j,k,m
!    integer                :: ci,ki
!    real(PMFDP)            :: v1,s1,l1,f1
!    ! --------------------------------------------------------------------------
!
!! shift accuvalue history
!    do i=1,hist_len-1
!        cvhist(:,i)         = cvhist(:,i+1)
!        epothist(i)         = epothist(i+1)
!        ersthist(i)         = ersthist(i+1)
!        ekinhist(i)         = ekinhist(i+1)
!        ekinvvhist(i)         = ekinvvhist(i+1)
!        ekinlfhist(i)         = ekinlfhist(i+1)
!        vhist(:,:,i)        = vhist(:,:,i+1)
!        zdhist(:,:,:,i)     = zdhist(:,:,:,i+1)
!        micfhist(:,i)       = micfhist(:,i+1)
!    end do
!
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len) = CVContext%CVsValues(ci)
!    end do
!    vhist(:,:,hist_len)    = Vel(:,:)
!
!! shift ene
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!select case(ftds_ekin_src)
!    case(1)
!        ekinhist(hist_len-1)    = KinEneVV - fekinaverage    ! shifted by -dt
!    case(2)
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
!        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
!    case(3)
!        ! kinetic part
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
!        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
!        ! correction
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) &
!                                + 1.0d0/8.0d0*(epothist(hist_len)+ersthist(hist_len) + epothist(hist_len-2)+ersthist(hist_len-2)) &
!                                - 2.0d0/8.0d0*(epothist(hist_len-1)+ersthist(hist_len-1))
!    case(4)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-3/2dt and t-5/2dt
!                                  ! t-dt, t-2dt, t-3dt
!        ekinhist(hist_len-2)    = 0.5d0*(ekinlfhist(hist_len-2)+ekinlfhist(hist_len-3))  &
!                                - 1.0d0/8.0d0*(ekinvvhist(hist_len-1)-2.0d0*ekinvvhist(hist_len-2)+ekinvvhist(hist_len-3))
!    case(5)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-1/2dt t-3/2dt and t-5/2dt t-7/2dt
!        ekinhist(hist_len-2)    = 1.0d0/10.0d0*(-ekinlfhist(hist_len-1)+6.0d0*ekinlfhist(hist_len-2)&
!                                                +6.0d0*ekinlfhist(hist_len-3)-ekinlfhist(hist_len-4))
!    case(6)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-1/2dt t-3/2dt and t-5/2dt t-7/2dt
!        ekinhist(hist_len-2)    = 1.0d0/10.0d0*(-ekinlfhist(hist_len-1)+6.0d0*ekinlfhist(hist_len-2)&
!                                                +6.0d0*ekinlfhist(hist_len-3)-ekinlfhist(hist_len-4))
!        ekinhist(hist_len-2)    = 0.5d0*(ekinhist(hist_len-2) + ekinvvhist(hist_len-2))
!    case(7)
!        ekinhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!    case(8)
!        ekinvvhist(hist_len-1)  = KinEneVV - fekinaverage    ! shifted by t-dt
!        ekinlfhist(hist_len-1)  = KinEneLF - fekinaverage   ! shifted by t-dt/2
!                                  ! t-3/2dt and t-5/2dt
!                                  ! t-dt, t-2dt, t-3dt
!        ekinhist(hist_len-2)    = (4.0d0/6.0d0)*(ekinlfhist(hist_len-2)+ekinlfhist(hist_len-3))  &
!                                - (1.0d0/6.0d0)*(ekinvvhist(hist_len-1)+ekinvvhist(hist_len-3))
!    case(9)
!        ekinhist(hist_len-2)    = KinEneV4 - fekinaverage
!    case default
!        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_force_3pV!')
!    end select
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!
!    do i=1,NumOfABFCVs
!        do j=1,NumOfLAtoms
!            do m=1,3
!                v1 = 0.0d0
!                do k=1,NumOfABFCVs
!                    ki = ABFCVList(k)%cvindx
!                    v1 = v1 + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
!                end do
!                zdhist(m,j,i,hist_len) = v1
!            end do
!        end do
!    end do
!
!! apply force filters
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode in abf_core_force_3pV!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len) = la(:)
!
!    if( frecord ) then
!        ! record time progress of data
!        call abf_accu_add_data_record_lf(cvhist(:,hist_len),micfhist(:,hist_len), &
!                                         epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!     end if
!
!! ABF part
!    if( fstep .ge. hist_len ) then
!        do i=1,NumOfABFCVs
!            f1 = 0.0d0
!            s1 = 0.0d0 ! zero - partly in f1
!            l1 = 0.0d0 ! zero
!            v1 = 0.0d0
!            do j=1,NumOfLAtoms
!                do m=1,3
!                    ! force part
!                    f1 = f1 + zdhist(m,j,i,hist_len-3) * &
!                         (  54.0d0*vhist(m,j,hist_len-0) -3875.0d0*vhist(m,j,hist_len-1) + 97875.0d0*vhist(m,j,hist_len-2) &
!                          - 97875.0d0*vhist(m,j,hist_len-3) + 3875.0d0*vhist(m,j,hist_len-4) - 54.0d0*vhist(m,j,hist_len-5))
!                    ! velocity part
!                    v1 = v1 + (     - zdhist(m,j,i,hist_len-1) + 8.0d0*zdhist(m,j,i,hist_len-2) &
!                               -8.0d0*zdhist(m,j,i,hist_len-4)       + zdhist(m,j,i,hist_len-5) ) * &
!                               (       - vhist(m,j,hist_len-1) + 7.0d0*vhist(m,j,hist_len-2) &
!                                 + 7.0d0*vhist(m,j,hist_len-3)       - vhist(m,j,hist_len-4))
!                end do
!            end do
!            pxi0(i) = (1.0d0/86520.0d0)*f1*ifdtx
!            pxi1(i) = s1
!            pxi2(i) = (1.0d0/144.0d0)*v1*ifdtx
!            pxi3(i) = l1
!        end do
!
!        ! write(7894,*) fstep, ekinhist(hist_len-2)
!
!        ! no filter
!        ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
!        call abf_core_register_rawdata(cvhist(:,hist_len-3),pxi0,pxi1,pxi2,pxi3,micfhist(:,hist_len-3), &
!                               epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))
!    end if
!
!    return
!
!end subroutine abf_core_force_3pV6
!
!!===============================================================================
!! Subroutine:  abf_core_force_3pF
!! this is leap-frog ABF version, simplified algorithm
!! ICF from forces + decomposition
!!===============================================================================
!
!subroutine abf_core_force_3pF()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!
!    implicit none
!    integer                :: i,j,k,m
!    integer                :: ci,ki
!    real(PMFDP)            :: v1,v2,s1,l1,f1
!    ! --------------------------------------------------------------------------
!
!! shift accuvalue history
!    do i=1,hist_len-1
!        cvhist(:,i)     = cvhist(:,i+1)
!        epothist(i)     = epothist(i+1)
!        ersthist(i)     = ersthist(i+1)
!        ekinhist(i)     = ekinhist(i+1)
!        vhist(:,:,i)    = vhist(:,:,i+1)
!        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
!        micfhist(:,i)   = micfhist(:,i+1)
!        fhist(:,:,i)    = fhist(:,:,i+1)
!        shist(:,:,i)    = shist(:,:,i+1)
!        lhist(:,:,i)    = lhist(:,:,i+1)
!    end do
!
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len) = CVContext%CVsValues(ci)
!    end do
!    vhist(:,:,hist_len)    = Vel(:,:)
!
!! shift ene
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!    select case(ftds_ekin_src)
!    case(1)
!        ekinhist(hist_len-1)    = KinEneVV - fekinaverage    ! shifted by -dt
!    case(2)
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
!        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
!    case(3)
!        ! kinetic part
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this is completed
!        ekinhist(hist_len-0)    =                      + 0.5d0*(KinEneLF - fekinaverage)    ! in t-dt/2, this will be completed in the next step
!        ! correction
!        ekinhist(hist_len-1)    = ekinhist(hist_len-1) &
!                                + 1.0d0/8.0d0*(epothist(hist_len)+ersthist(hist_len) + epothist(hist_len-2)+ersthist(hist_len-2)) &
!                                - 2.0d0/8.0d0*(epothist(hist_len-1)+ersthist(hist_len-1))
!    case default
!        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_force_3pF!')
!    end select
!
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!
!    do i=1,NumOfABFCVs
!        do j=1,NumOfLAtoms
!            do m=1,3
!                v1 = 0.0d0
!                do k=1,NumOfABFCVs
!                    ki = ABFCVList(k)%cvindx
!                    v1 = v1 + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
!                end do
!                zdhist(m,j,i,hist_len) = v1
!            end do
!        end do
!    end do
!
!! apply force filters
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode in abf_core_force_3pF!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len) = la(:)
!
!    fhist(:,:,hist_len)  = Frc(:,:) ! Frc must include the added ICF
!
!    if( frecord ) then
!        ! record time progress of data
!        call abf_accu_add_data_record_lf(cvhist(:,hist_len),micfhist(:,hist_len), &
!                                         epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!     end if
!
!! ABF part
!    if( fstep .ge. hist_len ) then
!        do i=1,NumOfABFCVs
!            f1 = 0.0d0
!            s1 = 0.0d0
!            l1 = 0.0d0
!            v1 = 0.0d0
!            v2 = 0.0d0
!            do j=1,NumOfLAtoms
!                do m=1,3
!                    ! force part
!                    f1 = f1 + zdhist(m,j,i,hist_len-1) * fhist(m,j,hist_len-1) * MassInv(j)
!                    s1 = s1 + zdhist(m,j,i,hist_len-1) * shist(m,j,hist_len-1) * MassInv(j)
!                    l1 = l1 + zdhist(m,j,i,hist_len-1) * lhist(m,j,hist_len-1) * MassInv(j)
!                    ! velocity part
!                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
!                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
!                end do
!            end do
!            pxi0(i) = f1
!            pxi1(i) = s1
!            pxi2(i) = 0.5d0*(v1+v2)*ifdtx
!            pxi3(i) = l1
!        end do
!        ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,licf,bicf,epot,erst,ekin)
!        call abf_core_register_rawdata(cvhist(:,hist_len-1),pxi0,pxi1,pxi2,pxi3,micfhist(:,hist_len-1), &
!                               epothist(hist_len-1),ersthist(hist_len-1),ekinhist(hist_len-1))
!
!    end if
!
!    return
!
!end subroutine abf_core_force_3pF

!===============================================================================
! Subroutine:  abf_core_register_rawdata
!===============================================================================

subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    real(PMFDP),intent(in)  :: cvs(:)
    real(PMFDP),intent(in)  :: ficf(:)
    real(PMFDP),intent(in)  :: sicf(:)
    real(PMFDP),intent(in)  :: vicf(:)
    real(PMFDP),intent(in)  :: bicf(:)
    real(PMFDP),intent(in)  :: epot
    real(PMFDP),intent(in)  :: erst
    real(PMFDP),intent(in)  :: ekin
    ! --------------------------------------------
    real(PMFDP)             :: etot
    ! --------------------------------------------------------------------------

    ! total ABF force
    pxia(:) = ficf(:) + sicf(:) + vicf(:) ! adaptive correction
    ! bicf  ! current bias

    ! total energy
    etot = epot + erst + ekin

    ! add data to accumulator
    ! subroutine abf_accu_add_data_online(cvs,gfx,epot,erst,ekin,etot)
    call abf_accu_add_data_online(cvs,pxia,bicf,epot,erst,ekin,etot)

    if( fentropy .and. fentdecomp ) then
        ! subroutine abf_accu_add_data_entropy_decompose(cvs,fx,sx,vx,lx,bx,epot,erst,ekin)
        call abf_accu_add_data_entropy_decompose(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    end if

end subroutine abf_core_register_rawdata

!===============================================================================
! subroutine:  abf_core_calc_Zmat
!===============================================================================

subroutine abf_core_calc_Zmat(ctx)

    use pmf_utils
    use abf_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfABFCVs
            cj = ABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then

        fzinv(:,:)  = fz(:,:)

        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzinv,NumOfABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if

        call dgetri(NumOfABFCVs,fzinv,NumOfABFCVs,indx,vv,NumOfABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_core_calc_Zmat!')
        end if
    else
        fzinv(1,1)  = 1.0d0/fz(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat

!===============================================================================

end module abf_core
