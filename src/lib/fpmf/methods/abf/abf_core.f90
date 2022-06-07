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

    ! if we have enough data - run ABF
    select case(fmode)
        ! 3-points
        case(1)
            call abf_core_force_3pV1
        case(2)
            call abf_core_force_3pF
        ! 5-points
        case(3)
            call abf_core_force_5pV1
        ! symmetry forms
        case(4)
            call abf_core_force_2pX
        case(5)
            call abf_core_force_2pH
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    ! the remaining
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
        case(1)
            ! ignored
        case(2)
            ! get forces from SHAKE - abf_core_force_3pF
            shist(:,:,hist_len) = SHAKEFrc(:,:)
        case(3,4,5,6)
            ! ignored
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_shake!')
    end select

end subroutine abf_core_shake

!===============================================================================
! Subroutine:  abf_core_update_history
! apply ABF force and update history buffers
!===============================================================================

subroutine abf_core_update_history_I()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m,ci,ki
    real(PMFDP)            :: v, ekin
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
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
    epothist(hist_len)     = PotEne - fepotaverage
    ersthist(hist_len)     = PMFEne
    select case(ftds_ekin_src)
        ! KE from velocities at full-step interpolated from velocities at half-step
        case(0)
            ekinhist(hist_len-1)    = KinEneLF - fekinaverage   ! shifted by -1/2dt
            ekin                    = ekinhist(hist_len-1)
        case(1)
            ekinhist(hist_len-1)    = KinEneVV - fekinaverage   ! shifted by -dt
            ekin                    = ekinhist(hist_len-1)
        case(2)
            ekinhist(hist_len-2)    = KinEneV4 - fekinaverage
            ekin                    = ekinhist(hist_len-2)
        case(3)
            ekinhist(hist_len-3)    = KinEneV6 - fekinaverage
            ekin                    = ekinhist(hist_len-3)
        ! KE from interpolated KE at half-step
        case(4)
            ekinlfhist(hist_len)    = KinEneLF - fekinaverage   ! shifted by -1/2dt
            ekinhist(hist_len-1)    = 0.5d0*(ekinlfhist(hist_len-0) + ekinlfhist(hist_len-1))
            ekin                    = ekinhist(hist_len-1)
        case(5)
            ekinlfhist(hist_len)    = KinEneLF - fekinaverage
            ekinhist(hist_len-2)    = (1.0d0/16.0d0)*(      -ekinlfhist(hist_len-0)+9.0d0*ekinlfhist(hist_len-1) &
                                                      +9.0d0*ekinlfhist(hist_len-2)      -ekinlfhist(hist_len-3))
            ekin                    = ekinhist(hist_len-2)
        case(6)
            ekinlfhist(hist_len)    = KinEneLF - fekinaverage
            ekinhist(hist_len-3)    = (1.0d0/256.0d0)*(  +3.0d0*ekinlfhist(hist_len-0)  -25.0d0*ekinlfhist(hist_len-1) &
                                                       +150.0d0*ekinlfhist(hist_len-2) +150.0d0*ekinlfhist(hist_len-3) &
                                                        -25.0d0*ekinlfhist(hist_len-4)   +3.0d0*ekinlfhist(hist_len-5))
            ekin                    = ekinhist(hist_len-3)
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_update_history_I!')
    end select

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
                call pmf_utils_exit(PMF_OUT,1, &
                     '[ABF] Not implemented extrapolation/interpolation mode in abf_core_update_history_I!')
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
                                         epothist(hist_len),ersthist(hist_len),ekin)
    end if

    return

end subroutine abf_core_update_history_I

!===============================================================================
! Subroutine:  abf_core_update_history_II
! apply ABF force and update history buffers
!===============================================================================

subroutine abf_core_update_history_II()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,ci
    real(PMFDP)            :: ekin
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        ekinlfhist(i)       = ekinlfhist(i+1)
        vhist(:,:,i)        = vhist(:,:,i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)

! shift ene
    epothist(hist_len)     = PotEne - fepotaverage
    ersthist(hist_len)     = PMFEne
    select case(ftds_ekin_src)
        case(0)
            ekinhist(hist_len-1)    = KinEneLF - fekinaverage   ! shifted by -1/2dt
            ekin                    = ekinhist(hist_len-1)
        ! KE from velocities at full-step interpolated from velocities at half-step
        case(1)
            ekinhist(hist_len-1)    = KinEneVV - fekinaverage   ! shifted by -dt
            ekin                    = ekinhist(hist_len-1)
        case(2)
            ekinhist(hist_len-2)    = KinEneV4 - fekinaverage
            ekin                    = ekinhist(hist_len-2)
        case(3)
            ekinhist(hist_len-3)    = KinEneV6 - fekinaverage
            ekin                    = ekinhist(hist_len-3)
        ! KE from interpolated KE at half-step
        case(4)
            ekinlfhist(hist_len)    = KinEneLF - fekinaverage   ! shifted by -1/2dt
            ekinhist(hist_len-1)    = 0.5d0*(ekinlfhist(hist_len-0) + ekinlfhist(hist_len-1))
            ekin                    = ekinhist(hist_len-1)
        case(5)
            ekinlfhist(hist_len)    = KinEneLF - fekinaverage
            ekinhist(hist_len-2)    = (1.0d0/16.0d0)*(      -ekinlfhist(hist_len-0)+9.0d0*ekinlfhist(hist_len-1) &
                                                      +9.0d0*ekinlfhist(hist_len-2)      -ekinlfhist(hist_len-3))
            ekin                    = ekinhist(hist_len-2)
        case(6)
            ekinlfhist(hist_len)    = KinEneLF - fekinaverage
            ekinhist(hist_len-3)    = (1.0d0/256.0d0)*(  +3.0d0*ekinlfhist(hist_len-0)  -25.0d0*ekinlfhist(hist_len-1) &
                                                       +150.0d0*ekinlfhist(hist_len-2) +150.0d0*ekinlfhist(hist_len-3) &
                                                        -25.0d0*ekinlfhist(hist_len-4)   +3.0d0*ekinlfhist(hist_len-5))
            ekin                    = ekinhist(hist_len-3)
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_update_history_II!')
    end select

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

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
                call pmf_utils_exit(PMF_OUT,1, &
                     '[ABF] Not implemented extrapolation/interpolation mode in abf_core_update_history_II!')
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
                                         epothist(hist_len),ersthist(hist_len),ekin)
    end if

    return

end subroutine abf_core_update_history_II

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

    call abf_core_update_history_I

    if( fstep .le. 2*hist_len ) return

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

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV1

!===============================================================================
! Subroutine:  abf_core_force_3pF
! this is leap-frog ABF version, simplified algorithm
! ICF from forces + decomposition
!===============================================================================

subroutine abf_core_force_3pF()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: v1,v2,s1,f1
    ! --------------------------------------------------------------------------

    call abf_core_update_history_I

! shift accuvalue history
    do i=1,hist_len-1
        fhist(:,:,i)    = fhist(:,:,i+1)
        shist(:,:,i)    = shist(:,:,i+1)
    end do

    fhist(:,:,hist_len) = Frc(:,:)     ! to be compatible with forces derived from velocities, which also contain the bias
    ! shist is added later

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        s1 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * fhist(m,j,hist_len-3) * MassInv(j)
                s1 = s1 + zdhist(m,j,i,hist_len-3) * shist(m,j,hist_len-3) * MassInv(j)
                ! velocity part
                v1 = v1 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-3)) * vhist(m,j,hist_len-2)
                v2 = v2 + (zdhist(m,j,i,hist_len-3)-zdhist(m,j,i,hist_len-4)) * vhist(m,j,hist_len-3)
            end do
        end do
        pxif(i) = f1
        pxis(i) = s1
        pxiv(i) = 0.5d0*(v1+v2)*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pF

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

    call abf_core_update_history_I

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        v4 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * (- vhist(m,j,hist_len-1) + 27.0d0*vhist(m,j,hist_len-2) &
                                                      - 27.0d0*vhist(m,j,hist_len-3) + vhist(m,j,hist_len-4))
                ! velocity part
                v1 = v1 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                v2 = v2 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-3)) * vhist(m,j,hist_len-2)
                v3 = v3 + (zdhist(m,j,i,hist_len-3)-zdhist(m,j,i,hist_len-4)) * vhist(m,j,hist_len-3)
                v4 = v4 + (zdhist(m,j,i,hist_len-4)-zdhist(m,j,i,hist_len-5)) * vhist(m,j,hist_len-4)
            end do
        end do
        pxif(i) = (1.0d0/24.0d0)*f1*ifdtx
        pxiv(i) = (1.0d0/16.0d0)*(-v1 + 9.0d0*v2 + 9.0d0*v3 -v4)*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_5pV1

!===============================================================================
! Subroutine:  abf_core_force_2pX
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_2pX()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m,ki
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

    call abf_core_update_history_II

! shift accuvalue history
    do i=1,hist_len-1
        xvhist(:,i)     = xvhist(:,i+1)
        xvhist2(:,i)    = xvhist2(:,i+1)
    end do

    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ki = ABFCVList(i)%cvindx
                v = v + CVContext%CVsDrvs(m,j,ki)*Vel(m,j)
            end do
        end do
        pxia(i) = v
    end do

    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            ! FIXME - for more CVs the interpolation must be sooner than this term?
            v = v + fzinv(i,j)*pxia(j)
        end do
        xvhist2(i,hist_len) = v

        select case(abf_p2_vx)
        case(2)
            xvhist(i,hist_len-1) = 0.5d0*(xvhist2(i,hist_len-0)+xvhist2(i,hist_len-1))
        case(3)
            xvhist(i,hist_len-1) = (1.0d0/8.0d0)*(+3.0d0*xvhist2(i,hist_len-0)+6.0d0*xvhist2(i,hist_len-1) &
                                                        -xvhist2(i,hist_len-2))
        case(4)
            xvhist(i,hist_len-2) = (1.0d0/16.0d0)*(      -xvhist2(i,hist_len-0)+9.0d0*xvhist2(i,hist_len-1)&
                                                   +9.0d0*xvhist2(i,hist_len-2)      -xvhist2(i,hist_len-3))
        case(6)
            xvhist(i,hist_len-3) = (1.0d0/256.0d0)*(  +3.0d0*xvhist2(i,hist_len-0) -25.0d0*xvhist2(i,hist_len-1)&
                                                    +150.0d0*xvhist2(i,hist_len-2)+150.0d0*xvhist2(i,hist_len-3)&
                                                     -25.0d0*xvhist2(i,hist_len-4)  +3.0d0*xvhist2(i,hist_len-5))
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pX!')
        end select
    end do

    if( fstep .le. 10 ) return

    do i=1,NumOfABFCVs
        select case(abf_p2_px)
        case(3)
            pxif(i) = 0.5d0*(xvhist(i,hist_len-5) - xvhist(i,hist_len-7))*ifdtx
        case(4)
            pxif(i) = (1.0d0/6.0d0)*( +2.0d0*xvhist(i,hist_len-5) + 3.0d0*xvhist(i,hist_len-6) &
                                      -6.0d0*xvhist(i,hist_len-7)      + xvhist(i,hist_len-8))*ifdtx
        case(5)
            pxif(i) = (1.0d0/12.0d0)*(      -xvhist(i,hist_len-4) + 8.0d0*xvhist(i,hist_len-5) &
                                      -8.0d0*xvhist(i,hist_len-7)      + xvhist(i,hist_len-8))*ifdtx
        case(7)
            pxif(i) = (1.0d0/60.0d0)*(        xvhist(i,hist_len-3)  -9.0d0*xvhist(i,hist_len-4) &
                                      +45.0d0*xvhist(i,hist_len-5) -45.0d0*xvhist(i,hist_len-7) &
                                       +9.0d0*xvhist(i,hist_len-8)        -xvhist(i,hist_len-9))*ifdtx
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_force_2pX!')
        end select
        pxiv(i) = 0.0d0
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-6),pxif,pxis,pxiv,micfhist(:,hist_len-6), &
                           epothist(hist_len-6),ersthist(hist_len-6),ekinhist(hist_len-6))

end subroutine abf_core_force_2pX

!===============================================================================
! Subroutine:  abf_core_force_2pH
! this is leap-frog ABF version, simplified algorithm
! symmetrical variant: 2+3
!===============================================================================

subroutine abf_core_force_2pH()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m,ki
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

    call abf_core_update_history_II

! shift accuvalue history
    do i=1,hist_len-1
        xvhist(:,i)         = xvhist(:,i+1)
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
    end do
    fzinvhist(:,:,hist_len) = fzinv(:,:)

    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
        case(3)
            ! -1
            pxia(i) = 0.5d0*(cvhist(i,hist_len-0)-cvhist(i,hist_len-2))*ifdtx
        case(4)
            ! -1
            pxia(i) = (1.0d0/6.0d0)*( 2.0d0*cvhist(i,hist_len-0)+3.0d0*cvhist(i,hist_len-1)&
                                     -6.0d0*cvhist(i,hist_len-2)      +cvhist(i,hist_len-3))*ifdtx
        case(5)
            ! -2
            pxia(i) = (1.0d0/12.0d0)*(      -cvhist(i,hist_len-0)+8.0d0*cvhist(i,hist_len-1)&
                                      -8.0d0*cvhist(i,hist_len-3)      +cvhist(i,hist_len-4))*ifdtx
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pH!')
        end select
    end do

    !write(478958,*) fstep-1,cvhist(1,hist_len-0),pxia(1),cvhist(1,hist_len-0)+cvhist(1,hist_len-2)

    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
        case(3,4)
            ! -1
            v = 0.0d0
            do j=1,NumOfABFCVs
                v = v + fzinvhist(i,j,hist_len-1)*pxia(j)
            end do
            xvhist(i,hist_len-1) = v
        case(5)
            ! -2
            v = 0.0d0
            do j=1,NumOfABFCVs
                v = v + fzinvhist(i,j,hist_len-2)*pxia(j)
            end do
            xvhist(i,hist_len-2) = v
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pH - B!')
        end select
    end do

    if( fstep .le. 10 ) return

    do i=1,NumOfABFCVs
        select case(abf_p2_px)
        case(3)
            pxif(i) = 0.5d0*(xvhist(i,hist_len-5) - xvhist(i,hist_len-7))*ifdtx
        case(4)
            pxif(i) = (1.0d0/6.0d0)*( +2.0d0*xvhist(i,hist_len-5) + 3.0d0*xvhist(i,hist_len-6) &
                                      -6.0d0*xvhist(i,hist_len-7)      + xvhist(i,hist_len-8))*ifdtx
        case(5)
            pxif(i) = (1.0d0/12.0d0)*(      -xvhist(i,hist_len-4) + 8.0d0*xvhist(i,hist_len-5) &
                                      -8.0d0*xvhist(i,hist_len-7)      + xvhist(i,hist_len-8))*ifdtx
        case(7)
            pxif(i) = (1.0d0/60.0d0)*(        xvhist(i,hist_len-3)  -9.0d0*xvhist(i,hist_len-4) &
                                      +45.0d0*xvhist(i,hist_len-5) -45.0d0*xvhist(i,hist_len-7) &
                                       +9.0d0*xvhist(i,hist_len-8)        -xvhist(i,hist_len-9))*ifdtx
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_force_2pX!')
        end select
        pxiv(i) = 0.0d0
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-6),pxif,pxis,pxiv,micfhist(:,hist_len-6), &
                           epothist(hist_len-6),ersthist(hist_len-6),ekinhist(hist_len-6))

end subroutine abf_core_force_2pH

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
    real(PMFDP)             :: etot, ekin_scaled
    ! --------------------------------------------------------------------------

    ! total ABF force
    pxia(:) = ficf(:) + sicf(:) + vicf(:) ! adaptive correction
    ! bicf  ! current bias

    ! scale ekin
    ekin_scaled = ekin * ftds_ekin_scale

    ! total energy
    etot = epot + erst + ekin_scaled

    ! add data to accumulator
    ! subroutine abf_accu_add_data_online(cvs,gfx,epot,erst,ekin,etot)
    call abf_accu_add_data_online(cvs,pxia,bicf,epot,erst,ekin_scaled,etot)

    if( fentropy .and. fentdecomp ) then
        ! subroutine abf_accu_add_data_entropy_decompose(cvs,fx,sx,vx,lx,bx,epot,erst,ekin)
        call abf_accu_add_data_entropy_decompose(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin_scaled)
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
