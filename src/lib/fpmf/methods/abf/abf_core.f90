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
        ! simplified
        case(1)
            call abf_core_force_3pV1
        case(2)
            call abf_core_force_3pF
        case(3)
            call abf_core_force_5pV1
        ! special cases
        case(4)
            call abf_core_force_2pV
        case(5)
            call abf_core_force_2pX
        case(6)
            call abf_core_force_3pV2
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    ! the remaining
    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    if( fupdate_abf ) then
        call abf_client_exchange_data(.false.)
    end if

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
        case(1,3,4,5)
            ! ignored
        case(2)
            ! get forces from SHAKE - abf_core_force_3pF
            shist(:,:,hist_len) = SHAKEFrc(:,:)
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
    integer                :: i,j,ci
    real(PMFDP)            :: ekin,eus
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        ekinlfhist(i)       = ekinlfhist(i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
    end do

! shift ene
    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne

    ekinlfhist(hist_len)    = KinEneLF - fekinaverage   ! shifted by -1/2dt

    select case(ftds_ekin_src)
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
            ekinhist(hist_len-1)    = 0.5d0*(ekinlfhist(hist_len-0) + ekinlfhist(hist_len-1))
            ekin                    = ekinhist(hist_len-1)
        case(5)
            ekinhist(hist_len-2)    = (1.0d0/16.0d0)*(      -ekinlfhist(hist_len-0)+9.0d0*ekinlfhist(hist_len-1) &
                                                      +9.0d0*ekinlfhist(hist_len-2)      -ekinlfhist(hist_len-3))
            ekin                    = ekinhist(hist_len-2)
        case(6)
            ekinhist(hist_len-3)    = (1.0d0/256.0d0)*(  +3.0d0*ekinlfhist(hist_len-0)  -25.0d0*ekinlfhist(hist_len-1) &
                                                       +150.0d0*ekinlfhist(hist_len-2) +150.0d0*ekinlfhist(hist_len-3) &
                                                        -25.0d0*ekinlfhist(hist_len-4)   +3.0d0*ekinlfhist(hist_len-5))
            ekin                    = ekinhist(hist_len-3)
        case(7)
            ekinhist(hist_len-1)    = KinEneV3 - fekinaverage
            ekin                    = ekinhist(hist_len-1)
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_update_history!')
    end select

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

! apply force filters
    la(:) = 0.0d0
    if( fapply_abf ) then
        if( fusmode ) then
            ! US-ABF
            call abf_core_get_us_bias(cvhist(:,hist_len),la,eus)
            ! since ersthist(hist_len) is already updated
            ! we can add the bias into energy PMFEne
            PMFEne = PMFEne + eus
        else
            ! regular ABF
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
                         '[ABF] Not implemented extrapolation/interpolation mode in abf_core_update_history!')
            end select
        end if

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

end subroutine abf_core_update_history

!===============================================================================
! Subroutine:  abf_core_get_us_bias
! get bias from US restraints
!===============================================================================

subroutine abf_core_get_us_bias(values,gfx,bene)

    use abf_dat

    implicit none
    real(PMFDP)     :: values(:)
    real(PMFDP)     :: gfx(:)
    real(PMFDP)     :: bene
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    bene = 0.0

    do i=1,NumOfABFCVs

        ABFCVList(i)%deviation = ABFCVList(i)%cv%get_deviation(values(i),ABFCVList(i)%target_value)

        ABFCVList(i)%energy = 0.5d0*ABFCVList(i)%force_constant*ABFCVList(i)%deviation**2
        bene = bene + ABFCVList(i)%energy

        gfx(i) = - ABFCVList(i)%force_constant*ABFCVList(i)%deviation
    end do

end subroutine abf_core_get_us_bias

!===============================================================================
! Subroutine:  abf_core_update_zdhist
!===============================================================================

subroutine abf_core_update_zdhist()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,k,m,ki
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

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

end subroutine abf_core_update_zdhist

!===============================================================================
! Subroutine:  abf_core_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV1()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1,v2,epot,erst
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
    end do

    vhist(:,:,hist_len) = Vel(:,:) * ftds_vel_scale
    call abf_core_update_zdhist

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

    select case(ftds_epot_src)
        case(1)
            epot = epothist(hist_len-3)
            erst = ersthist(hist_len-3)
        case(2)
            epot = 0.5d0*(epothist(hist_len-2)+epothist(hist_len-4))
            erst = 0.5d0*(ersthist(hist_len-3)+ersthist(hist_len-4))
        case(3)
            epot = (1.0d0/3.0d0)*(epothist(hist_len-2)+epothist(hist_len-3)+epothist(hist_len-4))
            erst = (1.0d0/3.0d0)*(ersthist(hist_len-2)+ersthist(hist_len-3)+ersthist(hist_len-4))
        case(4)
            epot = (1.0d0/35.0d0)*( -3.0d0*epothist(hist_len-1)+12.0d0*epothist(hist_len-2) &
                                   +17.0d0*epothist(hist_len-3)+12.0d0*epothist(hist_len-4) &
                                    -3.0d0*epothist(hist_len-5))
            erst = (1.0d0/35.0d0)*( -3.0d0*ersthist(hist_len-1)+12.0d0*ersthist(hist_len-2) &
                                   +17.0d0*ersthist(hist_len-3)+12.0d0*ersthist(hist_len-4) &
                                    -3.0d0*ersthist(hist_len-5))
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_epot_src in abf_core_force_2pV!')
    end select

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epot,erst,ekinhist(hist_len-3))

end subroutine abf_core_force_3pV1

!===============================================================================
! Subroutine:  abf_core_force_3pV2
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV2()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m,k,ki
    real(PMFDP)            :: f1,v1,epot,erst,v2,dx1,dx2
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        zdhist(:,:,:,i)     = zdhist(:,:,:,i+1)
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
        vhist(:,:,i)        = vhist(:,:,i+1)
        cdrvhist(:,:,:,i)   = cdrvhist(:,:,:,i+1)
        fhist(:,:,i)        = fhist(:,:,i+1)
    end do

    fzinvhist(:,:,hist_len) = fzinv(:,:)
    vhist(:,:,hist_len)     = Vel(:,:) * ftds_vel_scale
    cdrvhist(:,:,:,hist_len)= CVContext%CVsDrvs(:,:,:)
    fhist(:,:,hist_len)     = Frc(:,:)
    call abf_core_update_zdhist

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs

        f1 = 0.0d0
        v1 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * fhist(m,j,hist_len-3) * MassInv(j) ! * (vhist(m,j,hist_len-2) - vhist(m,j,hist_len-3))
                ! velocity part
                ki = ABFCVList(i)%cvindx
                v1 = v1 + (      -cdrvhist(m,j,ki,hist_len-1) +8.0d0*cdrvhist(m,j,ki,hist_len-2)   &
                           -8.0d0*cdrvhist(m,j,ki,hist_len-4)       +cdrvhist(m,j,ki,hist_len-5)) &
                        * (vhist(m,j,hist_len-2) + vhist(m,j,hist_len-3))
            end do
        end do


        ! velocity part II
        v2 = 0.0d0
        do k=1,NumOfABFCVs
            dx1 = ABFCVList(k)%cv%get_deviation(cvhist(k,hist_len-5),cvhist(k,hist_len-1))
            dx2 = ABFCVList(k)%cv%get_deviation(cvhist(k,hist_len-2),cvhist(k,hist_len-4))
            v2 = v2 + (      -fzinvhist(k,i,hist_len-1) +8.0d0*fzinvhist(k,i,hist_len-2)   &
                       -8.0d0*fzinvhist(k,i,hist_len-4)       +fzinvhist(k,i,hist_len-5))  &
                    * (dx1+8.0d0*dx2)

!            v2 = v2 + (      -fzinvhist(k,i,hist_len-1) +8.0d0*fzinvhist(k,i,hist_len-2)   &
!                       -8.0d0*fzinvhist(k,i,hist_len-4)       +fzinvhist(k,i,hist_len-5))  &
!                    * (      -cvhist(k,hist_len-1) +8.0d0*cvhist(k,hist_len-2)  &
!                       -8.0d0*cvhist(k,hist_len-4)       +cvhist(k,hist_len-5))
        end do

        pxif(i) = f1
        pxis(i) = (1.0d0/12.0d0/2.0d0)*v1*ifdtx
        pxiv(i) = (1.0d0/12.0d0/12.0d0)*v2*ifdtx*ifdtx

       ! write(65893,*) fstep, pxiv(i), 0.5d0*(p1+p2)*ifdtx
    end do

    do i=1,NumOfABFCVs
        v1 = 0.0d0
        do j=1,NumOfABFCVs
            v1 = v1 + fzinvhist(j,i,hist_len-3)*pxis(j)
        end do
        pxiv(i) = pxiv(i) + v1
    end do
    pxis(:) = 0.0d0

    ! write(4789568,*) fstep, pxif(1), pxiv(1)

!    write(78948,*) fstep, fzinvhist(1,1,hist_len), cvhist(1,hist_len), &
!                   (1.0d0/12.0d0)*(-cvhist(1,hist_len-1) +8.0d0*cvhist(1,hist_len-2) &
!                    -8.0d0*cvhist(1,hist_len-4) +cvhist(1,hist_len-5))*ifdtx, &
!                   (1.0d0/2.0d0)*(cvhist(1,hist_len-2)-cvhist(1,hist_len-4))*ifdtx, &
!                   (1.0d0/12.0d0)*( -fzinvhist(1,1,hist_len-1) +8.0d0*fzinvhist(1,1,hist_len-2)   &
!                       -8.0d0*fzinvhist(1,1,hist_len-4)       +fzinvhist(1,1,hist_len-5))*ifdtx, &
!                      (1.0d0/2.0d0)*(fzinvhist(1,1,hist_len-2)-fzinvhist(1,1,hist_len-4))*ifdtx


    select case(ftds_epot_src)
        case(1)
            epot = epothist(hist_len-3)
            erst = ersthist(hist_len-3)
        case(2)
            epot = 0.5d0*(epothist(hist_len-2)+epothist(hist_len-4))
            erst = 0.5d0*(ersthist(hist_len-3)+ersthist(hist_len-4))
        case(3)
            epot = (1.0d0/3.0d0)*(epothist(hist_len-2)+epothist(hist_len-3)+epothist(hist_len-4))
            erst = (1.0d0/3.0d0)*(ersthist(hist_len-2)+ersthist(hist_len-3)+ersthist(hist_len-4))
        case(4)
            epot = (1.0d0/35.0d0)*( -3.0d0*epothist(hist_len-1)+12.0d0*epothist(hist_len-2) &
                                   +17.0d0*epothist(hist_len-3)+12.0d0*epothist(hist_len-4) &
                                    -3.0d0*epothist(hist_len-5))
            erst = (1.0d0/35.0d0)*( -3.0d0*ersthist(hist_len-1)+12.0d0*ersthist(hist_len-2) &
                                   +17.0d0*ersthist(hist_len-3)+12.0d0*ersthist(hist_len-4) &
                                    -3.0d0*ersthist(hist_len-5))
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_epot_src in abf_core_force_2pV!')
    end select

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epot,erst,ekinhist(hist_len-3))

end subroutine abf_core_force_3pV2

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

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        fhist(:,:,i)    = fhist(:,:,i+1)
        shist(:,:,i)    = shist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
    end do
    fhist(:,:,hist_len) = Frc(:,:)     ! to be compatible with forces derived from velocities, which also contain the bias
    ! shist is added later
    vhist(:,:,hist_len) = Vel(:,:) * ftds_vel_scale
    call abf_core_update_zdhist

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

    ! write(4789,*) fstep-3,pxif(1), pxis(1), pxiv(1)

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

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
    end do

    vhist(:,:,hist_len) = Vel(:,:) * ftds_vel_scale
    call abf_core_update_zdhist

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
    integer                :: i,j,cidx
    real(PMFDP)            :: v,dx1,dx2,dx3
    ! --------------------------------------------------------------------------

    call abf_core_update_history

! shift accuvalue history
    do i=1,hist_len-1
        xvhist(:,i)         = xvhist(:,i+1)
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
    end do
    fzinvhist(:,:,hist_len) = fzinv(:,:)

    ! this algorithm is not suitable for periodic CVs
    ! this is tested in abf_init_print_summary

    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
        case(3)
            ! -1
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-2))
            pxia(i) = 0.5d0*dx1*ifdtx
            ! pxia(i) = 0.5d0*(cvhist(i,hist_len-0)-cvhist(i,hist_len-2))*ifdtx

            cidx = -1
        case(5)
            ! -2
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-4),cvhist(i,hist_len-0))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-3))
            pxia(i) = (1.0d0/12.0d0)*(dx1+8.0d0*dx2)*ifdtx
            !pxia(i) = (1.0d0/12.0d0)*(      -cvhist(i,hist_len-0)+8.0d0*cvhist(i,hist_len-1)&
            !                          -8.0d0*cvhist(i,hist_len-3)      +cvhist(i,hist_len-4))*ifdtx
            cidx = -2
        case(7)
            ! -3
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-6))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-5),cvhist(i,hist_len-1))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-4))
            pxia(i) = (1.0d0/60.0d0)*(dx1+9.0d0*dx2+45.0d0*dx3)*ifdtx
            cidx = -3
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pX!')
        end select
    end do

    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            v = v + fzinvhist(i,j,hist_len+cidx)*pxia(j)
        end do
        xvhist(i,hist_len+cidx) = v
    end do

    if( fstep .le. 2*hist_len ) return

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
! Subroutine:  abf_core_force_2pV
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_2pV()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m,vidx
    real(PMFDP)            :: v,epot,erst
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        xvhist(:,i)     = xvhist(:,i+1)
    end do
    vhist(:,:,hist_len) = Vel(:,:) * ftds_vel_scale
    call abf_core_update_zdhist

    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
        case(2)
            ! -1
            vint(:,:) = 0.5d0*(vhist(:,:,hist_len-0)+vhist(:,:,hist_len-1))
            vidx = -1
        case(3)
            ! -1
            vint(:,:) = (1.0d0/8.0d0)*(+3.0d0*vhist(:,:,hist_len-0)+6.0d0*vhist(:,:,hist_len-1) &
                                             -vhist(:,:,hist_len-2))
            vidx = -1
        case(4)
            ! -2
            vint(:,:) = (1.0d0/16.0d0)*(      -vhist(:,:,hist_len-0)+9.0d0*vhist(:,:,hist_len-1)&
                                        +9.0d0*vhist(:,:,hist_len-2)      -vhist(:,:,hist_len-3))
            vidx = -2
        case(6)
            ! -3
            vint(:,:) = (1.0d0/256.0d0)*(  +3.0d0*vhist(:,:,hist_len-0) -25.0d0*vhist(:,:,hist_len-1)&
                                         +150.0d0*vhist(:,:,hist_len-2)+150.0d0*vhist(:,:,hist_len-3)&
                                          -25.0d0*vhist(:,:,hist_len-4)  +3.0d0*vhist(:,:,hist_len-5))
            vidx = -3
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pV!')
        end select
    end do

    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zdhist(m,j,i,hist_len+vidx)*vint(m,j)
            end do
        end do
        xvhist(i,hist_len+vidx) = v
    end do

    if( fstep .le. 2*hist_len ) return

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
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_force_2pV!')
        end select
        pxiv(i) = 0.0d0
    end do


    select case(ftds_epot_src)
        case(1)
            epot = epothist(hist_len-6)
            erst = ersthist(hist_len-6)
        case(2)
            epot = 0.5d0*(epothist(hist_len-5)+epothist(hist_len-7))
            erst = 0.5d0*(ersthist(hist_len-5)+ersthist(hist_len-7))
        case(3)
            epot = (1.0d0/3.0d0)*(epothist(hist_len-5)+epothist(hist_len-6)+epothist(hist_len-7))
            erst = (1.0d0/3.0d0)*(ersthist(hist_len-5)+ersthist(hist_len-6)+ersthist(hist_len-7))
        case(4)
            epot = (1.0d0/35.0d0)*( -3.0d0*epothist(hist_len-4)+12.0d0*epothist(hist_len-5) &
                                   +17.0d0*epothist(hist_len-6)+12.0d0*epothist(hist_len-7) &
                                    -3.0d0*epothist(hist_len-8))
            erst = (1.0d0/35.0d0)*( -3.0d0*ersthist(hist_len-4)+12.0d0*ersthist(hist_len-5) &
                                   +17.0d0*ersthist(hist_len-6)+12.0d0*ersthist(hist_len-7) &
                                    -3.0d0*ersthist(hist_len-8))
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_epot_src in abf_core_force_2pV!')
    end select

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-6),pxif,pxis,pxiv,micfhist(:,hist_len-6), &
                           epot,erst,ekinhist(hist_len-6))

end subroutine abf_core_force_2pV

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
