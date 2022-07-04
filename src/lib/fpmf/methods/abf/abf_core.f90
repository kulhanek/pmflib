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
            call abf_core_force_3pV2
        case(3)
            call abf_core_force_3pV3
! not supported any more: 2, 3
        case(4)
            call abf_core_force_2pV
        case(5)
            call abf_core_force_2pX
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
    real(PMFDP)            :: epot,erst,ekin,eus
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        epotrwhist(i)       = epotrwhist(i+1)
        erstrwhist(i)       = erstrwhist(i+1)
        ekinlfhist(i)       = ekinlfhist(i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
    end do

! raw data
    epotrwhist(hist_len)    = PotEne - fepotaverage
    erstrwhist(hist_len)    = PMFEne
    ekinlfhist(hist_len)    = KinEne%KinEneLF - fekinaverage   ! shifted by -1/2dt

! process EPOT
    select case(ftds_epot_src)
        case(1)
            epot = epotrwhist(hist_len)
            erst = erstrwhist(hist_len)
            epothist(hist_len) = epot
            ersthist(hist_len) = erst
        case(2)
            epot = 0.5d0*(epotrwhist(hist_len-0)+epotrwhist(hist_len-2))
            erst = 0.5d0*(erstrwhist(hist_len-0)+erstrwhist(hist_len-2))
            epothist(hist_len-1) = epot
            ersthist(hist_len-1) = erst
        case(3)
            epot = (1.0d0/3.0d0)*(epotrwhist(hist_len-0)+epotrwhist(hist_len-1)+epotrwhist(hist_len-2))
            erst = (1.0d0/3.0d0)*(erstrwhist(hist_len-0)+erstrwhist(hist_len-1)+erstrwhist(hist_len-2))
            epothist(hist_len-1) = epot
            ersthist(hist_len-1) = erst
        case(4)
            epot = (1.0d0/35.0d0)*( -3.0d0*epotrwhist(hist_len-0)+12.0d0*epotrwhist(hist_len-1) &
                                   +17.0d0*epotrwhist(hist_len-2)+12.0d0*epotrwhist(hist_len-3) &
                                    -3.0d0*epotrwhist(hist_len-4))
            erst = (1.0d0/35.0d0)*( -3.0d0*erstrwhist(hist_len-0)+12.0d0*erstrwhist(hist_len-1) &
                                   +17.0d0*erstrwhist(hist_len-2)+12.0d0*erstrwhist(hist_len-3) &
                                    -3.0d0*erstrwhist(hist_len-4))
            epothist(hist_len-2) = epot
            ersthist(hist_len-2) = erst
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_epot_src in abf_core_update_history!')
    end select

! process EKIN
    select case(ftds_ekin_src)
        ! KE from velocities at full-step interpolated from velocities at half-step
        case(1)
            ekinhist(hist_len-1)    = KinEne%KinEneVV - fekinaverage   ! shifted by -dt
            ekin                    = ekinhist(hist_len-1)
        case(2)
            ekinhist(hist_len-2)    = KinEne%KinEneV4 - fekinaverage
            ekin                    = ekinhist(hist_len-2)
        case(3)
            ekinhist(hist_len-3)    = KinEne%KinEneV6 - fekinaverage
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
            ekinhist(hist_len-1)    = KinEne%KinEneV3 - fekinaverage
            ekin                    = ekinhist(hist_len-1)
        case(8)
            ekinhist(hist_len-2)    = KinEne%KinEneV5 - fekinaverage
            ekin                    = ekinhist(hist_len-2)
        case(9)
            ekinhist(hist_len-2)    = KinEne%KinEneVV - fekinaverage
            ekin                    = ekinhist(hist_len-2)
        case(10)
            ekinhist(hist_len-0)    = KinEne%KinEneVV - fekinaverage
            ekin                    = ekinhist(hist_len-0)
        case(11)
            ekinhist(hist_len-1)    = KinEne%KinEneHA - fekinaverage
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
        ! FIXME
        do i=1,abfaccu%tot_cvs
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
                                         epot,erst,ekin)
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
! Subroutine:  abf_core_update_cvder
!===============================================================================

subroutine abf_core_update_cvder()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m,ki
    ! --------------------------------------------------------------------------

    do i=1,NumOfABFCVs
        ki = ABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            do m=1,3
                cvderhist(m,j,i,hist_len) = CVContext%CVsDrvs(m,j,ki)
            end do
        end do
    end do

end subroutine abf_core_update_cvder

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
    real(PMFDP)            :: f1,v1,v2
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

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV1

!===============================================================================
! Subroutine:  abf_core_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV2()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1
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
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * (1.0d0/24.0d0)*( -1.0d0*vhist(m,j,hist_len-1) +27.0d0*vhist(m,j,hist_len-2) &
                                                                     -27.0d0*vhist(m,j,hist_len-3)  +1.0d0*vhist(m,j,hist_len-4) )
                ! velocity part
                v1 = v1 + (1.0d0/12.0d0)*(-1.0d0*zdhist(m,j,i,hist_len-1) +8.0d0*zdhist(m,j,i,hist_len-2)   &
                                          -8.0d0*zdhist(m,j,i,hist_len-4) +1.0d0*zdhist(m,j,i,hist_len-5))  &
                        * (1.0d0/16.0d0)*(-1.0d0*vhist(m,j,hist_len-1) +9.0d0*vhist(m,j,hist_len-2)         &
                                          +9.0d0*vhist(m,j,hist_len-3) -1.0d0*vhist(m,j,hist_len-4))
            end do
        end do
        pxif(i) = f1*ifdtx
        pxiv(i) = v1*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV2

!===============================================================================
! Subroutine:  abf_core_force_3pV3
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_force_3pV3()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1,v1
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
        do j=1,NumOfLAtoms
            do m=1,3
                ! force part
                f1 = f1 + zdhist(m,j,i,hist_len-3) * &
                          (1.0d0/1920.0d0)*(   +9.0d0*vhist(m,j,hist_len-0)  -125.0d0*vhist(m,j,hist_len-1) &
                                            +2250.0d0*vhist(m,j,hist_len-2) -2250.0d0*vhist(m,j,hist_len-3) &
                                             +125.0d0*vhist(m,j,hist_len-4)    -9.0d0*vhist(m,j,hist_len-5) )
                ! velocity part
                v1 = v1 + (1.0d0/60.0d0)*( +1.0d0*zdhist(m,j,i,hist_len-0)  -9.0d0*zdhist(m,j,i,hist_len-1)     &
                                          +45.0d0*zdhist(m,j,i,hist_len-2) -45.0d0*zdhist(m,j,i,hist_len-4)     &
                                           +9.0d0*zdhist(m,j,i,hist_len-5)  -1.0d0*zdhist(m,j,i,hist_len-6))    &
                        * (1.0d0/256.0d0)*(  +3.0d0*vhist(m,j,hist_len-0)  -25.0d0*vhist(m,j,hist_len-1)            &
                                          +150.0d0*vhist(m,j,hist_len-2) +150.0d0*vhist(m,j,hist_len-3)             &
                                           -25.0d0*vhist(m,j,hist_len-4)   +3.0d0*vhist(m,j,hist_len-5))
            end do
        end do
        pxif(i) = f1*ifdtx
        pxiv(i) = v1*ifdtx
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-3),pxif,pxis,pxiv,micfhist(:,hist_len-3), &
                           epothist(hist_len-3),ersthist(hist_len-3),ekinhist(hist_len-3))

end subroutine abf_core_force_3pV3

!===============================================================================
! Subroutine:  abf_core_force_2pX
! this is leap-frog ABF version, numerical via CV values
!===============================================================================

subroutine abf_core_force_2pX()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_accu

    implicit none
    integer                :: i,j,cidx,fidx
    real(PMFDP)            :: v,dx1,dx2,dx3,dx4,dx5
    ! --------------------------------------------------------------------------

    call abf_core_update_history

! shift accuvalue history
    do i=1,hist_len-1
        xphist(:,i)         = xphist(:,i+1)
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
    end do
    fzinvhist(:,:,hist_len) = fzinv(:,:)

! calculate CV velocity, consider CV periodicity
    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
    ! central differences
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
        case(9)
            ! -4
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-8),cvhist(i,hist_len-0))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-7))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-6),cvhist(i,hist_len-2))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-5))
            pxia(i) = (1.0d0/840.0d0)*(3.0d0*dx1+32.0d0*dx2+168.0d0*dx3+672.0d0*dx4)*ifdtx
            cidx = -4
    ! backward differences
        case(14)
            ! -1
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            pxia(i) = (1.0d0/6.0d0)*(2.0d0*dx1+5.0d0*dx2-dx3)*ifdtx
            cidx = -1
        case(15)
            ! -1
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-4))
            pxia(i) = (1.0d0/12.0d0)*(3.0d0*dx1+13.0d0*dx2-5.0d0*dx3+dx4)*ifdtx
            cidx = -1
        case(16)
            ! -1
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-4))
            dx5 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-4),cvhist(i,hist_len-5))
            pxia(i) = (1.0d0/60.0d0)*(12.0d0*dx1+77.0d0*dx2-43.0d0*dx3+17.0d0*dx4-3.0d0*dx5)*ifdtx
            cidx = -1
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pX!')
        end select
    end do

    if( abf_clear_shaken_cvvel ) then
        do i=abfaccu%tot_cvs+1,NumOfABFCVs
            pxia(i) = 0.0d0 ! reset SHAKEn velocities
        end do
    end if

! get CV momenta
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            v = v + fzinvhist(i,j,hist_len+cidx)*pxia(j)
        end do
        xphist(i,hist_len+cidx) = v
    end do

    if( fstep .le. 2*hist_len ) return

! get CV forces
    do i=1,NumOfABFCVs
        select case(abf_p2_px)
    ! central differences
        case(3)
            pxif(i) = 0.5d0*(xphist(i,hist_len-7) - xphist(i,hist_len-9))*ifdtx
            fidx = -8
        case(5)
            pxif(i) = (1.0d0/12.0d0)*(      -xphist(i,hist_len-6) + 8.0d0*xphist(i,hist_len-7) &
                                      -8.0d0*xphist(i,hist_len-9)      + xphist(i,hist_len-10))*ifdtx
            fidx = -8
        case(7)
            pxif(i) = (1.0d0/60.0d0)*(        xphist(i,hist_len-5)  -9.0d0*xphist(i,hist_len-6) &
                                      +45.0d0*xphist(i,hist_len-7) -45.0d0*xphist(i,hist_len-9) &
                                       +9.0d0*xphist(i,hist_len-10)       -xphist(i,hist_len-11))*ifdtx
            fidx = -8
        case(9)
            pxif(i) = (1.0d0/840.0d0)*( -3.0d0*xphist(i,hist_len-4) +32.0d0*xphist(i,hist_len-5) &
                                      -168.0d0*xphist(i,hist_len-6)+672.0d0*xphist(i,hist_len-7) &
                                      -672.0d0*xphist(i,hist_len-9)+168.0d0*xphist(i,hist_len-10) &
                                       -32.0d0*xphist(i,hist_len-11) +3.0d0*xphist(i,hist_len-12))*ifdtx
            fidx = -8
    ! backward differences
        case(14)
            pxif(i) = (1.0d0/6.0d0)*(+2.0d0*xphist(i,hist_len-1) +3.0d0*xphist(i,hist_len-2) &
                                     -6.0d0*xphist(i,hist_len-3) +1.0d0*xphist(i,hist_len-4))*ifdtx
            fidx = -2
        case(15)
            pxif(i) = (1.0d0/12.0d0)*( +3.0d0*xphist(i,hist_len-1) +10.0d0*xphist(i,hist_len-2) &
                                      -18.0d0*xphist(i,hist_len-3)  +6.0d0*xphist(i,hist_len-4) &
                                       -1.0d0*xphist(i,hist_len-5))*ifdtx
            fidx = -2
        case(16)
            pxif(i) = (1.0d0/60.0d0)*( +12.0d0*xphist(i,hist_len-1) +65.0d0*xphist(i,hist_len-2) &
                                      -120.0d0*xphist(i,hist_len-3) +60.0d0*xphist(i,hist_len-4) &
                                       -20.0d0*xphist(i,hist_len-5)  +3.0d0*xphist(i,hist_len-6))*ifdtx
            fidx = -2
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_force_2pX!')
        end select
        pxiv(i) = 0.0d0
    end do

    pxis(:) = 0.0d0
    if( abf_use_shaken_icf ) then
        do i=abfaccu%tot_cvs+1,NumOfABFCVs
            pxis(1) = pxis(1) + pxif(i)
        end do
    end if

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len+fidx),pxif,pxis,pxiv,micfhist(:,hist_len+fidx), &
                           epothist(hist_len+fidx),ersthist(hist_len+fidx),ekinhist(hist_len+fidx))

end subroutine abf_core_force_2pX

!===============================================================================
! Subroutine:  abf_core_force_2pV
! this is leap-frog ABF version, numerical via velocities
!===============================================================================

subroutine abf_core_force_2pV()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils

    implicit none
    integer                :: i,j,m,vidx
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
        cvderhist(:,:,:,i)  = cvderhist(:,:,:,i+1)
        vhist(:,:,i)        = vhist(:,:,i+1)
        xphist(:,i)         = xphist(:,i+1)
    end do
    vhist(:,:,hist_len)     = Vel(:,:) * ftds_vel_scale
    fzinvhist(:,:,hist_len) = fzinv(:,:)
    call abf_core_update_cvder

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
                v = v + cvderhist(m,j,i,hist_len+vidx)*vint(m,j)
            end do
        end do
        pxia(i) = v
    end do

    if( fdebug ) then
        write(14789,*) fstep, pxia
    end if

    if( abf_clear_shaken_cvvel ) then
        do i=abfaccu%tot_cvs+1,NumOfABFCVs
            pxia(i) = 0.0d0 ! reset SHAKEn velocities
        end do
    end if

! get CV momenta
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            v = v + fzinvhist(i,j,hist_len+vidx)*pxia(j)
        end do
        xphist(i,hist_len+vidx) = v
    end do

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs
        select case(abf_p2_px)
        case(3)
            pxif(i) = 0.5d0*(xphist(i,hist_len-5) - xphist(i,hist_len-7))*ifdtx
        case(4)
            pxif(i) = (1.0d0/6.0d0)*( +2.0d0*xphist(i,hist_len-5) + 3.0d0*xphist(i,hist_len-6) &
                                      -6.0d0*xphist(i,hist_len-7)      + xphist(i,hist_len-8))*ifdtx
        case(5)
            pxif(i) = (1.0d0/12.0d0)*(      -xphist(i,hist_len-4) + 8.0d0*xphist(i,hist_len-5) &
                                      -8.0d0*xphist(i,hist_len-7)      + xphist(i,hist_len-8))*ifdtx
        case(7)
            pxif(i) = (1.0d0/60.0d0)*(        xphist(i,hist_len-3)  -9.0d0*xphist(i,hist_len-4) &
                                      +45.0d0*xphist(i,hist_len-5) -45.0d0*xphist(i,hist_len-7) &
                                       +9.0d0*xphist(i,hist_len-8)        -xphist(i,hist_len-9))*ifdtx
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_force_2pV!')
        end select
        pxiv(i) = 0.0d0
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-6),pxif,pxis,pxiv,micfhist(:,hist_len-6), &
                           epothist(hist_len-6),ersthist(hist_len-6),ekinhist(hist_len-6))

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
    real(PMFDP)         :: v,logdet
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

!    write(7894,*) fstep, fz

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
      ! call abf_init_gpr_invK(fzinv,logdet)
    else
        fzinv(1,1)  = 1.0d0/fz(1,1)
    end if

!    write(7895,*) fstep, fzinv
!
!    do i=1,NumOfABFCVs
!        do j=1,NumOfABFCVs
!            v = 0.0d0
!            do k=1,NumOfABFCVs
!                v = v + fz(i,k)*fzinv(k,j)
!            end do
!            write(7896,*) fstep, i, j, v
!        end do
!    end do

    return

end subroutine abf_core_calc_Zmat

!===============================================================================
! Subroutine:  abf_init_gpr_invK
!===============================================================================

subroutine abf_init_gpr_invK(mat,logdet)

    use pmf_dat
    use abf_dat
    use pmf_utils

    implicit none
    real(PMFDP)                 :: mat(:,:)
    real(PMFDP)                 :: logdet
    integer                     :: gpr_len,i,alloc_failed,info,lwork,irank,gpr_rank
    real(PMFDP)                 :: minv, maxv, a_rcond, r_sigma, gpr_rcond, gpr_rsigma
    real(PMFDP),allocatable     :: sig(:)
    real(PMFDP),allocatable     :: u(:,:)
    real(PMFDP),allocatable     :: vt(:,:)
    real(PMFDP),allocatable     :: sig_plus(:,:)
    real(PMFDP),allocatable     :: temp_mat(:,:)
    real(PMFDP),allocatable     :: twork(:)
    integer,allocatable         :: iwork(:)
    ! --------------------------------------------------------------------------

    gpr_len = NumOfABFCVs
    gpr_rank = -1
    gpr_rcond = 1e-6
    gpr_rsigma = 0.0d0

! allocate arrays
    allocate(                           &
            sig(gpr_len),               &
            u(gpr_len,gpr_len),         &
            vt(gpr_len,gpr_len),        &
            sig_plus(gpr_len,gpr_len),  &
            temp_mat(gpr_len,gpr_len),  &
            iwork(8*gpr_len),           &
            twork(1),                   &
            stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory I for GPR arrays in abf_init_gpr_invK!')
    end if

    sig(:)          = 0.0d0
    u(:,:)          = 0.0d0
    vt(:,:)         = 0.0d0
    sig_plus(:,:)   = 0.0d0
    temp_mat(:,:)   = 0.0d0


! query work size
    lwork = -1
    call dgesdd('A', gpr_len, gpr_len, mat, gpr_len, sig, u, gpr_len, vt, gpr_len, twork, lwork, iwork, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to run SVD I in abf_init_gpr_invK!')
    end if


    lwork = int(twork(1)) + 1
    if( lwork < 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to illegal work size in abf_init_gpr_invK!')
    end if

    deallocate(twork)
    allocate(   twork(lwork),   &
                stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to allocate memory II for GPR arrays in abf_init_gpr_invK!')
    end if

! run SVD
    call dgesdd('A', gpr_len, gpr_len, mat, gpr_len, sig, u, gpr_len, vt, gpr_len, twork, lwork, iwork, info)

    if( info .ne. 0 )  then
        call pmf_utils_exit(PMF_OUT,1, &
            '[ABF] Unable to run SVD II in abf_init_gpr_invK!')
    end if

    deallocate(iwork,twork)

    if( fdebug ) then
        open(unit=7812,file='gpr-kernel.isigma',status='UNKNOWN')
        do i=1,gpr_len
            write(7812,*) i, 1.0d0/sig(i)
        end do
        close(7812)
    end if

! invert singular numbers
    maxv = sig(1)
    minv = sig(1)
    do i=1,gpr_len
        if( maxv .lt. sig(i) ) then
            maxv = sig(i)
        end if
        if( minv .gt. sig(i) ) then
            minv = sig(i)
        end if
    end do

    a_rcond = minv/maxv

    irank = 0
    logdet = 0.0d0  ! this is logarithm of the determinant of original matrix
    r_sigma = 0.0d0
    if( gpr_rank .gt. 1 ) then
        do i=1,gpr_len
            if( i .le. gpr_rank ) then
               r_sigma = sig(i)
               sig_plus(i,i) = 1.0d0/sig(i)
               logdet        = logdet +  log(sig(i))
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    else if( gpr_rsigma .ne. 0.0d0 ) then
        do i=1,gpr_len
            if( sig(i) .gt. gpr_rsigma ) then
               r_sigma = sig(i)
               sig_plus(i,i) = 1.0d0/sig(i)
               logdet        = logdet +  log(sig(i))
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    else
        do i=1,gpr_len
            if( sig(i) .gt. gpr_rcond*maxv ) then
               r_sigma = sig(i)
               sig_plus(i,i) = 1.0d0/sig(i)
               logdet        = logdet +  log(sig(i))
               irank = irank + 1
            else
               sig_plus(i,i) = 0.0d0
            end if
        end do
    end if

! build pseudoinverse: V*sig_plus*UT
    call dgemm('N', 'T', gpr_len, gpr_len, gpr_len, 1.0d0, sig_plus, gpr_len, u, gpr_len, 0.0d0, temp_mat, gpr_len)
    call dgemm('T', 'N', gpr_len, gpr_len, gpr_len, 1.0d0, vt, gpr_len, temp_mat, gpr_len, 0.0d0, mat, gpr_len)

    deallocate(sig,sig_plus,u,vt,temp_mat)

 !   write(PMF_OUT,10) gpr_len,irank,r_sigma,a_rcond

 10 format('    Size = ',I5,'; Rank = ',I5,'; Rank sigma = ',E10.5,'; Real rcond = ',E10.5)

end subroutine abf_init_gpr_invK

!===============================================================================

end module abf_core
