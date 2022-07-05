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

module abf_core_lf

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine abf_core_lf_main

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
            call abf_core_lf_force_3pV1
        case(2)
            call abf_core_lf_force_3pV2
        case(3)
            call abf_core_lf_force_3pV3
! not supported any more: 2, 3
        case(4)
            call abf_core_lf_force_2pV
        case(5)
            call abf_core_lf_force_2pX
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_lf_main!')
    end select

    ! the remaining
    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    if( fupdate_abf ) then
        call abf_client_exchange_data(.false.)
    end if

end subroutine abf_core_lf_main

!===============================================================================
! Subroutine:  abf_core_lf_force_3pV1
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_lf_force_3pV1()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

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

end subroutine abf_core_lf_force_3pV1

!===============================================================================
! Subroutine:  abf_core_lf_force_3pV2
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_lf_force_3pV2()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

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

end subroutine abf_core_lf_force_3pV2

!===============================================================================
! Subroutine:  abf_core_lf_force_3pV3
! this is leap-frog ABF version, simplified algorithm
! ICF from velocities + decomposition
!===============================================================================

subroutine abf_core_lf_force_3pV3()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

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

end subroutine abf_core_lf_force_3pV3

!===============================================================================
! Subroutine:  abf_core_lf_force_2pX
! this is leap-frog ABF version, numerical via CV values
!===============================================================================

subroutine abf_core_lf_force_2pX()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_accu
    use abf_core

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

end subroutine abf_core_lf_force_2pX

!===============================================================================
! Subroutine:  abf_core_lf_force_2pV
! this is leap-frog ABF version, numerical via velocities
!===============================================================================

subroutine abf_core_lf_force_2pV()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

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
            ! vint(:,:) = vhist(:,:,hist_len-1) ! 0.5d0*(vhist(:,:,hist_len-0)+vhist(:,:,hist_len-1))
            !vint(:,:) =  0.5d0*(vhist(:,:,hist_len-0)+vhist(:,:,hist_len-1))
            vint(:,:) =  vhist(:,:,hist_len-0)
            vidx =-1
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

    !if( fdebug ) then
        write(14789,*) fstep, pxia
    !end if

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

end subroutine abf_core_lf_force_2pV

!===============================================================================

end module abf_core_lf
