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
        case(1)
            call abf_core_lf_force_2pX
        case(2)
            call abf_core_lf_force_2pV
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
    integer                :: i,j,cidx
    real(PMFDP)            :: v,dx1,dx2,dx3,dx4,dx5,dx6
    ! --------------------------------------------------------------------------

    call abf_core_update_history_force

! shift accuvalue history
    do i=1,hist_len-1
        xphist(:,i)         = xphist(:,i+1)
        vhist(:,:,i)        = vhist(:,:,i+1)        ! we do not need this, but for fenthalpy_der > 0, it must be here
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
    end do
    vhist(:,:,hist_len)     = Vel(:,:)              ! t-dt/2
    fzinvhist(:,:,hist_len) = fzinv(:,:)

! Finite Difference Coefficients Calculator
! https://web.media.mit.edu/~crtaylor/calculator.html

! calculate CV velocity, consider CV periodicity
    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
    ! central differences
        case(3)
            ! -1,0,1
            ! f_x = (-1*f[i-1]+0*f[i+0]+1*f[i+1])/(2*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-2))
            pxia(i) = 0.5d0*dx1*ifdtx
            cidx = -1
        case(5)
            ! -2,-1,0,1,2
            ! f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-4),cvhist(i,hist_len-0))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-3))
            pxia(i) = (1.0d0/12.0d0)*(dx1+8.0d0*dx2)*ifdtx
            cidx = -2
        case(7)
            ! -3,-2,-1,0,1,2,3
            ! f_x = (-1*f[i-3]+9*f[i-2]-45*f[i-1]+0*f[i+0]+45*f[i+1]-9*f[i+2]+1*f[i+3])/(60*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-6))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-5),cvhist(i,hist_len-1))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-4))
            pxia(i) = (1.0d0/60.0d0)*(dx1+9.0d0*dx2+45.0d0*dx3)*ifdtx
            cidx = -3
        case(9)
            ! -4,-3,-2,-1,0,1,2,3,4
            ! f_x = (3*f[i-4]-32*f[i-3]+168*f[i-2]-672*f[i-1]+0*f[i+0]+672*f[i+1]-168*f[i+2]+32*f[i+3]-3*f[i+4])/(840*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-8),cvhist(i,hist_len-0))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-7))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-6),cvhist(i,hist_len-2))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-5))
            pxia(i) = (1.0d0/840.0d0)*(3.0d0*dx1+32.0d0*dx2+168.0d0*dx3+672.0d0*dx4)*ifdtx
            cidx = -4
        case(11)
            ! -5,-4,-3,-2,-1,0,1,2,3,4,5
            ! f_x = (-2*f[i-5]+25*f[i-4]-150*f[i-3]+600*f[i-2]-2100*f[i-1]+0*f[i+0]+2100*f[i+1]-600*f[i+2]+150*f[i+3]-25*f[i+4]+2*f[i+5])/(2520*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-10))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-9),cvhist(i,hist_len-1))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-8))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-7),cvhist(i,hist_len-3))
            dx5 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-4),cvhist(i,hist_len-6))
            pxia(i) = (1.0d0/2520.0d0)*(2.0d0*dx1+25.0d0*dx2+150.0d0*dx3+600.0d0*dx4+2100.0d0*dx5)*ifdtx
            cidx = -5
    ! backward differences
        case(13)
            ! -2,-1,0
            ! f_x = (1*f[i-2]-4*f[i-1]+3*f[i+0])/(2*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            pxia(i) = (1.0d0/2.0d0)*(3.0d0*dx1-1.0d0*dx2)*ifdtx
            cidx = 0
        case(14)
            ! -3,-2,-1,0
            ! f_x = (-2*f[i-3]+9*f[i-2]-18*f[i-1]+11*f[i+0])/(6*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            pxia(i) = (1.0d0/6.0d0)*(11.0d0*dx1-7.0d0*dx2+2.0d0*dx3)*ifdtx
            cidx = 0
        case(15)
            ! -4,-3,-2,-1,0
            ! f_x = (3*f[i-4]-16*f[i-3]+36*f[i-2]-48*f[i-1]+25*f[i+0])/(12*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-4))
            pxia(i) = (1.0d0/12.0d0)*(25.0d0*dx1-23.0d0*dx2+13.0d0*dx3-3.0d0*dx4)*ifdtx
            cidx = 0
        case(16)
            ! -5,-4,-3,-2,-1,0
            ! f_x = (-12*f[i-5]+75*f[i-4]-200*f[i-3]+300*f[i-2]-300*f[i-1]+137*f[i+0])/(60*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-4))
            dx5 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-4),cvhist(i,hist_len-5))
            pxia(i) = (1.0d0/60.0d0)*(137.0d0*dx1-163.0d0*dx2+137.0d0*dx3-63.0d0*dx4+12.0d0*dx5)*ifdtx
            cidx = 0
        case(17)
            ! -6,-5,-4,-3,-2,-1,0
            ! f_x = (10*f[i-6]-72*f[i-5]+225*f[i-4]-400*f[i-3]+450*f[i-2]-360*f[i-1]+147*f[i+0])/(60*1.0*h**1)
            dx1 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-0),cvhist(i,hist_len-1))
            dx2 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-1),cvhist(i,hist_len-2))
            dx3 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-2),cvhist(i,hist_len-3))
            dx4 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-3),cvhist(i,hist_len-4))
            dx5 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-4),cvhist(i,hist_len-5))
            dx6 = ABFCVList(i)%cv%get_deviation(cvhist(i,hist_len-5),cvhist(i,hist_len-6))
            pxia(i) = (1.0d0/60.0d0)*(147.0d0*dx1-213.0d0*dx2+237.0d0*dx3-163.0d0*dx4+62.0d0*dx5-10.0d0*dx6)*ifdtx
            cidx = 0
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_force_2pX!')
        end select
    end do

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
            ! -1,0,1
            ! f_x = (-1*f[i-1]+0*f[i+0]+1*f[i+1])/(2*1.0*h**1)
            pxif(i) = 0.5d0*(xphist(i,hist_len-9) - xphist(i,hist_len-11))*ifdtx
            hist_fidx = -10
        case(5)
            ! -2,-1,0,1,2
            ! f_x = (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*h**1)
            pxif(i) = (1.0d0/12.0d0)*( -       xphist(i,hist_len-8)  + 8.0d0*xphist(i,hist_len-9) &
                                       - 8.0d0*xphist(i,hist_len-11) +       xphist(i,hist_len-12))*ifdtx
            hist_fidx = -10
        case(7)
            ! -3,-2,-1,0,1,2,3
            ! f_x = (-1*f[i-3]+9*f[i-2]-45*f[i-1]+0*f[i+0]+45*f[i+1]-9*f[i+2]+1*f[i+3])/(60*1.0*h**1)
            pxif(i) = (1.0d0/60.0d0)*(          xphist(i,hist_len-7)  -  9.0d0*xphist(i,hist_len-8)  &
                                       + 45.0d0*xphist(i,hist_len-9)  - 45.0d0*xphist(i,hist_len-11) &
                                       +  9.0d0*xphist(i,hist_len-12) -        xphist(i,hist_len-13))*ifdtx
            hist_fidx = -10
        case(9)
            ! -4,-3,-2,-1,0,1,2,3,4
            ! f_x = (3*f[i-4]-32*f[i-3]+168*f[i-2]-672*f[i-1]+0*f[i+0]+672*f[i+1]-168*f[i+2]+32*f[i+3]-3*f[i+4])/(840*1.0*h**1)
            pxif(i) = (1.0d0/840.0d0)*( -3.0d0*xphist(i,hist_len-6)  +  32.0d0*xphist(i,hist_len-7)  &
                                      -168.0d0*xphist(i,hist_len-8)  + 672.0d0*xphist(i,hist_len-9)  &
                                      -672.0d0*xphist(i,hist_len-11) + 168.0d0*xphist(i,hist_len-12) &
                                       -32.0d0*xphist(i,hist_len-13) +   3.0d0*xphist(i,hist_len-14))*ifdtx
            hist_fidx = -10
        case(11)
            ! -5,-4,-3,-2,-1,0,1,2,3,4,5
            ! f_x = (-2*f[i-5]+25*f[i-4]-150*f[i-3]+600*f[i-2]-2100*f[i-1]+0*f[i+0]+2100*f[i+1]-600*f[i+2]+150*f[i+3]-25*f[i+4]+2*f[i+5])/(2520*1.0*h**1)
            pxif(i) = (1.0d0/2520.0d0)*( +    2.0d0*xphist(i,hist_len-5)  -   25.0d0*xphist(i,hist_len-6)  &
                                         +  150.0d0*xphist(i,hist_len-7)  -  600.0d0*xphist(i,hist_len-8)  &
                                         + 2100.0d0*xphist(i,hist_len-9)  - 2100.0d0*xphist(i,hist_len-11) &
                                         +  600.0d0*xphist(i,hist_len-12) -  150.0d0*xphist(i,hist_len-13) &
                                         +   25.0d0*xphist(i,hist_len-14) -    2.0d0*xphist(i,hist_len-15))*ifdtx
            hist_fidx = -10
    ! backward differences
    ! offset -2 to get fenthlapy_der
        case(13)
            ! -2,-1,0
            ! f_x = (1*f[i-2]-4*f[i-1]+3*f[i+0])/(2*1.0*h**1)
            pxif(i) = (1.0d0/2.0d0) *(+3.0d0*xphist(i,hist_len-0-2) -4.0d0*xphist(i,hist_len-1-2) &
                                      +1.0d0*xphist(i,hist_len-2-2))*ifdtx
            hist_fidx = -2
        case(14)
            ! -3,-2,-1,0
            ! f_x = (-2*f[i-3]+9*f[i-2]-18*f[i-1]+11*f[i+0])/(6*1.0*h**1)
            pxif(i) = (1.0d0/6.0d0) *(+11.0d0*xphist(i,hist_len-0-2) -18.0d0*xphist(i,hist_len-1-2) &
                                       +9.0d0*xphist(i,hist_len-2-2)  -2.0d0*xphist(i,hist_len-3-2))*ifdtx
            hist_fidx = -2
        case(15)
            ! -4,-3,-2,-1,0
            ! f_x = (3*f[i-4]-16*f[i-3]+36*f[i-2]-48*f[i-1]+25*f[i+0])/(12*1.0*h**1)
            pxif(i) = (1.0d0/12.0d0)*(+25.0d0*xphist(i,hist_len-0-2) -48.0d0*xphist(i,hist_len-1-2) &
                                      +36.0d0*xphist(i,hist_len-2-2) -16.0d0*xphist(i,hist_len-3-2) &
                                       +3.0d0*xphist(i,hist_len-4-2))*ifdtx
            hist_fidx = -2
        case(16)
            ! -5,-4,-3,-2,-1,0
            ! f_x = (-12*f[i-5]+75*f[i-4]-200*f[i-3]+300*f[i-2]-300*f[i-1]+137*f[i+0])/(60*1.0*h**1)
            pxif(i) = (1.0d0/60.0d0)*(+137.0d0*xphist(i,hist_len-0-2) -300.0d0*xphist(i,hist_len-1-2) &
                                      +300.0d0*xphist(i,hist_len-2-2) -200.0d0*xphist(i,hist_len-3-2) &
                                       +75.0d0*xphist(i,hist_len-4-2)  -12.0d0*xphist(i,hist_len-5-2))*ifdtx
            hist_fidx = -2
        case(17)
            ! -6,-5,-4,-3,-2,-1,0
            ! f_x = (10*f[i-6]-72*f[i-5]+225*f[i-4]-400*f[i-3]+450*f[i-2]-360*f[i-1]+147*f[i+0])/(60*1.0*h**1)
            pxif(i) = (1.0d0/60.0d0)*(+147.0d0*xphist(i,hist_len-0-2) -360.0d0*xphist(i,hist_len-1-2) &
                                      +450.0d0*xphist(i,hist_len-2-2) -400.0d0*xphist(i,hist_len-3-2) &
                                      +225.0d0*xphist(i,hist_len-4-2)  -72.0d0*xphist(i,hist_len-5-2) &
                                       +10.0d0*xphist(i,hist_len-6-2))*ifdtx
            hist_fidx = -2
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_force_2pX!')
        end select
    end do

    icfhist(:,hist_len+hist_fidx) = pxif(:)

    if( mod(fstep,ficfsample) .eq. 0 ) then

!        write(1486,*) fstep

        ! register data
        call abf_accu_add_data_online(cvhist(:,hist_len+hist_fidx),icfhist(:,hist_len+hist_fidx),&
                                      micfhist(:,hist_len+hist_fidx),fdetzhist(hist_len+hist_fidx),fziihist(:,hist_len+hist_fidx))
    end if

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
    use abf_accu

    implicit none
    integer                :: i,j,m,vidx
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history_force

    ! shift history buffers
    do i=1,hist_len-1
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
        cvderhist(:,:,:,i)  = cvderhist(:,:,:,i+1)
        vhist(:,:,i)        = vhist(:,:,i+1)
        xphist(:,i)         = xphist(:,i+1)
    end do
    vhist(:,:,hist_len)     = Vel(:,:)             ! t-dt/2
    fzinvhist(:,:,hist_len) = fzinv(:,:)
    call abf_core_update_cvder

    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
        case(2)
            ! -1
            vint(:,:) =  0.5d0*(vhist(:,:,hist_len-0)+vhist(:,:,hist_len-1))
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

! get CV momenta
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            v = v + fzinvhist(i,j,hist_len+vidx)*pxia(j)
        end do
        xphist(i,hist_len+vidx) = v
    end do

    if( fstep .le. 2*hist_len ) return

    hist_fidx = -6
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
    end do

    icfhist(:,hist_len+hist_fidx) = pxif(:)

    if( mod(fstep,ficfsample) .eq. 0 ) then
        ! register data
        call abf_accu_add_data_online(cvhist(:,hist_len+hist_fidx),icfhist(:,hist_len+hist_fidx),&
                                      micfhist(:,hist_len+hist_fidx),fdetzhist(hist_len+hist_fidx),fziihist(:,hist_len+hist_fidx))
    end if

end subroutine abf_core_lf_force_2pV

!===============================================================================
! Subroutine:  abf_core_lf_get_icfp
! this is leap-frog ABF version
!===============================================================================

subroutine abf_core_lf_get_icfp()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use abf_core

    implicit none
    integer                :: i,j,m
    real(PMFDP)            :: f1
    ! --------------------------------------------------------------------------

    ! shift history buffers
    do i=1,hist_len-1
        fhist(:,:,i)    = fhist(:,:,i+1)
        icfphist(:,i)   = icfphist(:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
    end do
                                        ! at this moment, Frc contains ABF bias
    fhist(:,:,hist_len) = Frc(:,:)      ! to be compatible with forces derived from velocities, which also contain the bias
    call abf_core_update_zdhist

    icfphist(:,hist_len) = 0.0d0

    if( fstep .le. 2*hist_len ) return

    if( fenthalpy_der .eq. 1 ) then
        do i=1,NumOfABFCVs
            f1 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    !                  t                        t
                    f1 = f1 + zdhist(m,j,i,hist_len) * fhist(m,j,hist_len) * MassInv(j)
                end do
            end do
            ! remove bias
            icfphist(i,hist_len) = f1 - micfhist(i,hist_len)
        end do
    else if ( fenthalpy_der .eq. 2 ) then
        do i=1,NumOfABFCVs
            f1 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    !                  t-dt                    t-dt/2               t-2*dt/2
                    f1 = f1 + zdhist(m,j,i,hist_len-1) * (vhist(m,j,hist_len) - vhist(m,j,hist_len-1))
                end do
            end do
            ! remove bias
            icfphist(i,hist_len-1) = f1*ifdtx - micfhist(i,hist_len-1)
        end do
    else if ( fenthalpy_der .eq. 3 ) then
        do i=1,NumOfABFCVs
            f1 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f1 = f1 + zdhist(m,j,i,hist_len-2) * (       - vhist(m,j,hist_len-0) + 27.0d0*vhist(m,j,hist_len-1) &
                                                          - 27.0d0*vhist(m,j,hist_len-2)        + vhist(m,j,hist_len-3))
                end do
            end do
            ! remove bias
            icfphist(i,hist_len-2) = (1.0d0/24.0d0)*f1*ifdtx - micfhist(i,hist_len-2)
        end do
    end if

end subroutine abf_core_lf_get_icfp

!===============================================================================
! Subroutine:  abf_core_lf_register_ekin
! this is leap-frog ABF version
!===============================================================================

subroutine abf_core_lf_register_ekin()

    use pmf_dat
    use abf_dat
    use abf_accu
    use pmf_utils

    implicit none
    integer     :: i
    real(PMFDP) :: lepot, lerst, lekin
    ! --------------------------------------------------------------------------

    if( .not. (fenthalpy .or. fentropy) ) return

    ! kinetic energy at this point is in the same time as potential energy

! shift accuvalue history
    do i=1,hist_len-1
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        ekinlfhist(i)       = ekinlfhist(i+1)
        enevalidhist(i)     = enevalidhist(i+1)
        epvhist(i)          = epvhist(i+1)
        volhist(i)          = volhist(i+1)
    end do

! raw data
    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = RstEne
    ekinlfhist(hist_len)    = KinEne%KinEneLF - fekinaverage   ! shifted by +1/2dt

    select case(finclude_pv)
        case(0)
            epvhist(hist_len)       = 0.0d0
        case(1)
            epvhist(hist_len)       = p0VEne
        case(2)
            epvhist(hist_len)       = pVEne
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unknown pV mode in abf_core_lf_register_ekin!')
    end select

    enevalidhist(hist_len)  = KinEne%Valid
    volhist(hist_len)       = fbox_volume

    ! get ICF-P
    call abf_core_lf_get_icfp

!    write(14789,*) fstep, pVEne, fbox_volume

! process EKIN
    select case(ftds_ekin_src)
        case(1)
            ekinhist(hist_len)      = KinEne%KinEneVV - fekinaverage
        case(2)
            ekinhist(hist_len)      = 0.5d0*(ekinlfhist(hist_len-0) + ekinlfhist(hist_len-1))
        case(3)
            ekinhist(hist_len)      = KinEne%KinEneHA - fekinaverage
        case(4)
            ekinhist(hist_len)      = ekinlfhist(hist_len) ! shifted by +1/2dt
        case(5)
            call abf_core_lf_register_ekin_v5
        case(6)
            call abf_core_lf_register_ekin_v7
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_lf_register_ekin!')
    end select

    if( fstep .le. 2*hist_len ) return

    if( enevalidhist(hist_len) ) fene_step = fene_step + 1

    if( .not. ( (mod(fene_step,fenesample) .eq. 0) .and. enevalidhist(hist_len+hist_fidx) ) ) return

    lepot = epothist(hist_len+hist_fidx)
    lerst = ersthist(hist_len+hist_fidx)
    lekin = ekinhist(hist_len+hist_fidx)

! https://web.media.mit.edu/~crtaylor/calculator.html

    select case(fepotsmooth)
        case(0)
            ! nothing to be here
        case(5)
            lepot = (  - 1.0d0 * epothist(hist_len+hist_fidx-2) &
                       + 4.0d0 * epothist(hist_len+hist_fidx-1) &
                       + 4.0d0 * epothist(hist_len+hist_fidx+1) &
                       - 1.0d0 * epothist(hist_len+hist_fidx+2) &
                    ) / 6.0d0
        case(7)
            lepot = (  + 1.0d0 * epothist(hist_len+hist_fidx-3) &
                       - 6.0d0 * epothist(hist_len+hist_fidx-2) &
                       +15.0d0 * epothist(hist_len+hist_fidx-1) &
                       +15.0d0 * epothist(hist_len+hist_fidx+1) &
                       - 6.0d0 * epothist(hist_len+hist_fidx+2) &
                       + 1.0d0 * epothist(hist_len+hist_fidx+3) &
                    ) / 20.0d0
        case(9)
            lepot = (  - 1.0d0 * epothist(hist_len+hist_fidx-4) &
                       + 8.0d0 * epothist(hist_len+hist_fidx-3) &
                       -28.0d0 * epothist(hist_len+hist_fidx-2) &
                       +56.0d0 * epothist(hist_len+hist_fidx-1) &
                       +56.0d0 * epothist(hist_len+hist_fidx+1) &
                       -28.0d0 * epothist(hist_len+hist_fidx+2) &
                       + 8.0d0 * epothist(hist_len+hist_fidx+3) &
                       - 1.0d0 * epothist(hist_len+hist_fidx+4) &
                    ) / 70.0d0
        case(11)
            ! -5,-4,-3,-2,-1,1,2,3,4,5
            ! f_ = (1*f[i-5]-10*f[i-4]+45*f[i-3]-120*f[i-2]+210*f[i-1]+210*f[i+1]-120*f[i+2]+45*f[i+3]-10*f[i+4]+1*f[i+5])/(252*1.0*h**0)
            lepot = (  +   1.0d0 * epothist(hist_len+hist_fidx-5) &
                       -  10.0d0 * epothist(hist_len+hist_fidx-4) &
                       +  45.0d0 * epothist(hist_len+hist_fidx-3) &
                       - 120.0d0 * epothist(hist_len+hist_fidx-2) &
                       + 210.0d0 * epothist(hist_len+hist_fidx-1) &
                       + 210.0d0 * epothist(hist_len+hist_fidx+1) &
                       - 120.0d0 * epothist(hist_len+hist_fidx+2) &
                       +  45.0d0 * epothist(hist_len+hist_fidx+3) &
                       -  10.0d0 * epothist(hist_len+hist_fidx+4) &
                       +   1.0d0 * epothist(hist_len+hist_fidx+5) &
                    ) / 252.0d0

        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported fepotsmooth mode in abf_core_lf_register_ekin!')
    end select


  !  write(1487,*) fstep

    ! register data
    call abf_accu_add_data_energy(cvhist(:,hist_len+hist_fidx), &
                  icfhist(:,hist_len+hist_fidx), micfhist(:,hist_len+hist_fidx), icfphist(:,hist_len+hist_fidx), &
                  lepot, lerst, lekin, &
                  epvhist(hist_len+hist_fidx),volhist(hist_len+hist_fidx))

end subroutine abf_core_lf_register_ekin

!===============================================================================
! Subroutine:  abf_core_lf_register_ekin_v5
! this is leap-frog ABF version
!===============================================================================

subroutine abf_core_lf_register_ekin_v5()

    use pmf_dat
    use abf_dat
    use abf_accu
    use pmf_utils

    implicit none
    integer     :: i, m
    real(PMFDP) :: ekinvv,ekinv5,v2,v5
    ! --------------------------------------------------------------------------

    ! old VV Ekin
    ekinvv = 0.0d0
    do i=1,NumOfLAtoms
        v2 = 0.0d0
        do m=1,3
            v2 = v2 + vhist(m,i,hist_len-2)**2
        end do
        ekinvv = ekinvv + Mass(i)*v2
    end do
    ekinvv = 0.5d0*ekinvv

    ! new V5 Ekin
    ekinv5 = 0.0d0
    do i=1,NumOfLAtoms
        v2 = 0.0d0
        do m=1,3
            !            t-dt/2
            v5 = -       vhist(m,i,hist_len-0) + 9.0d0*vhist(m,i,hist_len-1) &
                 + 9.0d0*vhist(m,i,hist_len-2) -       vhist(m,i,hist_len-3)
            v5 = (1.0d0 / 16.0d0) * v5
            v2 = v2 + v5**2
        end do
        ekinv5 = ekinv5 + Mass(i)*v2
    end do
    ekinv5 = 0.5d0*ekinv5

  !  write(7894,*) ekinvv, ekinv5, ekinv5

   ekinhist(hist_len) = KinEne%KinEneVV - fekinaverage
   ekinhist(hist_len-2) = ekinhist(hist_len-2) - ekinvv + ekinv5

end subroutine abf_core_lf_register_ekin_v5

!===============================================================================
! Subroutine:  abf_core_lf_register_ekin_v7
! this is leap-frog ABF version
!===============================================================================

subroutine abf_core_lf_register_ekin_v7()

    use pmf_dat
    use abf_dat
    use abf_accu
    use pmf_utils

    implicit none
    integer     :: i, m
    real(PMFDP) :: ekinvv,ekinv7,v2,v7
    ! --------------------------------------------------------------------------

    ! old VV Ekin
    ekinvv = 0.0d0
    do i=1,NumOfLAtoms
        v2 = 0.0d0
        do m=1,3
            v2 = v2 + vhist(m,i,hist_len-3)**2
        end do
        ekinvv = ekinvv + Mass(i)*v2
    end do
    ekinvv = 0.5d0*ekinvv

    ! new V7 Ekin
    ekinv7 = 0.0d0
    do i=1,NumOfLAtoms
        v2 = 0.0d0
        do m=1,3
            !            t-dt/2
            v7 = +  3.0d0*vhist(m,i,hist_len-0) - 25.0d0*vhist(m,i,hist_len-1) &
                 +150.0d0*vhist(m,i,hist_len-2) +150.0d0*vhist(m,i,hist_len-3) &
                 - 25.0d0*vhist(m,i,hist_len-4) +  3.0d0*vhist(m,i,hist_len-5)
            v7 = (1.0d0 / 256.0d0) * v7
            v2 = v2 + v7**2
        end do
        ekinv7 = ekinv7 + Mass(i)*v2
    end do
    ekinv7 = 0.5d0*ekinv7

   ! write(7894,*) ekinvv, ekinv5, ekinv5

   ekinhist(hist_len) = KinEne%KinEneVV - fekinaverage
   ekinhist(hist_len-3) = ekinhist(hist_len-3) - ekinvv + ekinv7

end subroutine abf_core_lf_register_ekin_v7

!===============================================================================

end module abf_core_lf
