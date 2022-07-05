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

module abf_core_vv

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_core_vv_main
! this is velocity Verlet ABF version
!===============================================================================

subroutine abf_core_vv_main

    use abf_trajectory
    use abf_restart
    use abf_output
    use abf_client
    use abf_dat
    use pmf_utils
    ! --------------------------------------------------------------------------

    ! if we have enough data - run ABF
    select case(fmode)
! not supported any more: 2, 3
        case(1)
            call abf_core_vv_force_3pV1
        case(3)
            call abf_core_vv_force_2pF
        case(4)
            call abf_core_vv_force_2pV
        case(5)
            call abf_core_vv_force_2pX
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_vv_main!')
    end select

    ! the remaining
    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    if( fupdate_abf ) then
        call abf_client_exchange_data(.false.)
    end if

end subroutine abf_core_vv_main

!===============================================================================
! Subroutine:  abf_core_vv_shake
! forces from SHAKE
!===============================================================================

subroutine abf_core_vv_shake

    use abf_dat
    use pmf_utils
    ! --------------------------------------------------------------------------

    select case(fmode)
        case(1,2,4,5)
            ! ignored
        case(3)
            ! get forces from SHAKE
            shist(:,:,hist_len) = SHAKEFrc(:,:)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_vv_shake!')
    end select

  !  write(74896,*) fstep, SHAKEFrc

end subroutine abf_core_vv_shake

!===============================================================================
! Subroutine:  abf_core_vv_rattle
! forces from RATTLE
!===============================================================================

subroutine abf_core_vv_rattle

    use abf_dat
    use pmf_utils
    ! --------------------------------------------------------------------------

    select case(fmode)
        case(1,2,4,5)
            ! ignored
        case(3)
            ! get forces from RATTLE
          !  shist(:,:,hist_len-1) = shist(:,:,hist_len-1) + 0.5d0*RATTLEFrc(:,:)
            rhist(:,:,hist_len) = RATTLEFrc(:,:)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_vv_rattle!')
    end select

   ! write(74895,*) fstep, RATTLEFrc

end subroutine abf_core_vv_rattle

!===============================================================================
! Subroutine:  abf_core_vv_force_3pV1
! this is velocity verlet ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_vv_force_3pV1()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

    implicit none
    integer                :: i,j,m,fidx
    real(PMFDP)            :: f1,v1
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
    end do

    vhist(:,:,hist_len-1) = Vel(:,:) * ftds_vel_scale  ! shifted by -dt

    call abf_core_update_zdhist

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        v1 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                select case(abf_p2_px)
                    case(3)
                        f1 = f1 + zdhist(m,j,i,hist_len-2) &
                                * (1.0d0/2.0d0) * (vhist(m,j,hist_len-1) - vhist(m,j,hist_len-3))
                        v1 = v1 +  (1.0d0/2.0d0)*(zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-3)) &
                                                * vhist(m,j,hist_len-2)
                        fidx = -2
                    case(5)
                        f1 = f1 + zdhist(m,j,i,hist_len-3) &
                                * (1.0d0/2.0d0) * (vhist(m,j,hist_len-2) - vhist(m,j,hist_len-4))
                        v1 = v1 + (1.0d0/12.0d0)*(-1.0d0*zdhist(m,j,i,hist_len-1) +8.0d0*zdhist(m,j,i,hist_len-2) &
                                                  -8.0d0*zdhist(m,j,i,hist_len-4) +1.0d0*zdhist(m,j,i,hist_len-5)) &
                                                * vhist(m,j,hist_len-3)
                        fidx = -3
                    case(7)
                        f1 = f1 + zdhist(m,j,i,hist_len-4) &
                                * (1.0d0/2.0d0) * (vhist(m,j,hist_len-3) - vhist(m,j,hist_len-5))
                        v1 = v1 + (1.0d0/60.0d0)*( +1.0d0*zdhist(m,j,i,hist_len-1)  -9.0d0*zdhist(m,j,i,hist_len-2) &
                                                  +45.0d0*zdhist(m,j,i,hist_len-3) -45.0d0*zdhist(m,j,i,hist_len-5) &
                                                   +9.0d0*zdhist(m,j,i,hist_len-6)  -1.0d0*zdhist(m,j,i,hist_len-7)) &
                                                * vhist(m,j,hist_len-4)
                        fidx = -4
                    case default
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_vv_force_2pF!')
                end select
            end do
        end do
        pxif(i) = f1*ifdtx
        pxis(i) = 0.0d0
        pxiv(i) = v1*ifdtx

        write(78948+i,*) fstep, pxif(i), pxis(i), pxiv(i)
    end do

   ! write(78948,*) fstep, pxif(1), pxis(1)

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len+fidx),pxif,pxis,pxiv,micfhist(:,hist_len+fidx), &
                           epothist(hist_len+fidx),ersthist(hist_len+fidx),ekinhist(hist_len+fidx))

end subroutine abf_core_vv_force_3pV1

!===============================================================================
! Subroutine:  abf_core_vv_force_2pF
! this is velocity verlet ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_vv_force_2pF()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

    implicit none
    integer                :: i,j,m,fidx
    real(PMFDP)            :: f1,s1,r1,v1
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        fhist(:,:,i)    = fhist(:,:,i+1)
        shist(:,:,i)    = shist(:,:,i+1)
        rhist(:,:,i)    = rhist(:,:,i+1)
    end do

    vhist(:,:,hist_len-1) = Vel(:,:) * ftds_vel_scale  ! shifted by -dt
    fhist(:,:,hist_len)   = Frc(:,:)
    ! shist is updated latter

    call abf_core_update_zdhist

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs
        f1 = 0.0d0
        s1 = 0.0d0
        r1 = 0.0d0
        v1 = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                select case(abf_p2_px)
                    case(3)
                        f1 = f1 + zdhist(m,j,i,hist_len-2) * fhist(m,j,hist_len-2) * MassInv(j)
                        s1 = s1 + zdhist(m,j,i,hist_len-2) * shist(m,j,hist_len-2) * MassInv(j)
                        r1 = r1 + zdhist(m,j,i,hist_len-2) * rhist(m,j,hist_len-2) * MassInv(j)
                        v1 = v1 +  (1.0d0/2.0d0)*(zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-3)) &
                                                * vhist(m,j,hist_len-2)
                        fidx = -2
                    case(5)
                        f1 = f1 + zdhist(m,j,i,hist_len-3) * fhist(m,j,hist_len-3) * MassInv(j)
                        s1 = s1 + zdhist(m,j,i,hist_len-3) * shist(m,j,hist_len-3) * MassInv(j)
                        r1 = r1 + zdhist(m,j,i,hist_len-3) * rhist(m,j,hist_len-3) * MassInv(j)
                        v1 = v1 + (1.0d0/12.0d0)*(-1.0d0*zdhist(m,j,i,hist_len-1) +8.0d0*zdhist(m,j,i,hist_len-2) &
                                                  -8.0d0*zdhist(m,j,i,hist_len-4) +1.0d0*zdhist(m,j,i,hist_len-5)) &
                                                * vhist(m,j,hist_len-3)
                        fidx = -3
                    case(7)
                        f1 = f1 + zdhist(m,j,i,hist_len-4) * fhist(m,j,hist_len-4) * MassInv(j)
                        s1 = s1 + zdhist(m,j,i,hist_len-4) * shist(m,j,hist_len-4) * MassInv(j)
                        r1 = r1 + zdhist(m,j,i,hist_len-4) * rhist(m,j,hist_len-4) * MassInv(j)
                        v1 = v1 + (1.0d0/60.0d0)*( +1.0d0*zdhist(m,j,i,hist_len-1)  -9.0d0*zdhist(m,j,i,hist_len-2) &
                                                  +45.0d0*zdhist(m,j,i,hist_len-3) -45.0d0*zdhist(m,j,i,hist_len-5) &
                                                   +9.0d0*zdhist(m,j,i,hist_len-6)  -1.0d0*zdhist(m,j,i,hist_len-7)) &
                                                * vhist(m,j,hist_len-4)
                        fidx = -4
                    case default
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_vv_force_2pF!')
                end select
            end do
        end do
        pxif(i) = f1
        pxis(i) = s1
        pxiv(i) = v1*ifdtx
        ! write(78948+i,*) fstep, f1, s1, r1, v1*ifdtx
    end do

    ! write(78948,*) fstep, pxif(1), pxis(1)

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len+fidx),pxif,pxis,pxiv,micfhist(:,hist_len+fidx), &
                           epothist(hist_len+fidx),ersthist(hist_len+fidx),ekinhist(hist_len+fidx))

end subroutine abf_core_vv_force_2pF

!===============================================================================
! Subroutine:  abf_core_vv_force_2pV
! this is velocity Verlet ABF version, numerical via velocities
!===============================================================================

subroutine abf_core_vv_force_2pV()

    use pmf_dat
    use pmf_cvs
    use abf_dat
    use pmf_utils
    use abf_core

    implicit none
    integer                :: i,j,m,fidx
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

    ! update core history and apply bias
    call abf_core_update_history

    ! shift history buffers
    do i=1,hist_len-1
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
        cvderhist(:,:,:,i)  = cvderhist(:,:,:,i+1)
        xphist(:,i)         = xphist(:,i+1)
    end do
    fzinvhist(:,:,hist_len) = fzinv(:,:)
    call abf_core_update_cvder

    ! velocities are in t-dt
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + cvderhist(m,j,i,hist_len-1)*Vel(m,j)* ftds_vel_scale
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
            v = v + fzinvhist(i,j,hist_len-1)*pxia(j)
        end do
        xphist(i,hist_len-1) = v
    end do

    if( fstep .le. 2*hist_len ) return

    do i=1,NumOfABFCVs
        select case(abf_p2_px)
    ! central differences
        case(3)
            pxif(i) = 0.5d0*(xphist(i,hist_len-1) - xphist(i,hist_len-3))*ifdtx
            fidx = -2
        case(5)
            pxif(i) = (1.0d0/12.0d0)*(      -xphist(i,hist_len-1) + 8.0d0*xphist(i,hist_len-2) &
                                      -8.0d0*xphist(i,hist_len-4)      + xphist(i,hist_len-5))*ifdtx
            fidx = -3
        case(7)
            pxif(i) = (1.0d0/60.0d0)*(        xphist(i,hist_len-1)  -9.0d0*xphist(i,hist_len-2) &
                                      +45.0d0*xphist(i,hist_len-3) -45.0d0*xphist(i,hist_len-5) &
                                       +9.0d0*xphist(i,hist_len-6)        -xphist(i,hist_len-7))*ifdtx
            fidx = -4
     ! central+backward differences
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
     ! backward differences
        case(23)
            pxif(i) = (1.0d0/2.0d0)*(+3.0d0*xphist(i,hist_len-1) -4.0d0*xphist(i,hist_len-2) &
                                     +1.0d0*xphist(i,hist_len-3))*ifdtx
            fidx = -1
        case(24)
            pxif(i) = (1.0d0/6.0d0)*(+11.0d0*xphist(i,hist_len-1) -18.0d0*xphist(i,hist_len-2) &
                                      +9.0d0*xphist(i,hist_len-3)  -2.0d0*xphist(i,hist_len-4))*ifdtx
            fidx = -1
        case(25)
            pxif(i) = (1.0d0/12.0d0)*(+25.0d0*xphist(i,hist_len-1) -48.0d0*xphist(i,hist_len-2) &
                                      +36.0d0*xphist(i,hist_len-3) -16.0d0*xphist(i,hist_len-4) &
                                       +3.0d0*xphist(i,hist_len-5))*ifdtx
            fidx = -1
        case(26)
            pxif(i) = (1.0d0/60.0d0)*(+137.0d0*xphist(i,hist_len-1) -300.0d0*xphist(i,hist_len-2) &
                                      +300.0d0*xphist(i,hist_len-3) -200.0d0*xphist(i,hist_len-4) &
                                       +75.0d0*xphist(i,hist_len-5)  -12.0d0*xphist(i,hist_len-6))*ifdtx
            fidx = -1
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_vv_force_2pV!')
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

end subroutine abf_core_vv_force_2pV

!===============================================================================
! Subroutine:  abf_core_vv_force_2pX
! this is velocity verlet ABF version, numerical via CV values
!===============================================================================

subroutine abf_core_vv_force_2pX()

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

end subroutine abf_core_vv_force_2pX

!===============================================================================

end module abf_core_vv
