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
        case(4)
            call abf_core_vv_force_2pV
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
    vhist(:,:,hist_len-1)   = Vel(:,:) * ftds_vel_scale ! velocities are in t-dt
    fzinvhist(:,:,hist_len) = fzinv(:,:)
    call abf_core_update_cvder

    do i=1,NumOfABFCVs
        select case(abf_p2_vx)
        case(2)
            ! -1
            vint(:,:) =  vhist(:,:,hist_len-1)
            vidx =-1
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_vx in abf_core_vv_force_2pV!')
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
    !    write(14789,*) fstep, pxia
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
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented abf_p2_px in abf_core_vv_force_2pV!')
        end select
        pxiv(i) = 0.0d0
    end do

    ! subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)
    call abf_core_register_rawdata(cvhist(:,hist_len-6),pxif,pxis,pxiv,micfhist(:,hist_len-6), &
                           epothist(hist_len-6),ersthist(hist_len-6),ekinhist(hist_len-6))

end subroutine abf_core_vv_force_2pV

!===============================================================================

end module abf_core_vv
