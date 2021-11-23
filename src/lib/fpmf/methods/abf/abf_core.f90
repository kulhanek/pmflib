!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

    select case(fmode)
        case(1)
            ! simplified ABF algorithm
            call abf_core_force_2p
        case(2)
            ! simplified ABF algorithm
            call abf_core_force_3pA
        case(3)
            ! simplified ABF algorithm
            call abf_core_force_3pB
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    call abf_client_exchange_data(.false.)

end subroutine abf_core_main

!===============================================================================
! Subroutine:  abf_core_force_2p
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_2p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,v1,v2,f,etot
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ersthist(i)     = ersthist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)

! shift epot ene

    epothist(hist_len) = PotEne - fepotaverage
    ersthist(hist_len) = PMFEne

! shift etot ene
    select case(fekinsrc)
        case(0)
            ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
        case(1)
            ekinhist(hist_len-1)    = KinEneVV - fekinaverage  ! shifted by -dt
        case(2)
            ekinhist(hist_len)      = KinEneH - fekinaverage   ! shifted by -dt/2
    end select

    ! write(6587,*) fstep,KinEne, KinEneVV, KinEneH

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
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
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

    if( fstep .ge. hist_len ) then

        do i=1,NumOfABFCVs
            f  = 0.0d0
            v1 = 0.0d0
            v2 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f = f + zdhist(m,j,i,hist_len-1)*(vhist(m,j,hist_len-0)-vhist(m,j,hist_len-1))
                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                end do
            end do
            pxi0(i) = (f + 0.5d0*(v1+v2)) / fdtx
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)  ! unbiased estimate

        ! add data to accumulator
        select case(fekinsrc)
            case(0,1)
                etot = epothist(hist_len-1) + ersthist(hist_len-1) + ekinhist(hist_len-1)
            case(2)
                etot = epothist(hist_len-1) + ersthist(hist_len-1) + 0.5d0*(ekinhist(hist_len) + ekinhist(hist_len-1))
        end select

!        if( cvhist(1,hist_len-1) > 15 ) then
        ! write(784,*) fstep-1, cvhist(:,hist_len-1), pxi0, f/ fdtx, v1/ fdtx, v2/ fdtx, etot
!        end if

        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epothist(hist_len-1),ersthist(hist_len-1),etot)

        ! call abf_accu_add_data_record(cvhist(:,1),fzinv0,pxi0,pxi1,epothist(1),ersthist(1),ekinhist(1))
    end if

    return

end subroutine abf_core_force_2p

!===============================================================================
! Subroutine:  abf_core_force_3pA
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_3pA()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,v1,v2,f,etot
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ersthist(i)     = ersthist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        xhist(:,:,i)    = xhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)
    xhist(:,:,hist_len)    = Crd(:,:)

! shift epot ene

    epothist(hist_len) = PotEne - fepotaverage
    ersthist(hist_len) = PMFEne

! shift etot ene
    select case(fekinsrc)
        case(0)
            ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
        case(1)
            ekinhist(hist_len-1)    = KinEneVV - fekinaverage  ! shifted by -dt
        case(2)
            ekinhist(hist_len)      = KinEneH - fekinaverage   ! shifted by -dt/2
    end select

    ! write(6587,*) fstep,KinEne, KinEneVV, KinEneH

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
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
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

    if( fstep .ge. hist_len ) then

        do i=1,NumOfABFCVs
            f  = 0.0d0
            v1 = 0.0d0
            v2 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f = f + zdhist(m,j,i,hist_len-1) &
                      * (xhist(m,j,hist_len-0) - 2.0d0 * xhist(m,j,hist_len-1) + xhist(m,j,hist_len-2)) / fdtx
                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                end do
            end do
            pxi0(i) = (f + 0.5d0*(v1+v2)) / fdtx
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)  ! unbiased estimate

        ! add data to accumulator
        select case(fekinsrc)
            case(0,1)
                etot = epothist(hist_len-1) + ersthist(hist_len-1) + ekinhist(hist_len-1)
            case(2)
                etot = epothist(hist_len-1) + ersthist(hist_len-1) + 0.5d0*(ekinhist(hist_len) + ekinhist(hist_len-1))
        end select

        ! write(784,*) fstep-1, cvhist(:,hist_len-1), pxi0, f/ fdtx, v1/ fdtx, v2/ fdtx, etot

        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epothist(hist_len-1),ersthist(hist_len-1),etot)

        ! call abf_accu_add_data_record(cvhist(:,1),fzinv0,pxi0,pxi1,epothist(1),ersthist(1),ekinhist(1))
    end if

    return

end subroutine abf_core_force_3pA

!===============================================================================
! Subroutine:  abf_core_force_3p
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_3pB()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,f,etot
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ersthist(i)     = ersthist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        xhist(:,:,i)    = xhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)
    xhist(:,:,hist_len)    = Crd(:,:)

! shift epot ene

    epothist(hist_len) = PotEne - fepotaverage
    ersthist(hist_len) = PMFEne

! shift etot ene
    select case(fekinsrc)
        case(0)
            ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
        case(1)
            ekinhist(hist_len-1)    = KinEneVV - fekinaverage  ! shifted by -dt
        case(2)
            ekinhist(hist_len)      = KinEneH - fekinaverage   ! shifted by -dt/2
    end select

    ! write(6587,*) fstep,KinEne, KinEneVV, KinEneH

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
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
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

    if( fstep .ge. hist_len ) then

        do i=1,NumOfABFCVs
            f  = 0.0d0
            v  = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f = f + zdhist(m,j,i,hist_len-1) &
                      * (xhist(m,j,hist_len) - 2.0d0 * xhist(m,j,hist_len-1) + xhist(m,j,hist_len-2)) / fdtx
                    ! velocity part
                    v = v + 0.25d0*(zdhist(m,j,i,hist_len) - zdhist(m,j,i,hist_len-2)) &
                      * (vhist(m,j,hist_len) + vhist(m,j,hist_len-1))
                end do
            end do
            pxi0(i) = (f + v) / fdtx
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)  ! unbiased estimate

        ! add data to accumulator
        select case(fekinsrc)
            case(0,1)
                etot = epothist(hist_len-1) + ersthist(hist_len-1) + ekinhist(hist_len-1)
            case(2)
                etot = epothist(hist_len-1) + ersthist(hist_len-1) + 0.5d0*(ekinhist(hist_len) + ekinhist(hist_len-1))
        end select

        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epothist(hist_len-1),ersthist(hist_len-1),etot)

        ! call abf_accu_add_data_record(cvhist(:,1),fzinv0,pxi0,pxi1,epothist(1),ersthist(1),ekinhist(1))
    end if

    return

end subroutine abf_core_force_3pB

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
