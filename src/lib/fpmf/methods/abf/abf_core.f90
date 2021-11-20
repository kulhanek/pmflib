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
            ! original ABF algorithm
            call abf_core_force_4p
        case(3)
            ! LPF+SG ABF algorithm
            call abf_core_force_lpf_sg
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
    real(PMFDP)            :: v,e,etot,zc
    ! --------------------------------------------------------------------------

! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)

    do i=1,NumOfBiasedABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,2) = CVContext%CVsValues(ci)
    end do

! shift epot ene
    epothist(1) = epothist(2)
    epothist(2) = PotEne - fepotaverage

    ersthist(1) = ersthist(2)
    ersthist(2) = PMFEne

! shift etot ene
    ekinhist(1) = KinEne - fekinaverage ! shifted by -dt

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)
    call abf_core_calc_Zmat_all(CVContext)

    do i=1,NumOfBiasedABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfBiasedABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfBiasedABFCVs
        v = 0.0d0
        e = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! zd0 in t-dt
                ! Vel in t-1/2dt
                ! v0 (OldVel) in t-3/2dt
                ! a <- Vel(m,j)-v0(m,j))/fdtx in t-dt, do not use forces (Frc) because of SHAKE
                v = v + zd0(m,j,i)*(Vel(m,j)-v0(m,j))
                ! zd1 in t
                ! zd0 in t-dt
                ! vel in t-1/2dt
                e = e + (zd1(m,j,i)-zd0(m,j,i))* Vel(m,j)
            end do
        end do
        pxi0(i) = v / fdtx ! in t-dt
        pxip(i) = e / fdtx ! in t-1/2dt
    end do

    if( fstep .ge. 4 ) then

        ! complete ICF in t-dt
        ! pxi0 in t-dt
        ! pxi1 - old ABF forces in t-dt
        ! pxip in t-1/2dt
        ! pxim in t-3/2dt
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))

        ! total ABF force
        pxi0(:) = pxi0(:) + pxim(:)     ! biased estimate
        pxi0(:) = pxi0(:) - pxi1(:)     ! unbiased estimate

        !write(789,*) fstep-1,fzdet0,fzdetA0,fzdetA0/fzdet0

        ! add data to accumulator
        etot = epothist(1) + ersthist(1) + ekinhist(1)

        zc = sqrt(fzdetall0/fzdet0)
        ! write(4589,*) fstep-1, zc, fzdet0, fzdetall0

        call abf_accu_add_data_online(cvhist(:,1),pxi0,epothist(1),ersthist(1),etot,zc)

        ! call abf_accu_add_data_record(cvhist(:,1),fzinv0,pxi0,pxi1,epothist(1),ersthist(1),ekinhist(1))
    end if

    ! backup to the next step
    zd0         = zd1
    pxim        = pxip
    v0          = Vel
    fzinv0      = fzinv
    fzdet0      = fzdet
    fzdetall0   = fzdetall

    ! apply ABF bias
    la(:) = 0.0d0

    ! apply force filters
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvhist(:,2),la)
            case(1)
                call abf_accu_get_data_lramp(cvhist(:,2),la)
            case(2)
                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                    call abf_accu_get_data_ksmooth(cvhist(:,2),la)
                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
            case(3)
                call abf_accu_get_data_lsmooth(cvhist(:,2),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfBiasedABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if

    ! keep ABF forces to subtract them in the next step
    pxi1 = la

    return

end subroutine abf_core_force_2p

!===============================================================================
! Subroutine:  abf_core_force_4p
! this is leap-frog ABF version, original implementation
!===============================================================================

subroutine abf_core_force_4p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer     :: i,j,k,m
    integer     :: ci,ck
    real(PMFDP) :: v
    ! --------------------------------------------------------------------------

! calculate acceleration in time t for all pmf atoms
    do i=1,NumOfLAtoms
        a1(:,i) = MassInv(i)*Frc(:,i)
    end do

! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)
    cvhist(:,2) = cvhist(:,3)
    cvhist(:,3) = cvhist(:,4)

    do i=1,NumOfBiasedABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,4) = CVContext%CVsValues(ci)
    end do

! apply force filters
    la(:) = 0.0d0
    if( fapply_abf ) then
    ! calculate abf force to be applied -------------
        select case(feimode)
                case(0)
                    call abf_accu_get_data(cvhist(:,4),la)
                case(1)
                    call abf_accu_get_data_lramp(cvhist(:,4),la)
                case(2)
                    call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                        call abf_accu_get_data_ksmooth(cvhist(:,4),la)
                    call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
                case(3)
                    call abf_accu_get_data_lsmooth(cvhist(:,4),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select
    ! project abf force along coordinate ------------
        do i=1,NumOfBiasedABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                a1(:,j) = a1(:,j) + la(i) * MassInv(j) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
   end if

! rest of ABF stuff -----------------------------

    ! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

    ! pxip = zd0(t-dt)*[v(t-dt/2)/2 - dt*a1(t)/12]
    do i=1,NumOfBiasedABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) - fdtx*a1(m,j)/12.0)
            end do
        end do
        pxip(i) = v
    end do

    ! ZD0(3, NumOfLAtoms, NumOfBiasedABFCVs)   <-- 1/dt * m_ksi grad ksi(r0)
    do i=1,NumOfBiasedABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfBiasedABFCVs
                    ck = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k) * CVContext%CVsDrvs(m,j,ck)
                end do
                zd0(m,j,i) = v / fdtx
            end do
        end do
    end do

    ! pxim = zd0(t)*[v(t-dt/2)/2 + dt*a0(t-dt)/12]
    do i=1,NumOfBiasedABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) + fdtx*a0(m,j)/12.0)
            end do
        end do
        pxim(i) = v
    end do

    !a0 <-- a1
    a0(:,:) = a1(:,:)

    ! pxi0 <--- pxi0 + pxip
    pxi0(:) = pxi0(:) + pxip(:)

    ! update accumulator - we need at least four samples
    if( fstep .ge. 4 ) then
        ! calculate coordinate values at time t-3/2dt
        do i=1,NumOfBiasedABFCVs
            cvave(i) = ABFCVList(i)%cv%get_average_value(cvhist(i,2),cvhist(i,3))
        end do

        ! add data to accumulator
        ! FIXME
        ! call abf_accu_add_data_online(cvave(:),pxi0,0.0d0,0.0d0,0.0d0)
    end if

    ! pxi0 <--- -pxip + pxim + pxi1 - la/2
    pxi0(:) = -pxip(:) + pxim(:) + pxi1(:) - 0.5d0*la(:)

    ! pxi1 <--- -pxim - la/2
    pxi1(:) = -pxim(:) - 0.5d0*la(:)

    ! project updated acceleration back
    do i=1,NumOfLAtoms
        Frc(:,i) = a1(:,i)*Mass(i)
    end do

    return

end subroutine abf_core_force_4p

!===============================================================================
! Subroutine:  abf_core_force_lpf_sg
! this is leap-frog ABF version
! low-pass filter + SG filter for differentiation
!===============================================================================

subroutine abf_core_force_lpf_sg()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i, j, ci
    real(PMFDP)            :: v, epot, erst, ekin, etot
    ! --------------------------------------------------------------------------

!    if( mod(fstep,10000) .eq. 0 ) then
!        rewind(789)
!        do c=0.87,3.05,0.001
!            cvcur(1) = c
!            call abf_accu_get_data(cvcur,pxi0)
!            call abf_accu_get_data_lsmooth(cvcur,pxi1)
!            write(789,*) c, pxi0(1), pxi1(1)
!        end do
!    end if

! get current value of CVs
    do i=1,NumOfBiasedABFCVs
        ci = ABFCVList(i)%cvindx
        cvcur(i)  = CVContext%CVsValues(ci)
    end do

! apply force filters
    pxi1(:) = 0.0d0
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvcur,pxi1)
            case(1)
                call abf_accu_get_data_lramp(cvcur,pxi1)
            case(2)
                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                    call abf_accu_get_data_ksmooth(cvcur,pxi1)
                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
            case(3)
                call abf_accu_get_data_lsmooth(cvcur,pxi1)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode in abf_core_force_lpf_sg!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfBiasedABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + pxi1(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if

! get zvinv
    call abf_core_calc_Zmat(CVContext)

! update history buffers
    do i=1,hist_len-1
        xihist(i,:)     = xihist(i+1,:)
        micfhist(i,:)   = micfhist(i+1,:)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        zinvhist(i,:,:) = zinvhist(i+1,:,:)
    end do

    xihist(hist_len,:)      = cvcur(:)
    micfhist(hist_len,:)    = pxi1(:)
    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    ekinhist(hist_len)      = KinEne - fekinaverage   ! in t-dt
    zinvhist(hist_len,:,:)  = fzinv(:,:)

    if( fstep .lt. hist_len ) return

    if( fsgsmoothall ) then
        do i=1,NumOfBiasedABFCVs
            cvcur(i) = dot_product(sg_c0(:),xihist(1:hist_len-1,i))
            cv1dr(i) = dot_product(sg_c1(:),xihist(1:hist_len-1,i))
            cv2dr(i) = dot_product(sg_c2(:),xihist(1:hist_len-1,i))
            pxi1(i)  = dot_product(sg_c0(:),micfhist(1:hist_len-1,i))
            do j=1,NumOfBiasedABFCVs
                fzinv(i,j)     = dot_product(sg_c0(:),zinvhist(1:hist_len-1,i,j))
                fzinv0(i,j)    = dot_product(sg_c1(:),zinvhist(1:hist_len-1,i,j))
            end do
        end do

        epot = dot_product(sg_c0(:),epothist(1:hist_len-1))
        erst = dot_product(sg_c0(:),ersthist(1:hist_len-1))
        ekin = dot_product(sg_c0(:),ekinhist(2:hist_len))
        etot = epot + erst + ekin

        do i=1,NumOfBiasedABFCVs
            v = 0.0d0
            do j=1,NumOfBiasedABFCVs
                v = v +  fzinv0(j,i)*cv1dr(j) + fzinv(j,i)*cv2dr(j)
            end do
            pxi0(i) = v
        end do

        pxi0(:) = pxi0(:) - pxi1(:)

        ! add data to accumulator
        ! FIXME call abf_accu_add_data_online(cvcur,pxi0,epot,erst,etot)
    else
        do i=1,NumOfBiasedABFCVs
            cv1dr(i) = dot_product(sg_c1(:),xihist(1:hist_len-1,i))
            cv2dr(i) = dot_product(sg_c2(:),xihist(1:hist_len-1,i))
            do j=1,NumOfBiasedABFCVs
                fzinv(i,j)     = dot_product(sg_c0(:),zinvhist(1:hist_len-1,i,j))
                fzinv0(i,j)    = dot_product(sg_c1(:),zinvhist(1:hist_len-1,i,j))
            end do
        end do

        epot = epothist(hist_len/2)     ! hist_len is +1 bigger than fsgframelen
        erst = ersthist(hist_len/2)
        ekin = ekinhist(hist_len/2+1)   ! delayed
        etot = epot + erst + ekin

        do i=1,NumOfBiasedABFCVs
            v = 0.0d0
            do j=1,NumOfBiasedABFCVs
                v = v +  fzinv0(j,i)*cv1dr(j) + fzinv(j,i)*cv2dr(j)
            end do
            pxi0(i) = v
        end do

        pxi0(:)     = pxi0(:) - micfhist(hist_len/2,:)
        cvcur(:)    = xihist(hist_len/2,:)

        ! add data to accumulator
        ! FIXME call abf_accu_add_data_online(cvcur,pxi0,epot,erst,etot)
    end if

!    call abf_accu_add_data_record(cvhist(:,hist_len-1),fzinv0,pxi0,pxi1, &
!         epothist(hist_len-1),ersthist(1),ekinhist(1))

end subroutine abf_core_force_lpf_sg

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
    do i=1,NumOfBiasedABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfBiasedABFCVs
            cj = ABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfBiasedABFCVs .gt. 1) then

        fzinv(:,:)  = fz(:,:)

        call dgetrf(NumOfBiasedABFCVs,NumOfBiasedABFCVs,fzinv,NumOfBiasedABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if

        fzdet = 1.0d0
        ! and finally determinant
        do i=1,NumOfBiasedABFCVs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * fzinv(i,i)
            else
                fzdet = fzdet * fzinv(i,i)
            end if
        end do

        call dgetri(NumOfBiasedABFCVs,fzinv,NumOfBiasedABFCVs,indx,vv,NumOfBiasedABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_core_calc_Zmat!')
        end if
    else
        fzdet       = fz(1,1)
        fzinv(1,1)  = 1.0d0/fz(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat

!===============================================================================
! subroutine:  abf_core_calc_Zmat_all
!===============================================================================

subroutine abf_core_calc_Zmat_all(ctx)

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
            fzall(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fzall(i,j) = fzall(i,j) + MassInv(k)*dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
        end do
    end do

    ! and get determinant - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then

        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzall,NumOfABFCVs,indxall,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat_all!')
        end if

        fzdetall = 1.0d0
        ! and finally determinant
        do i=1,NumOfABFCVs
            if( indxall(i) .ne. i ) then
                fzdetall = - fzdetall * fzall(i,i)
            else
                fzdetall = fzdetall * fzall(i,i)
            end if
        end do
    else
        fzdetall = fzall(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat_all

!===============================================================================

end module abf_core
