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
            ! GPr ABF
            call abf_core_force_gpr
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
    use abf_output

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,e
    ! --------------------------------------------------------------------------

! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,2) = CVContext%CVsValues(ci)
    end do

! shift epot ene
    epothist(1) = epothist(2)
    if( fenthalpy ) then
        epothist(2) = PotEne + PMFEne - fepotaverage
    else
        epothist(2) = 0.0d0
    end if

! shift etot ene
    etothist(1) = etothist(2) + KinEne - fekinaverage
    if( fentropy ) then
        etothist(2) = PotEne + PMFEne - fepotaverage
    else
        etothist(2) = 0.0d0
    end if

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat

    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfABFCVs
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
        pxi0(:) = pxi0(:) - pxi1(:)
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))

        ! total ABF force
        pxi1(:) = pxi0(:) + pxim(:)

        ! write(456,*) fstep, epothist0 + ekinhist1 + ersthist0

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,1),pxi1,epothist(1),etothist(1))
    end if

    ! backup to the next step
    zd0  = zd1
    pxim = pxip
    v0   = Vel

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
    use abf_output

    implicit none
    integer     :: i,j,k,m
    integer     :: ci,ck
    real(PMFDP) :: v,avg_epot,avg_etot
    ! --------------------------------------------------------------------------

! calculate acceleration in time t for all pmf atoms
    do i=1,NumOfLAtoms
        a1(:,i) = MassInv(i)*Frc(:,i)
    end do

! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)
    cvhist(:,2) = cvhist(:,3)
    cvhist(:,3) = cvhist(:,4)

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,4) = CVContext%CVsValues(ci)
    end do

! shift epot ene
    epothist(1) = epothist(2)
    epothist(2) = epothist(3)
    epothist(3) = epothist(4)
    if( fenthalpy ) then
        epothist(4) = PotEne + PMFEne - fepotaverage
    else
        epothist(4) = 0.0d0
    end if

    etothist(1) = etothist(2)
    etothist(2) = etothist(3)
    etothist(3) = etothist(4) + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist(4) = PotEne + PMFEne - fepotaverage
    else
        etothist(4) = 0.0d0
    end if

! calculate abf force to be applied -------------
    select case(feimode)
            case(0)
                call abf_accu_get_data(cvhist(:,4),la)
            case(1)
                call abf_accu_get_data_lramp(cvhist(:,4),la)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
    end select

    ! apply force filters
    if( .not. fapply_abf ) then
        ! we do not want to apply ABF force - set the forces to zero
        la(:) = 0.0d0
    end if

    ! project abf force along coordinate ------------
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            a1(:,j) = a1(:,j) + la(i) * MassInv(j) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! rest of ABF stuff -----------------------------

    ! calculate Z matrix and its inverse
    call abf_core_calc_Zmat

    ! pxip = zd0(t-dt)*[v(t-dt/2)/2 - dt*a1(t)/12]
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) - fdtx*a1(m,j)/12.0)
            end do
        end do
        pxip(i) = v
    end do

    ! ZD0(3, NumOfLAtoms, NumOfABFCVs)   <-- 1/dt * m_ksi grad ksi(r0)
    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ck = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k) * CVContext%CVsDrvs(m,j,ck)
                end do
                zd0(m,j,i) = v / fdtx
            end do
        end do
    end do

    ! pxim = zd0(t)*[v(t-dt/2)/2 + dt*a0(t-dt)/12]
    do i=1,NumOfABFCVs
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
        do i=1,NumOfABFCVs
            avg_values(i) = ABFCVList(i)%cv%get_average_value(cvhist(i,2),cvhist(i,3))
        end do

        avg_epot = 0.5d0*(epothist(2) + epothist(3))
        avg_etot = 0.5d0*(etothist(2) + etothist(3))

        ! add data to accumulator
        call abf_accu_add_data_online(avg_values,pxi0(:),avg_epot,avg_etot)
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
! Subroutine:  abf_core_force_gpr
! using Gaussian Process Regression
!===============================================================================

subroutine abf_core_force_gpr()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use abf_output

    implicit none
    integer     :: i,j,ci,dt_index
    real(PMFDP) :: mean,etot_dt_index
    ! --------------------------------------------------------------------------

    dt_index = gpr_len/2+1

! shift history buffers
    do i=1,hist_len-1
        cvhist(:,i) = cvhist(:,i+1)
        epothist(i) = epothist(i+1)
        etothist(i) = etothist(i+1)
    end do
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do

! update new history values - energy
    if( fenthalpy ) then
        epothist(hist_len) = PotEne + PMFEne - fepotaverage
    else
        epothist(hist_len) = 0.0d0
    end if

    etothist(hist_len-1) = etothist(hist_len-1) + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist(hist_len) = PotEne + PMFEne - fepotaverage
    else
        etothist(hist_len) = 0.0d0
    end if

! get first derivative of CV values in time
    if( fstep .ge. gpr_len ) then

        do i=1,NumOfABFCVs
            ! calculate mean value
            mean = 0.0d0
            do j=1,gpr_len
                mean = mean + cvhist(i,hist_len-j+1)
            end do
            mean = mean / real(gpr_len,PMFDP)

            ! shift data
            do j=1,gpr_len
                gpr_model(j) = cvhist(i,hist_len-j+1) - mean
            end do

            ! solve GPR
            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)

            if( gpr_info .ne. 0 ) then
                ! throw error
                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model for CV time derivatives!')
            end if

            ! calculate CV derivative in time - derivative is shift invariant
            pxi0(i) = dot_product(gpr_model,gpr_kdf)
        end do

        ! shift history buffer
        do i=1,gpr_len-1
            pcvhist(:,i) = pcvhist(:,i+1)
        end do

        ! get Zmat
        call abf_core_calc_Zmat

        ! calculate CV momenta
        do i=1,NumOfABFCVs
            do j=1,NumOfABFCVs
                pcvhist(i,gpr_len) = fzinv(i,j) * pxi0(j)
            end do
        end do

        ! write(78947,*) pcvhist(:,5)
    end if

! get first derivative of CV momenta
    if( fstep .ge. hist_len ) then

        do i=1,NumOfABFCVs

            ! input data
            do j=1,gpr_len
                gpr_model(j) = pcvhist(i,gpr_len-j+1)
            end do

            ! solve GPR
            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)

            if( gpr_info .ne. 0 ) then
                ! throw error
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to solve GPR model for CV momenta time derivatives!')
            end if

            ! calculate mean force - derivative is shift invariant
            pxi0(i) = dot_product(gpr_model,gpr_kdf)
        end do

        ! subtract biasing force
        la(:) = 0.0d0
        if( fapply_abf ) then
            ! get ABF force
            select case(feimode)
                case(0)
                    call abf_accu_get_data(cvhist(:,dt_index),la)
                case(1)
                    call abf_accu_get_data_lramp(cvhist(:,dt_index),la)
                case default
                    call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
            end select
        end if
        pxi1 = pxi0 - la

        ! smooth etot
        if( fsmoothetot ) then
            ! calculate mean value
            mean = 0.0d0
            do j=1,gpr_len
                mean = mean + etothist(gpr_len-j+1)
            end do
            mean = mean / real(gpr_len,PMFDP)

            ! shift data
            do j=1,gpr_len
                gpr_model(j) = etothist(gpr_len-j+1) - mean
            end do

            ! solve GPR
            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)

            if( gpr_info .ne. 0 ) then
                ! throw error
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to solve GPR model for CV time derivatives!')
            end if

            ! calculate CV derivative in time - derivative is shift invariant
            etot_dt_index = dot_product(gpr_model,gpr_kff) + mean
        else
            etot_dt_index = etothist(dt_index)
        end if

        ! write(789,*) cvhist(:,dt_index),pxi1,epothist(dt_index),etothist(dt_index),pxi0,etot_dt_index
        ! write(704,*) fstep-hist_len+dt_index, pxi0


        ! record the data
        call abf_accu_add_data_online(cvhist(:,dt_index),pxi1,epothist(dt_index),etot_dt_index)

    end if

! calculate abf force to be applied -------------
    if( fapply_abf ) then
        ! get ABF force
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvhist(:,hist_len),la)
            case(1)
                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
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

end subroutine abf_core_force_gpr

!===============================================================================
! subroutine:  abf_core_calc_Zmat
!===============================================================================

subroutine abf_core_calc_Zmat()

    use pmf_utils
    use abf_dat

    implicit none
    integer         :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfABFCVs
            cj = ABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            fzinv(i,j) = fz(i,j)            ! we need this for LAPACK
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then
        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzinv,NumOfABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_calc_Zmat!')
        end if

        call dgetri(NumOfABFCVs,fzinv,NumOfABFCVs,indx,vv,NumOfABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_calc_Zmat!')
        end if
    else
        fzinv(1,1)=1.0d0/fz(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat

!===============================================================================

end module abf_core
