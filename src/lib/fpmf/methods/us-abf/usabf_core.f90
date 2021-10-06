!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module usabf_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  usabf_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine usabf_core_main

    use usabf_trajectory
    use usabf_restart
    use usabf_output
    use usabf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    select case(fmode)
        case(1)
            ! simplified ABF algorithm - 2-points
            call usabf_core_force_2p
        case(2)
            ! 7-points
            call usabf_core_force_5p_kin
        case(3)
            ! 10-points
            call usabf_core_force_10p
        case(4)
            ! via GPR
            call usabf_core_force_gpr
        case(5)
            ! via GPR
            call usabf_core_force_gpr2
        case default
            call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Not implemented fmode in usabf_core_main!')
    end select

    call usabf_output_write
    call usabf_trajectory_write_snapshot
    call usabf_restart_update

end subroutine usabf_core_main

!===============================================================================
! Subroutine:  usabf_core_get_us_bias
! get bias from restraints
!===============================================================================

subroutine usabf_core_get_us_bias(values,gfx,bene)

    use usabf_dat

    implicit none
    real(PMFDP)     :: values(:)
    real(PMFDP)     :: gfx(:)
    real(PMFDP)     :: bene
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    bene = 0.0

    do i=1,NumOfUSABFCVs

        USABFCVList(i)%deviation = USABFCVList(i)%cv%get_deviation(values(i),USABFCVList(i)%target_value)

        USABFCVList(i)%energy = 0.5d0*USABFCVList(i)%force_constant*USABFCVList(i)%deviation**2
        bene = bene + USABFCVList(i)%energy

        gfx(i) = - USABFCVList(i)%force_constant*USABFCVList(i)%deviation
    end do

end subroutine usabf_core_get_us_bias

!===============================================================================
! Subroutine:  usabf_core_force_2p
! this is leap-frog ABF version
!===============================================================================

subroutine usabf_core_force_2p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,e
    ! --------------------------------------------------------------------------

    ! shift accuvalue history
    cvhist(:,1) = cvhist(:,2)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
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
    call usabf_core_calc_Zmat(CVContext)

    do i=1,NumOfUSABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfUSABFCVs
                    ki = USABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfUSABFCVs
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

        ! get sum
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))
        pxi0(:) = pxi0(:) + pxim(:)

        ! correct for applied bias
        call usabf_accu_get_mbcf(cvhist(:,1),la)
        pxi1(:) = pxi0(:) - la(:)

        write(489,*) fstep-1, pxi0

        ! add data to accumulator
        call usabf_accu_add_data_online(cvhist(:,1),pxi1,epothist(1),etothist(1))
    end if

    ! backup to the next step
    zd0  = zd1
    pxim = pxip
    v0   = Vel

    ! calculate bias force to be applied
    call usabf_core_get_us_bias(cvhist(:,2),la,TotalUSABFEnergy)
    call usabf_accu_add_mbcf(cvhist(:,2),la)

    ! project abf force along coordinate
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    return

end subroutine usabf_core_force_2p

!===============================================================================
! Subroutine:  usabf_core_force_5p_kin
! 5-points from CV momenta
!===============================================================================

subroutine usabf_core_force_5p_kin()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer     :: i,j,ci,m
    real(PMFDP) :: invh,etot3,v
    ! --------------------------------------------------------------------------

    invh = 1.0d0 / (12.0d0 * fdtx)

! shift accuvalue history
    do i=1,5
        cvhist(:,i) = cvhist(:,i+1)
    end do
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist(i,6) = CVContext%CVsValues(ci)
    end do

! shift epot ene
    do i=1,5
        epothist(i) = epothist(i+1)
    end do
    if( fenthalpy ) then
        epothist(6) = PotEne + PMFEne - fepotaverage
    else
        epothist(6) = 0.0d0
    end if

! shift etot ene
    do i=1,5
        etothist(i) = etothist(i+1)
    end do
    etothist(5) = etothist(5) + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist(6) = PotEne + PMFEne - fepotaverage
    else
        etothist(6) = 0.0d0
    end if

! get CV momenta from velocities
    if( fstep .ge. 2 ) then

        do i=1,NumOfUSABFCVs
            ci = USABFCVList(i)%cvindx
            v = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    v = v + 0.5d0*cvcontex0%CVsDrvs(m,j,ci)*(Vel(m,j) + v0(m,j)) ! in t-dt
                end do
            end do
            pxi0(i) = v
        end do

        ! shift history buffer
        do i=1,4
            pcvhist(:,i) = pcvhist(:,i+1)
        end do

        ! get Zmat
        call usabf_core_calc_Zmat(cvcontex0)

        ! calculate CV momenta
        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                pcvhist(i,5) = fzinv(i,j) * pxi0(j)     ! in t-dt
            end do
        end do
    end if

    v0 = Vel                ! backup velocities
    cvcontex0 = CVContext   ! backup context

    if( fstep .ge. 6 ) then
        ! calculated ICF
        do i=1,NumOfUSABFCVs
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            pxi0(i) = (pcvhist(i,1) - 8.0d0*pcvhist(i,2) + 8.0d0*pcvhist(i,4) - pcvhist(i,5)) * invh
        end do

        ! subtract biasing force
        call usabf_accu_get_mbcf(cvhist(:,3),la)
        pxi1 = pxi0 - la

        ! smooth etot
        if( fsmoothetot ) then
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            etot3 = ( - 3.0d0*etothist(1) + 12.0d0*etothist(2) + 17.0d0*etothist(3)  &
                     + 12.0d0*etothist(4) -  3.0d0*etothist(5))/35.0d0
        else
            etot3 = etothist(3)
        end if

        ! write(790,*) cvhist(:,3),pxi0,epothist(3),etothist(3)

        write(790,*) fstep-1-1-1, pxi0

        ! record the data
        call usabf_accu_add_data_online(cvhist(:,3),pxi1,epothist(3),etot3)
    end if

    ! calculate bias force to be applied
    call usabf_core_get_us_bias(cvhist(:,6),la,TotalUSABFEnergy)
    call usabf_accu_add_mbcf(cvhist(:,6),la)

    ! project abf force along coordinate
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

end subroutine usabf_core_force_5p_kin

!===============================================================================
! Subroutine:  usabf_core_force_7p
! 7-points from CV values only
!===============================================================================

subroutine usabf_core_force_7p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer     :: i,j,ci
    real(PMFDP) :: invh,etot3,bene
    ! --------------------------------------------------------------------------

    invh = 1.0d0 / (12.0d0 * fdtx)

    ! shift accuvalue history
    cvhist(:,3) = cvhist(:,4)
    cvhist(:,4) = cvhist(:,5)
    cvhist(:,5) = cvhist(:,6)
    cvhist(:,6) = cvhist(:,7)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist(i,7) = CVContext%CVsValues(ci)
    end do

! get US force to be applied --------------------
    call usabf_core_get_us_bias(cvhist(:,7),la,TotalUSABFEnergy)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! shift epot ene
    epothist(3) = epothist(4)
    epothist(4) = epothist(5)
    epothist(5) = epothist(6)
    epothist(6) = epothist(7)
    if( fenthalpy ) then
        epothist(7) = PotEne + PMFEne - fepotaverage
    else
        epothist(7) = 0.0d0
    end if

    etothist(3) = etothist(4)
    etothist(4) = etothist(5)
    etothist(5) = etothist(6)
    etothist(6) = etothist(7) + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist(7) = PotEne + PMFEne - fepotaverage
    else
        etothist(7) = 0.0d0
    end if

    if( fstep .ge. 5 ) then
        ! calculate CV momenta
        pcvhist(:,1) = pcvhist(:,2)
        pcvhist(:,2) = pcvhist(:,3)
        pcvhist(:,3) = pcvhist(:,4)
        pcvhist(:,4) = pcvhist(:,5)

        call usabf_core_calc_Zmat(CVContext)

        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
                pcvhist(i,5) = fzinv(i,j) * (cvhist(j,3) - 8.0d0*cvhist(j,4) + 8.0d0*cvhist(j,6) - cvhist(j,7)) * invh
            end do
        end do

        ! write(78948,*) pcvhist(:,5)
    end if

    if( fstep .ge. 7 ) then
        ! calculated ICF
        do i=1,NumOfUSABFCVs
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            pxi0(i) = (pcvhist(i,1) - 8.0d0*pcvhist(i,2) + 8.0d0*pcvhist(i,4) - pcvhist(i,5)) * invh
        end do

        ! substract biasing force
        call usabf_core_get_us_bias(cvhist(:,3),la,bene)
        pxi1 = pxi0 - la

        ! smooth etot
        if( fsmoothetot ) then
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            etot3 = ( - 3.0d0*etothist(1) + 12.0d0*etothist(2) + 17.0d0*etothist(3)  &
                     + 12.0d0*etothist(4) -  3.0d0*etothist(5))/35.0d0
        else
            etot3 = etothist(3)
        end if

        ! write(790,*) cvhist(:,3),pxi0,epothist(3),etothist(3)

        ! record the data
        call usabf_accu_add_data_online(cvhist(:,3),pxi1,epothist(3),etot3)
    end if

end subroutine usabf_core_force_7p

!===============================================================================
! Subroutine:  usabf_core_force_10p
! 10-points from CV values only
!===============================================================================

subroutine usabf_core_force_10p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer     :: i,j,ci
    real(PMFDP) :: invh,etot4,bene
    ! --------------------------------------------------------------------------

    invh = 1.0d0 / (252.0d0 * fdtx)

    ! shift accuvalue history
    cvhist(:,4) = cvhist(:,5)
    cvhist(:,5) = cvhist(:,6)
    cvhist(:,6) = cvhist(:,7)
    cvhist(:,7) = cvhist(:,8)
    cvhist(:,8) = cvhist(:,9)
    cvhist(:,9) = cvhist(:,10)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist(i,10) = CVContext%CVsValues(ci)
    end do

    ! get US force to be applied --------------------
    call usabf_core_get_us_bias(cvhist(:,10),la,TotalUSABFEnergy)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! shift epot ene
    epothist(3) = epothist(4)
    epothist(4) = epothist(5)
    epothist(5) = epothist(6)
    epothist(6) = epothist(7)
    epothist(7) = epothist(8)
    epothist(8) = epothist(9)
    epothist(9) = epothist(10)
    if( fenthalpy ) then
        epothist(10) = PotEne + PMFEne - fepotaverage
    else
        epothist(10) = 0.0d0
    end if

    etothist(1) = etothist(2)
    etothist(2) = etothist(3)
    etothist(3) = etothist(4)
    etothist(4) = etothist(5)
    etothist(5) = etothist(6)
    etothist(6) = etothist(7)
    etothist(7) = etothist(8)
    etothist(8) = etothist(9)
    etothist(9) = etothist(10) + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist(10) = PotEne + PMFEne - fepotaverage
    else
        etothist(10) = 0.0d0
    end if

    if( fstep .ge. 7 ) then
        ! calculate CV momenta
        pcvhist(:,1) = pcvhist(:,2)
        pcvhist(:,2) = pcvhist(:,3)
        pcvhist(:,3) = pcvhist(:,4)
        pcvhist(:,4) = pcvhist(:,5)
        pcvhist(:,5) = pcvhist(:,6)
        pcvhist(:,6) = pcvhist(:,7)

        call usabf_core_calc_Zmat(CVContext)

        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
                pcvhist(i,7) = fzinv(i,j) * ( 22.0d0*cvhist(j,4) - 67.0d0*cvhist(j,5) - 58.0d0*cvhist(j,6) &
                                            + 58.0d0*cvhist(j,8) + 67.0d0*cvhist(j,9) - 22.0d0*cvhist(j,10)) * invh
            end do
        end do

    end if

    if( fstep .ge. 10 ) then
        ! calculated ICF
        do i=1,NumOfUSABFCVs
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            pxi0(i) = ( 22.0d0*pcvhist(i,1) - 67.0d0*pcvhist(i,2) - 58.0d0*pcvhist(i,3) &
                      + 58.0d0*pcvhist(i,5) + 67.0d0*pcvhist(i,6) - 22.0d0*pcvhist(i,7)) * invh
        end do

        ! substract biasing force
        call usabf_core_get_us_bias(cvhist(:,4),la,bene)
        pxi1 = pxi0 - la

        ! smooth etot
        if( fsmoothetot ) then
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            etot4 = ( 5.0d0*etothist(1) - 30.0d0*etothist(2) + 75.0d0*etothist(3) + 131.0d0*etothist(4) &
                   + 75.0d0*etothist(5) - 30.0d0*etothist(6) +  5.0d0*etothist(7))/231.0d0
        else
            etot4 = etothist(4)
        end if

        ! record the data
        call usabf_accu_add_data_online(cvhist(:,4),pxi1,epothist(4),etot4)
    end if

end subroutine usabf_core_force_10p

!===============================================================================
! Subroutine:  usabf_core_force_gpr
! using Gaussian Process Regression
!===============================================================================

subroutine usabf_core_force_gpr()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer     :: i,j,ci,dt_index
    real(PMFDP) :: bene,mean,etot_dt_index
    ! --------------------------------------------------------------------------

    dt_index = gpr_len/2+1

! shift history buffers
    do i=1,hist_len-1
        cvhist(:,i) = cvhist(:,i+1)
        epothist(i) = epothist(i+1)
        etothist(i) = etothist(i+1)
    end do

! update new history values - CVs
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do

! apply US bias
    ! get US force to be applied --------------------
    call usabf_core_get_us_bias(cvhist(:,hist_len),la,TotalUSABFEnergy)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
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

        do i=1,NumOfUSABFCVs
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
        call usabf_core_calc_Zmat(CVContext)

        ! calculate CV momenta
        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                pcvhist(i,gpr_len) = fzinv(i,j) * pxi0(j)
            end do
        end do

        ! write(78947,*) pcvhist(:,5)
    end if

! get first derivative of CV momenta
    if( fstep .ge. hist_len ) then

        do i=1,NumOfUSABFCVs

            ! input data
            do j=1,gpr_len
                gpr_model(j) = pcvhist(i,gpr_len-j+1)
            end do

            ! solve GPR
            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)

            if( gpr_info .ne. 0 ) then
                ! throw error
                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model for CV momenta time derivatives!')
            end if

            ! calculate mean force - derivative is shift invariant
            pxi0(i) = dot_product(gpr_model,gpr_kdf)
        end do

        ! substract biasing force
        call usabf_core_get_us_bias(cvhist(:,dt_index),la,bene)
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
                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model for CV time derivatives!')
            end if

            ! calculate CV derivative in time - derivative is shift invariant
            etot_dt_index = dot_product(gpr_model,gpr_kff) + mean
        else
            etot_dt_index = etothist(dt_index)
        end if

        ! write(789,*) cvhist(:,dt_index),pxi1,epothist(dt_index),etothist(dt_index),pxi0,etot_dt_index

        ! record the data
        call usabf_accu_add_data_online(cvhist(:,dt_index),pxi1,epothist(dt_index),etot_dt_index)

    end if

end subroutine usabf_core_force_gpr

!===============================================================================
! Subroutine:  usabf_core_force_gpr2
! using Gaussian Process Regression
!===============================================================================

subroutine usabf_core_force_gpr2()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer     :: i,j,m,ci,dt_index
    real(PMFDP) :: bene,mean,etot_dt_index,v
    ! --------------------------------------------------------------------------

    dt_index = gpr_len/2+1

! shift history buffers
    do i=1,hist_len-1
        cvhist(:,i) = cvhist(:,i+1)
        epothist(i) = epothist(i+1)
        etothist(i) = etothist(i+1)
    end do

! update new history values - CVs
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do

! apply US bias
    ! get US force to be applied --------------------
    call usabf_core_get_us_bias(cvhist(:,hist_len),la,TotalUSABFEnergy)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
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
    if( fstep .ge. 2 ) then

        do i=1,NumOfUSABFCVs
            ci = USABFCVList(i)%cvindx
            v = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    v = v + 0.5d0*cvcontex0%CVsDrvs(m,j,ci)*(Vel(m,j) + v0(m,j))
                end do
            end do
            pxi0(i) = v
        end do

        ! shift history buffer
        do i=1,gpr_len-1
            pcvhist(:,i) = pcvhist(:,i+1)
        end do

        ! get Zmat
        call usabf_core_calc_Zmat(cvcontex0)

        ! calculate CV momenta
        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                pcvhist(i,gpr_len) = fzinv(i,j) * pxi0(j)
            end do
        end do
    end if

    v0 = Vel                ! backup velocities
    cvcontex0 = CVContext   ! backup context

! get first derivative of CV momenta
    if( fstep .ge. hist_len ) then

        do i=1,NumOfUSABFCVs

            ! input data
            do j=1,gpr_len
                gpr_model(j) = pcvhist(i,gpr_len-j+1)
            end do

            ! solve GPR
            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)

            if( gpr_info .ne. 0 ) then
                ! throw error
                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model for CV momenta time derivatives!')
            end if

            ! calculate mean force - derivative is shift invariant
            pxi0(i) = dot_product(gpr_model,gpr_kdf)
        end do

        ! substract biasing force
        call usabf_core_get_us_bias(cvhist(:,dt_index),la,bene)
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
                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model for CV time derivatives!')
            end if

            ! calculate CV derivative in time - derivative is shift invariant
            etot_dt_index = dot_product(gpr_model,gpr_kff) + mean
        else
            etot_dt_index = etothist(dt_index)
        end if

        ! write(790,*) cvhist(:,dt_index),pxi1,epothist(dt_index),etothist(dt_index),pxi0,etot_dt_index

        ! record the data
        call usabf_accu_add_data_online(cvhist(:,dt_index),pxi1,epothist(dt_index),etot_dt_index)

    end if

end subroutine usabf_core_force_gpr2

!===============================================================================
! subroutine:  usabf_core_calc_Zmat
!===============================================================================

subroutine usabf_core_calc_Zmat(ctx)

    use pmf_utils
    use usabf_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfUSABFCVs
            cj = USABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
            fzinv(i,j) = fz(i,j)            ! we need this for LAPACK
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfUSABFCVs .gt. 1) then
        call dgetrf(NumOfUSABFCVs,NumOfUSABFCVs,fzinv,NumOfUSABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[TABF] LU decomposition failed in usabf_calc_Zmat!')
        end if

        call dgetri(NumOfUSABFCVs,fzinv,NumOfUSABFCVs,indx,vv,NumOfUSABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Matrix inversion failed in usabf_calc_Zmat!')
        end if
    else
        fzinv(1,1)=1.0d0/fz(1,1)
    end if

    return

end subroutine usabf_core_calc_Zmat

!===============================================================================

end module usabf_core
