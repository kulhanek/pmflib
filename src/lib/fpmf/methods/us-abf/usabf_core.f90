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
            call usabf_core_force_7p
        case(3)
            ! 10-points
            call usabf_core_force_10p
        case default
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Not implemented fmode in usabf_core_main!')
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
    use pmf_dat
    use usabf_cvs
    use pmf_cvs

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
    cvhist0(:) = cvhist1(:)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist1(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    if( fenthalpy ) then
        epothist1 = PotEne + PMFEne - fepotaverage
    else
        epothist1 = 0.0d0
    end if

    ! shift etot ene
    etothist0 = etothist1 + KinEne - fekinaverage
    if( fentropy ) then
        etothist1 = PotEne + PMFEne - fepotaverage
    else
        etothist1 = 0.0d0
    end if

    ! calculate Z matrix and its inverse
    call usabf_core_calc_Zmat

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
        icf2(:) = pxi0(:) - pxi1(:)

        ! add data to accumulator
        call usabf_accu_add_data_online(cvhist0,icf2,epothist0,etothist0)
    end if

    ! backup to the next step
    zd0  = zd1
    pxim = pxip
    v0   = Vel

    ! get us force to be applied --------------------
    call usabf_core_get_us_bias(cvhist1(:),la,TotalUSABFEnergy)

    ! project us force along coordinate
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! keep US forces to subtract them in the next step
    pxi1 = la

    return

end subroutine usabf_core_force_2p

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
    real(PMFDP) :: invh,etot2,bene
    ! --------------------------------------------------------------------------

    invh = 1.0d0 / (12.0d0 * fdtx)

    ! shift accuvalue history
    cvhist2(:) = cvhist3(:)
    cvhist3(:) = cvhist4(:)
    cvhist4(:) = cvhist5(:)
    cvhist5(:) = cvhist6(:)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist6(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist2 = epothist3
    epothist3 = epothist4
    epothist4 = epothist5
    epothist5 = epothist6
    if( fenthalpy ) then
        epothist6 = PotEne + PMFEne - fepotaverage
    else
        epothist6 = 0.0d0
    end if

    etothist2 = etothist3
    etothist3 = etothist4
    etothist4 = etothist5
    etothist5 = etothist6 + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist6 = PotEne + PMFEne - fepotaverage
    else
        etothist6 = 0.0d0
    end if

    write(1489,*) PotEne, KinEne, etothist5

    if( fstep .ge. 5 ) then
        ! calculate CV momenta
        pcvhist0(:) = pcvhist1(:)
        pcvhist1(:) = pcvhist2(:)
        pcvhist2(:) = pcvhist3(:)
        pcvhist3(:) = pcvhist4(:)

        call usabf_core_calc_Zmat

        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
                pcvhist4(i) = fzinv(i,j) * (cvhist2(j) - 8.0d0*cvhist3(j) + 8.0d0*cvhist5(j) - cvhist6(j)) * invh
            end do
        end do

        ! write(78948,*) pcvhist4
    end if

    if( fstep .ge. 7 ) then
        ! calculated ICF
        do i=1,NumOfUSABFCVs
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            icf2(i) = (pcvhist0(i) - 8.0d0*pcvhist1(i) + 8.0d0*pcvhist3(i) - pcvhist4(i)) * invh
        end do

        ! substract biasing force
        call usabf_core_get_us_bias(cvhist2(:),la,bene)
        pxi0 = icf2 - la

        ! smooth etot
        if( fsmoothetot ) then
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            etot2 = ( - 3.0d0*etothist0 + 12.0d0*etothist1 + 17.0d0*etothist2  &
                     + 12.0d0*etothist3 -  3.0d0*etothist4)/35.0d0
        else
            etot2 = etothist2
        end if

        ! record the data
        call usabf_accu_add_data_online(cvhist2,pxi0,epothist2,etot2)
    end if

    ! get US force to be applied --------------------
    call usabf_core_get_us_bias(cvhist6(:),la,TotalUSABFEnergy)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

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
    real(PMFDP) :: invh,etot3,bene
    ! --------------------------------------------------------------------------

    invh = 1.0d0 / (252.0d0 * fdtx)

    ! shift accuvalue history
    cvhist3(:) = cvhist4(:)
    cvhist4(:) = cvhist5(:)
    cvhist5(:) = cvhist6(:)
    cvhist6(:) = cvhist7(:)
    cvhist7(:) = cvhist8(:)
    cvhist8(:) = cvhist9(:)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvhist9(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist2 = epothist3
    epothist3 = epothist4
    epothist4 = epothist5
    epothist5 = epothist6
    epothist6 = epothist7
    epothist7 = epothist8
    epothist8 = epothist9
    if( fenthalpy ) then
        epothist9 = PotEne + PMFEne - fepotaverage
    else
        epothist9 = 0.0d0
    end if

    etothist0 = etothist1
    etothist1 = etothist2
    etothist2 = etothist3
    etothist3 = etothist4
    etothist4 = etothist5
    etothist5 = etothist6
    etothist6 = etothist7
    etothist7 = etothist8
    etothist8 = etothist9 + KinEne - fekinaverage  ! kinetic energy is delayed by dt
    if( fentropy ) then
        etothist9 = PotEne + PMFEne - fepotaverage
    else
        etothist9 = 0.0d0
    end if

    if( fstep .ge. 7 ) then
        ! calculate CV momenta
        pcvhist0(:) = pcvhist1(:)
        pcvhist1(:) = pcvhist2(:)
        pcvhist2(:) = pcvhist3(:)
        pcvhist3(:) = pcvhist4(:)
        pcvhist4(:) = pcvhist5(:)
        pcvhist5(:) = pcvhist6(:)

        call usabf_core_calc_Zmat

        do i=1,NumOfUSABFCVs
            do j=1,NumOfUSABFCVs
                ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
                pcvhist6(i) = fzinv(i,j) * ( 22.0d0*cvhist3(j) - 67.0d0*cvhist4(j) - 58.0d0*cvhist5(j) &
                                           + 58.0d0*cvhist7(j) + 67.0d0*cvhist8(j) - 22.0d0*cvhist9(j)) * invh
            end do
        end do

        ! write(78948,*) pcvhist4
    end if

    if( fstep .ge. 10 ) then
        ! calculated ICF
        do i=1,NumOfUSABFCVs
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            icf3(i) = ( 22.0d0*pcvhist0(i) - 67.0d0*pcvhist1(i) - 58.0d0*pcvhist2(i) &
                      + 58.0d0*pcvhist4(i) + 67.0d0*pcvhist5(i) - 22.0d0*pcvhist6(i)) * invh
        end do

        ! substract biasing force
        call usabf_core_get_us_bias(cvhist3(:),la,bene)
        pxi0 = icf3 - la

        ! smooth etot
        if( fsmoothetot ) then
            ! https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
            etot3 = ( 5.0d0*etothist0 - 30.0d0*etothist1 + 75.0d0*etothist2 + 131.0d0*etothist3 &
                   + 75.0d0*etothist4 - 30.0d0*etothist5 +  5.0d0*etothist6)/231.0d0
        else
            etot3 = etothist3
        end if

        ! record the data
        call usabf_accu_add_data_online(cvhist3,pxi0,epothist3,etot3)
    end if

    ! get US force to be applied --------------------
    call usabf_core_get_us_bias(cvhist9(:),la,TotalUSABFEnergy)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

end subroutine usabf_core_force_10p

!===============================================================================
! subroutine:  usabf_core_calc_Zmat
!===============================================================================

subroutine usabf_core_calc_Zmat()

    use pmf_utils
    use usabf_dat

    implicit none
    integer         :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfUSABFCVs
            cj = USABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
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
