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
            ! standard ABF algorithm
            call usabf_core_force_4p
        case(2)
            ! numerical differentiation
            call usabf_core_force_2p
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

subroutine usabf_core_get_us_bias(values,gfx)

    use usabf_dat
    use pmf_dat
    use usabf_cvs
    use pmf_cvs

    implicit none
    real(PMFDP)     :: values(:)
    real(PMFDP)     :: gfx(:)
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    TotalUSABFEnergy = 0.0

    do i=1,NumOfUSABFCVs

        USABFCVList(i)%deviation = USABFCVList(i)%cv%get_deviation(values(i),USABFCVList(i)%target_value)

        USABFCVList(i)%energy = 0.5d0*USABFCVList(i)%force_constant*USABFCVList(i)%deviation**2
        TotalUSABFEnergy = TotalUSABFEnergy + USABFCVList(i)%energy

        gfx(i) = - USABFCVList(i)%force_constant*USABFCVList(i)%deviation
    end do

end subroutine usabf_core_get_us_bias

!===============================================================================
! Subroutine:  usabf_core_force_4p
! this is leap-frog ABF version, original implementation
!===============================================================================

subroutine usabf_core_force_4p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use usabf_dat
    use usabf_accu
    use usabf_output

    implicit none
    integer     :: i,j,k,m
    integer     :: ci,ck
    real(PMFDP) :: v,avg_epot,avg_ekin,avg_erst
    ! --------------------------------------------------------------------------

    ! calculate acceleration in time t for all pmf atoms
    do i=1,NumOfLAtoms
        a1(:,i) = MassInv(i)*Frc(:,i)
    end do

    ! shift accuvalue history
    cvaluehist0(:) = cvaluehist1(:)
    cvaluehist1(:) = cvaluehist2(:)
    cvaluehist2(:) = cvaluehist3(:)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvaluehist3(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    epothist1 = epothist2
    epothist2 = epothist3
    if( fenthalpy .or. fentropy ) then
        epothist3 = PotEne + fepotoffset
    else
        epothist3 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    ekinhist1 = ekinhist2
    ekinhist2 = ekinhist3
    if( fentropy ) then
        ekinhist3 = KinEne + fekinoffset
    else
        ekinhist3 = 0.0d0
    end if

    ! shift erst ene
    ersthist0 = ersthist1
    ersthist1 = ersthist2
    ersthist2 = ersthist3
    if( fentropy ) then
        ersthist3 = PMFEne
    else
        ersthist3 = 0.0d0
    end if

    ! get us force to be applied --------------------
    call usabf_core_get_us_bias(cvaluehist3(:),la)

    ! project abf force along coordinate ------------
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            a1(:,j) = a1(:,j) + la(i) * MassInv(j) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! rest of ABF stuff -----------------------------

    ! calculate Z matrix and its inverse
    call usabf_core_calc_Zmat

    ! pxip = zd0(t-dt)*[v(t-dt/2)/2 - dt*a1(t)/12]
    do i=1,NumOfUSABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) - fdtx*a1(m,j)/12.0)
            end do
        end do
        pxip(i) = v
    end do

    ! ZD0(3, NumOfLAtoms, NumOfUSABFCVs)   <-- 1/dt * m_ksi grad ksi(r0)
    do i=1,NumOfUSABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfUSABFCVs
                    ck = USABFCVList(k)%cvindx
                    v = v + fzinv(i,k) * CVContext%CVsDrvs(m,j,ck)
                end do
                zd0(m,j,i) = v / fdtx
            end do
        end do
    end do

    ! pxim = zd0(t)*[v(t-dt/2)/2 + dt*a0(t-dt)/12]
    do i=1,NumOfUSABFCVs
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
        do i=1,NumOfUSABFCVs
            avg_values(i) = USABFCVList(i)%cv%get_average_value(cvaluehist1(i),cvaluehist2(i))
        end do

        avg_epot = 0.5d0*(epothist1 + epothist2) ! t - 3/2*dt
        avg_ekin = 0.5d0*(ekinhist2 + ekinhist3) ! t - 1/2*dt; ekin already shifted by -dt
        avg_erst = 0.5d0*(ersthist1 + ersthist2) ! t - 3/2*dt

        ! add data to accumulator
        call usabf_accu_add_data_online(cvaluehist0,pxi0(:),avg_epot,avg_ekin,avg_erst)
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

end subroutine usabf_core_force_4p

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
    cvaluehist0(:) = cvaluehist1(:)

    ! save coordinate value to history
    do i=1,NumOfUSABFCVs
        ci = USABFCVList(i)%cvindx
        cvaluehist1(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    if( fenthalpy .or. fentropy ) then
        epothist1 = PotEne + fepotoffset
    else
        epothist1 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    if( fentropy ) then
        ekinhist1 = KinEne + fekinoffset
    else
        ekinhist1 = 0.0d0
    end if

    ! shift erst ene
    ersthist0 = ersthist1
    if( fentropy ) then
        ersthist1 = PMFEne
    else
        ersthist1 = 0.0d0
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
        pxi0(:) = pxi0(:) - pxi1(:)
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))

        ! get sum
        pxi0(:) = pxi0(:) + pxim(:)

        ! add data to accumulator
        call usabf_accu_add_data_online(cvaluehist0,pxi0(:),epothist0,ekinhist1,ersthist0)
    end if

    ! backup to the next step
    zd0  = zd1
    pxim = pxip
    v0   = Vel

    ! get us force to be applied --------------------
    call usabf_core_get_us_bias(cvaluehist1(:),la)

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
