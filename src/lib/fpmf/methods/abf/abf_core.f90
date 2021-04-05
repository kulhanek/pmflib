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
            ! standard algorithm
            call abf_core_force
        case(2)
            ! numerical differentiation
            call abf_core_force_numeric
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    call abf_client_exchange_data(.false.)

end subroutine abf_core_main

!===============================================================================
! Subroutine:  abf_core_force
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine abf_core_force()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accumulator
    use abf_output

    implicit none
    integer     :: i,j,k,m
    integer     :: ci,ck
    real(PMFDP) :: v,avg_etot
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
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvaluehist3(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    epothist1 = epothist2
    epothist2 = epothist3
    if( faccuepot ) then
        epothist3 = PotEne
    else
        epothist3 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    ekinhist1 = ekinhist2
    ekinhist2 = ekinhist3
    if( faccuekin ) then
        ekinhist3 = KinEne
    else
        ekinhist3 = 0.0d0
    end if

    ! calculate abf force to be applied -------------
    select case(feimode)
        case(1)
            call abf_accumulator_get_data1(cvaluehist3(:),la)
        case(2)
            call abf_accumulator_get_data2(cvaluehist3(:),la)
        case(3)
            call abf_accumulator_get_data3(cvaluehist3(:),la)
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
            avg_values(i) = ABFCVList(i)%cv%get_average_value(cvaluehist1(i),cvaluehist2(i))
        end do

        avg_etot =            0.5d0*(epothist1 + epothist2) ! t - 3/2*dt
        avg_etot = avg_etot + 0.5d0*(ekinhist2 + ekinhist3) ! t - 1/2*dt; ekin already shifted by -dt
        avg_etot = avg_etot - ftotoffset

        ! add data to accumulator
        select case(feimode)
            case(1)
                call abf_accumulator_add_data12(avg_values,avg_etot,pxi0(:))
            case(2)
                call abf_accumulator_add_data12(avg_values,avg_etot,pxi0(:))
            case(3)
                call abf_accumulator_add_data3(avg_values,avg_etot,pxi0(:))
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! write icf
        call abf_output_write_icf(avg_values,pxi0(:))
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

end subroutine abf_core_force

!===============================================================================
! Subroutine:  abf_core_force_numeric
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine abf_core_force_numeric()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accumulator
    use abf_output

    implicit none
    integer                :: i,j,m
    integer                :: ci
    real(PMFDP)            :: v,avg_etot
    ! --------------------------------------------------------------------------

    ! shift accuvalue history
    cvaluehist0(:) = cvaluehist1(:)

    ! save coordinate value to history
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvaluehist1(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    if( faccuepot ) then
        epothist1 = PotEne
    else
        epothist1 = 0.0d0
    end if

    ! shift ekin ene
    if( faccuekin ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] faccuekin cannot be on in abf_core_force_numeric!')
    end if

    ! calculate abf force to be applied -------------
    select case(feimode)
        case(1)
            call abf_accumulator_get_data1(cvaluehist1(:),la)
        case(2)
            call abf_accumulator_get_data2(cvaluehist1(:),la)
        case(3)
            call abf_accumulator_get_data3(cvaluehist1(:),la)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
    end select

    ! apply force filters
    if( .not. fapply_abf ) then
        ! we do not want to apply ABF force - set the forces to zero
        la(:) = 0.0d0
    end if

    ! project abf force along coordinate
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! calculate Z matrix and its inverse
    call abf_core_calc_Zmat

    ! pxip
    do i=1,NumOfABFCVs
        v = 0.0d0
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + CVContext%CVsDrvs(m,j,ci)*Vel(m,j)
            end do
        end do
        pxi0(i) = v
    end do

    ! mult by fzinv
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            v = v + fzinv(i,j) * pxi0(j)
        end do
        pxip(i) = v
    end do

    pxip(:) = pxip(:)

    ! update accumulator - we need at least two samples
    if( fstep .ge. 2 ) then
        ! final derivatives by central differences
        do i=1,NumOfABFCVs
            pxi0(i) = (pxip(i) - pxim(i))/fdtx - 0.5*(la(i) + pxi1(i))
        end do
        pxi1(:) = la(:)
        ! calculate coordinate values at time t-dt/2
        do i=1,NumOfABFCVs
            avg_values(i) = ABFCVList(i)%cv%get_average_value(cvaluehist0(i),cvaluehist1(i))
        end do

        avg_etot = 0.5d0*(epothist0 + epothist1) - ftotoffset

        ! add data to accumulator
        select case(feimode)
            case(1)
                call abf_accumulator_add_data12(avg_values,avg_etot,pxi0(:))
            case(2)
                call abf_accumulator_add_data12(avg_values,avg_etot,pxi0(:))
            case(3)
                call abf_accumulator_add_data3(avg_values,avg_etot,pxi0(:))
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! write icf
        call abf_output_write_icf(avg_values,pxi0(:))
    end if

    pxim(:) = pxip(:)

    return

end subroutine abf_core_force_numeric

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
