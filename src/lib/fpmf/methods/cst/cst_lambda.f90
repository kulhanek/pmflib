!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module cst_lambda

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_lambda_calculate
!===============================================================================

subroutine cst_lambda_calculate

    use pmf_utils
    use cst_dat
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    ! call pmf_timers_start_timer(PMFLIB_CST_LAMBDA_TIMER)

    select case(flambdasolver)
        case(CON_LAMSOL_MD)
            ! noting to do
        case(CON_LAMSOL_SIMPLE)
            call cst_lambda_calculate_simple()
        case(CON_LAMSOL_FULL)
            call cst_lambda_calculate_full()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] LAMBDA solver is not implemented in cst_lambda_calculate!')
    end select

    ! call pmf_timers_stop_timer(PMFLIB_CST_LAMBDA_TIMER)

end subroutine cst_lambda_calculate

!===============================================================================
! Subroutine:  cst_lambda_calculate_simple
!===============================================================================

subroutine cst_lambda_calculate_simple

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,j,k,m,ci
    real(PMFDP)         :: f1,lp,lk1,lk2,k1,k2
    ! -----------------------------------------------------------------------------

    if ( NumOfAllCONs .ne. 1 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CST] Only one CV supported in cst_lambda_calculate_simple!')
    end if

! at t - force part
    i = 1   ! CV index
    ci = CONList(i)%cvindx
    f1 = 0.0d0
    do j=1,CONList(i)%cv%natoms
        k = CONList(i)%cv%lindexes(j)
        do m=1,3
            ! force part
            f1 = f1 + CVContext%CVsDrvs(m,k,ci) * Frc(m,k) * MassInv(k)
        end do
    end do
    lamphist(i,hist_len) = f1

! at t and t-dt/2 - kinetic part
    i = 1   ! CV index
    ci = CONList(i)%cvindx
    k1 = 0.0d0
    do j=1,CONList(i)%cv%natoms
        k = CONList(i)%cv%lindexes(j)
        do m=1,3
            ! force part
            k1 = k1 + CVContext%CVsDrvs(m,k,ci) * Vel(m,k)
        end do
    end do
    lamk1hist(i,hist_len) = k1

    cvderhist(:,:,:,hist_len) = CVContext%CVsDrvs(:,:,:)

! at t-dt and t-dt/2 - kinetic part
    i = 1   ! CV index
    ci = CONList(i)%cvindx
    k2 = 0.0d0
    do j=1,CONList(i)%cv%natoms
        k = CONList(i)%cv%lindexes(j)
        do m=1,3
            ! force part
            k2 = k2 + cvderhist(m,k,ci,hist_len-1) * Vel(m,k)
        end do
    end do
    lamk2hist(i,hist_len) = k2

! zmat
    call cst_shake_lambda_jacobian

! at t-dt
    lp = -lamphist(i,hist_len-1)
    lk1 = (lamk1hist(i,hist_len)   - lamk2hist(i,hist_len)) * ifdtx
    lk2 = (lamk1hist(i,hist_len-1) - lamk2hist(i,hist_len-1)) * ifdtx
    lambdaThist(i,hist_len-1) =  (lp - 0.5d0*(lk1 + lk2)) / jac(i,i)

end subroutine cst_lambda_calculate_simple

!===============================================================================
! Subroutine:  cst_lambda_calculate_full
!===============================================================================

subroutine cst_lambda_calculate_full

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,j,k,m,ci,info
    real(PMFDP)         :: f1,lp,lk1,lk2,k1,k2
    ! -----------------------------------------------------------------------------


! at t - force part
    do i = 1,NumOfAllCONs
        ci = CONList(i)%cvindx
        f1 = 0.0d0
        do j=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(j)
            do m=1,3
                ! force part
                f1 = f1 + CVContext%CVsDrvs(m,k,ci) * Frc(m,k) * MassInv(k)
            end do
        end do
        lamphist(i,hist_len) = f1
    end do

! at t and t-dt/2 - kinetic part
    do i = 1,NumOfAllCONs
        ci = CONList(i)%cvindx
        k1 = 0.0d0
        do j=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(j)
            do m=1,3
                ! force part
                k1 = k1 + CVContext%CVsDrvs(m,k,ci) * Vel(m,k)
            end do
        end do
        lamk1hist(i,hist_len) = k1
    end do

    cvderhist(:,:,:,hist_len) = CVContext%CVsDrvs(:,:,:)

! at t-dt and t-dt/2 - kinetic part
    do i = 1,NumOfAllCONs
        ci = CONList(i)%cvindx
        k2 = 0.0d0
        do j=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(j)
            do m=1,3
                ! force part
                k2 = k2 + cvderhist(m,k,ci,hist_len-1) * Vel(m,k)
            end do
        end do
        lamk2hist(i,hist_len) = k2
    end do

! at t-dt

! zmat
    call cst_shake_lambda_jacobian

! cv
    do i = 1,NumOfAllCONs
        lp = -lamphist(i,hist_len-1)
        lk1 = (lamk1hist(i,hist_len)   - lamk2hist(i,hist_len)) * ifdtx
        lk2 = (lamk1hist(i,hist_len-1) - lamk2hist(i,hist_len-1)) * ifdtx
        cv(i) = lp - 0.5d0*(lk1 + lk2)
    end do

! linear equations
     if( NumOfAllCONs .gt. 1 ) then
        indx(:) = 0
        call dgetrf(NumOfAllCONs,NumOfAllCONs,jac,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_rattlev_calculate_ma!')
        end if
        call dgetrs('N',NumOfAllCONs,1,jac,NumOfAllCONs,indx,cv,NumOfAllCONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] Solution of LE failed in cst_rattlev_calculate_ma!')
        end if
     else
        cv(1) = cv(1) / jac(1,1)
     end if

    lambdaThist(:,hist_len-1) =  cv(:)

end subroutine cst_lambda_calculate_full

!===============================================================================
! Subroutine:  cst_shake_lambda_jacobian
!===============================================================================

subroutine cst_shake_lambda_jacobian

    use pmf_dat
    use cst_dat
    use cst_constraints

    implicit none
    integer                :: i,ci,j,cj,k
    real(PMFDP)            :: jacv
    ! --------------------------------------------------------------------------

    ! complete Jacobian matrix
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfAllCONs
            cj = CONList(j)%cvindx
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                jacv = jacv + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            jac(i,j)=jacv
        end do
    end do

end subroutine cst_shake_lambda_jacobian

!===============================================================================

end module cst_lambda

