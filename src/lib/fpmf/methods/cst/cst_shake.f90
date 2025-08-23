!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_shake

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_shake_calculate
!===============================================================================

subroutine cst_shake_calculate

    use pmf_utils
    use cst_dat
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_timers_start_timer(PMFLIB_CST_SHAKE_TIMER)

    select case(fshakesolver)
        case(CON_SHAKESOL_FM)
            call cst_shake_calculate_fm()   ! fixed shake: JAC(0,0)
        case(CON_SHAKESOL_MM)
            call cst_shake_calculate_mm()   ! mixed shake: JAC(0,P)
        case(CON_SHAKESOL_NM)
            call cst_shake_calculate_nm()   ! Newton-Raphson shake: JAC(P,P)
        case(CON_SHAKESOL_NMSVD)
            call cst_shake_calculate_nm_svd()   ! Newton-Raphson shake: JAC(P,P)
        case(CON_SHAKESOL_NMSVD_P)
            call cst_shake_calculate_nm_svd_p()   ! Newton-Raphson shake: JAC(P,P) + ContextP
        case(CON_SHAKESOL_DI)
            call cst_shake_calculate_di()   ! mixed shake: JAC(0,P) - diagonal solver
        case(CON_SHAKESOL_DIWG)
            call cst_shake_calculate_diwg() ! mixed shake: JAC(0,P) - diagonal solver with initial guess
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] SHAKE solver is not implemented in cst_shake_calculate!')
    end select

    call pmf_timers_stop_timer(PMFLIB_CST_SHAKE_TIMER)

end subroutine cst_shake_calculate

!===============================================================================
! Subroutine:  cst_shake_calculate_fm
!===============================================================================

subroutine cst_shake_calculate_fm

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,info,ci
    logical             :: done
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    lambdax(:) = 0.0d0

    ! calculate Jacobian matrix ------------------------
    call cst_shake_calc_jacobian_fm ! it calculates jac(0,0)

    if ( NumOfAllCONs .gt. 1 ) then
        ! LU decomposition
        indx(:) = 0
        call dgetrf(NumOfAllCONs,NumOfAllCONs,jac,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CST] LU decomposition failed in cst_shake_calculate_fm!')
        end if
    end if

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        if ( NumOfAllCONs .gt. 1 ) then
            ! solve LE
            call dgetrs('N',NumOfAllCONs,1,jac,NumOfAllCONs,indx,cv,NumOfAllCONs,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1, &
                                 '[CST] Solution of LE failed in cst_shake_calculate_fm!')
            end if
        else
            cv(1)=cv(1)/jac(1,1)
        end if

        ! correct lambda vector
        lambdax(:) = lambdax(:) + cv(:)

        ! calculate new position vector
        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_fm!')
    end if

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_fm

!===============================================================================
! Subroutine:  cst_shake_calculate_mm
!===============================================================================

subroutine cst_shake_calculate_mm

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,info,ci
    logical             :: done
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    lambdax(:) = 0.0d0

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        ! calculate Jacobian matrix
        call cst_shake_calc_jacobian_mm ! it calculates jac(0,P)

        if ( NumOfAllCONs .gt. 1 ) then
            ! LU decomposition
            indx(:) = 0
            call dgetrf(NumOfAllCONs,NumOfAllCONs,jac,NumOfAllCONs,indx,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                                 '[CST] LU decomposition failed in cst_shake_calculate_mm!')
            end if
            ! solve LE
            call dgetrs('N',NumOfAllCONs,1,jac,NumOfAllCONs,indx,cv,NumOfAllCONs,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1, &
                                 '[CST] Solution of LE failed in cst_shake_calculate_mm!')
            end if
        else
            cv(1)=cv(1)/jac(1,1)
        end if

        ! correct lambda vector
        lambdax(:) = lambdax(:) + cv(:)

        ! calculate new position vector
        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_mm!')
    end if

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_mm

!===============================================================================
! Subroutine:  cst_shake_calculate_nm
!===============================================================================

subroutine cst_shake_calculate_nm

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,info,ci
    logical             :: done
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    lambdax(:) = 0.0d0

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        ! calculate Jacobian matrix
        call cst_shake_calc_jacobian_nm ! it calculates jac(P,P)

        if ( NumOfAllCONs .gt. 1 ) then
            ! LU decomposition
            indx(:) = 0
            call dgetrf(NumOfAllCONs,NumOfAllCONs,jac,NumOfAllCONs,indx,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                                 '[CST] LU decomposition failed in cst_shake_calculate_nm!')
            end if
            ! solve LE
            call dgetrs('N',NumOfAllCONs,1,jac,NumOfAllCONs,indx,cv,NumOfAllCONs,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1, &
                                 '[CST] Solution of LE failed in cst_shake_calculate_nm!')
            end if
        else
            cv(1)=cv(1)/jac(1,1)
        end if

        ! correct lambda vector
        lambdax(:) = lambdax(:) + cv(:)

        ! calculate new position vector
        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_nm!')
    end if

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_nm

!===============================================================================
! Subroutine:  cst_shake_calculate_nm_svd
!===============================================================================

subroutine cst_shake_calculate_nm_svd

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,info,ci,orank
    logical             :: done
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    lambdax(:) = 0.0d0

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        ! calculate Jacobian matrix
        call cst_shake_calc_jacobian_nm ! it calculates jac(P,P)

        if ( NumOfAllCONs .gt. 1 ) then
            ! SVD decomposition
            call dgelss(NumOfAllCONs,NumOfAllCONs,1,jac,NumOfAllCONs,cv,NumOfAllCONs,vv,frcond,orank,work,lwork,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                                 '[CST] SVD decomposition failed in cst_calculate_lambda_nm_svd!')
            end if
        else
            cv(1)=cv(1)/jac(1,1)
        end if

        ! correct lambda vector
        lambdax(:) = lambdax(:) + cv(:)

        ! calculate new position vector
        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_nm_svd!')
    end if

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_nm_svd

!===============================================================================
! Subroutine:  cst_shake_calculate_nm_svd_P
!===============================================================================

subroutine cst_shake_calculate_nm_svd_P

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,info,ci,orank
    logical             :: done
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    lambdax(:) = 0.0d0

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        ! calculate Jacobian matrix
        call cst_shake_calc_jacobian_nm ! it calculates jac(P,P)

        if ( NumOfAllCONs .gt. 1 ) then
            ! SVD decomposition
            call dgelss(NumOfAllCONs,NumOfAllCONs,1,jac,NumOfAllCONs,cv,NumOfAllCONs,vv,frcond,orank,work,lwork,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                                 '[CST] SVD decomposition failed in cst_calculate_lambda_nm_svd!')
            end if
        else
            cv(1)=cv(1)/jac(1,1)
        end if

        ! correct lambda vector
        lambdax(:) = lambdax(:) + cv(:)

        ! calculate new position vector
        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContextP%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_nm_svd!')
    end if

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_nm_svd_P

!===============================================================================
! Subroutine:  cst_shake_calculate_di
! JAC - only diagonal elements are considered - no matrix algebra is necessary
!===============================================================================

subroutine cst_shake_calculate_di

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,ci
    real(PMFDP)         :: jacv
    logical             :: done
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    lambdax(:) = 0.0d0

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx

            ! calculate diagonal value
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                jacv = jacv - MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,ci))
            end do

            ! solve LE
            cv(i)=cv(i)/jacv

            ! correct lambda vector
            lambdax(i) = lambdax(i) + cv(i)

            ! calculate new position vector
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_di!')
    end if

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_di

!===============================================================================
! Subroutine:  cst_shake_calculate_diwg
! JAC - only diagonal elements are considered - no matrix algebra is necessary
! take lambda guess from the previous step
!===============================================================================

subroutine cst_shake_calculate_diwg

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,k,ci
    real(PMFDP)         :: jacv
    logical             :: done
    logical,save        :: initialized_lambda = .false. ! static
    real(PMFDP)         :: invn,dfsiter1,dfsiter2
    ! -----------------------------------------------------------------------------

    if( initialized_lambda ) then
        cv(:) = lambdax(:)
        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do
    else
        lambdax(:) = 0.0d0
    end if

! do step
    do fsiter=1,fmaxiter

        ! go through constraint list and calculate first derivative and constraint values at CrdP and cv
        call cst_constraints_calc_fdxp

        do i=1,NumOfAllCONs
            ci = CONList(i)%cvindx

            ! calculate diagonal value
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                jacv = jacv - MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,ci))
            end do

            ! solve LE
            cv(i)=cv(i)/jacv

            ! correct lambda vector
            lambdax(i) = lambdax(i) + cv(i)

            ! calculate new position vector
            do k=1,NumOfLAtoms
                CrdP(:,k) = CrdP(:,k) + MassInv(k)*cv(i)*CVContext%CVsDrvs(:,k,ci)
            end do
        end do

        ! check convergence criteria in lambdax
        done = .true.
        do i=1,NumOfAllCONs
            if( abs(cv(i)*isfdts) .gt. flambdatol ) done = .false.
        end do

        if( done ) exit

    end do

    if( fsiter .eq. fmaxiter ) then
        call pmf_utils_exit(PMF_OUT,1, &
                         '[CST] Maximum number of iterations in lambda calculation exceeded in cst_shake_calculate_di!')
    end if

    initialized_lambda = .true.

! update stats about iterations
    nsupdates = nsupdates + 1.0d0
    invn = 1.0d0 / nsupdates
    dfsiter1 = fsiter - mfsiter
    mfsiter  = mfsiter  + dfsiter1 * invn
    dfsiter2 = fsiter - mfsiter
    m2fsiter = m2fsiter + dfsiter1 * dfsiter2

end subroutine cst_shake_calculate_diwg

!===============================================================================
! Subroutine:  cst_shake_calc_jacobian_fm
!===============================================================================

subroutine cst_shake_calc_jacobian_fm

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
                jacv = jacv - MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            jac(i,j)=jacv
        end do
    end do

end subroutine cst_shake_calc_jacobian_fm

!===============================================================================
! Subroutine:  cst_shake_calc_jacobian_mm
!===============================================================================

subroutine cst_shake_calc_jacobian_mm

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
                jacv = jacv - MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,cj))
            end do
            jac(i,j)=jacv
        end do
    end do

end subroutine cst_shake_calc_jacobian_mm

!===============================================================================
! Subroutine:  cst_shake_calc_jacobian_nm
!===============================================================================

subroutine cst_shake_calc_jacobian_nm

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
                jacv = jacv - MassInv(k)*dot_product(CVContextP%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,cj))
            end do
            jac(i,j)=jacv
        end do
    end do

end subroutine cst_shake_calc_jacobian_nm

!===============================================================================

end module cst_shake

