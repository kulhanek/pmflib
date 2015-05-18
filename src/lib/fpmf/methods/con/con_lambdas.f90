!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module con_lambdas

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  con_lambdas_calculate
!===============================================================================

subroutine con_lambdas_calculate

    use pmf_utils
    use con_dat

    implicit none
    ! --------------------------------------------------------------------------

    select case(flambdasolver)
    case(CON_LS_NM)
        call con_lambdas_calculate_nm(fliter)
    case(CON_LS_CM)
        call con_lambdas_calculate_cm(fliter)
    case default
        call pmf_utils_exit(PMF_OUT,1,'[CON] solver is not implemented!')
    end select

end subroutine con_lambdas_calculate

!===============================================================================
! Subroutine:  con_lambdas_calculate_nm
!===============================================================================

subroutine con_lambdas_calculate_nm(iter)

 use pmf_dat
 use pmf_utils
 use con_dat
 use con_constraints

 implicit none
 integer            :: iter  ! number of iteration to achieve flambdatol
 ! -----------------------------------------------
 integer            :: i,k,info,ci
 real(PMFDP)        :: isfdt
 logical            :: done
 ! -----------------------------------------------------------------------------

 lambda(:) = 0.0d0

 select case(fintalg)
    case(IA_LEAP_FROG)
        isfdt = 1.0d0/(fdt*fdt) * PMF_L2CL
    case(IA_VEL_VERLET)
        isfdt = 2.0d0/(fdt*fdt) * PMF_L2CL
    case default
        call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm!')
 end select

 ! do Newton step -----------------------------------------------------------
 do iter=1,fmaxiter

    ! calculate Jacobian matrix ------------------------
    call con_lambdas_calc_jacobian ! it calculates jac and cv

    ! DEBUG
    ! write(*,*) 'A', iter, cv(1)

    if ( NumOfCONs .gt. 1 ) then
        ! LU decomposition
        indx(:) = 0
        call dgetrf(NumOfCONs,NumOfCONs,jac,NumOfCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CON] LU decomposition failed in con_calculate_lambda_nm!')
        end if
        call dgetrs('N',NumOfCONs,1,jac,NumOfCONs,indx,cv,NumOfCONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                             '[CON] Solution of LE failed in con_calculate_lambda_nm!')
        end if
    else
        cv(1)=cv(1)/jac(1,1)
    end if

    ! DEBUG
    ! write(*,*) 'B', iter, cv(1)

    ! correct lambda vector -------------------------
    lambda = lambda + cv

    ! calculate new position vector
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do k=1,NumOfLAtoms
            CrdP(:,k) = CrdP(:,k) - MassInv(k)*cv(i)*CVContextP%CVsDrvs(:,k,ci)
        end do
    end do

    ! check convergence criteria in lambdax
    done = .true.
    do i=1,NumOfCONs
        if( abs(cv(i)*isfdt) .gt. flambdatol ) done = .false.
    end do

    if( done ) exit

 end do

 ! DEBUG
 !write(*,*)

 do i=1,NumOfCONs
    lambda(i) = lambda(i)*isfdt
 end do

 ! final derivatives and values
 call con_constraints_calc_fdxp

 if( iter .eq. fmaxiter ) then
    call pmf_utils_exit(PMF_OUT,1, &
                     '[CON] Maximum number of iterations in lambda calculation exceeded!')
 end if 

 return

end subroutine con_lambdas_calculate_nm

!===============================================================================
! Subroutine:  con_lambdas_calculate_cm
!===============================================================================

subroutine con_lambdas_calculate_cm(iter)

 use pmf_utils
 use pmf_dat
 use con_dat
 use con_constraints

 implicit none
 integer            :: iter  ! number of iteration to achieve flambdatol
 ! -----------------------------------------------
 integer            :: i,k,info,ci
 real(PMFDP)        :: isfdt
 logical            :: done
 ! -----------------------------------------------------------------------------

 lambda(:) = 0.0d0

 select case(fintalg)
    case(IA_LEAP_FROG)
        isfdt = 1.0d0/(fdt*fdt) * PMF_L2CL
    case(IA_VEL_VERLET)
        isfdt = 2.0d0/(fdt*fdt) * PMF_L2CL
    case default
        call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm!')
 end select

 ! calculate Jacobian matrix ------------------------
 call con_lambdas_calc_jacobian ! it calculates jac and cv

 if ( NumOfCONs .gt. 1 ) then
    ! LU decomposition
    indx(:) = 0
    call dgetrf(NumOfCONs,NumOfCONs,jac,NumOfCONs,indx,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                            '[CON] LU decomposition failed in con_calculate_lambda_cm!')
    end if
 end if

 ! do Newton step -----------------------------------------------------------
 do iter=1,fmaxiter

    if ( NumOfCONs .gt. 1 ) then
        call dgetrs('N',NumOfCONs,1,jac,NumOfCONs,indx,cv,NumOfCONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                             '[CON] Solution of LE failed in con_calculate_lambda_cm!')
        end if
    else
        cv(1)=cv(1)/jac(1,1)
    end if

    ! correct lambda vector -------------------------
    lambda = lambda + cv

    ! calculate new position vector
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do k=1,NumOfLAtoms
            CrdP(:,k) = CrdP(:,k) - MassInv(k)*cv(i)*CVContextP%CVsDrvs(:,k,ci)
        end do
    end do

    ! check convergence criteria in lambdax
    done = .true.
    do i=1,NumOfCONs
        if( abs(cv(i)*isfdt) .gt. flambdatol ) done = .false.
    end do

    if( done ) exit

 end do

 do i=1,NumOfCONs
    lambda(i) = lambda(i)*isfdt
 end do

 ! final derivatives
 call con_constraints_calc_fdxp

 if( iter .ge. fmaxiter ) then
    call pmf_utils_exit(PMF_OUT,1, &
                     '[CON] Maximum number of iterations in lambda calculation exceeded!')
 end if 

 return

end subroutine con_lambdas_calculate_cm

!===============================================================================
! Subroutine:  con_lambdas_calc_jacobian
!===============================================================================

subroutine con_lambdas_calc_jacobian

 use pmf_dat
 use con_dat
 use con_constraints

 implicit none
 integer                :: i,ci,j,cj,k
 real(PMFDP)            :: jacv
 ! ------------------------------------------------------------------------------

 ! go throught constraint list and calculate first derivative matrixes and constraint values
 call con_constraints_calc_fdxp

 ! complete Jacobian matrix
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx
    do j=1,NumOfCONs
        cj = CONList(j)%cvindx
        jacv = 0.0d0
        do k=1,NumOfLAtoms
            jacv = jacv - MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
        end do
        jac(i,j)=jacv
    end do
 end do

end subroutine con_lambdas_calc_jacobian

!===============================================================================

end module con_lambdas

