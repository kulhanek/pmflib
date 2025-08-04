!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

module cst_rattlev

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_rattlev_calculate
!===============================================================================

subroutine cst_rattlev_calculate

    use pmf_utils
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    select case(frattlesolver)
        case(CON_RATTLESOL_MA)
            call cst_rattlev_calculate_ma()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] RATTLE-V solver is not implemented in cst_rattlev_calculate!')
    end select

end subroutine cst_rattlev_calculate

!===============================================================================
! Subroutine:  cst_rattlev_calculate_ma
!===============================================================================

subroutine cst_rattlev_calculate_ma

    use pmf_utils
    use pmf_dat
    use cst_dat
    use cst_constraints

    implicit none
    integer         :: i,ci,k,info
    real(PMFDP)     :: tmp
    real(PMFDP)     :: invn,dfriter1,dfriter2
    ! --------------------------------------------------------------------------

    ! XP and VP

! get derivatives
    call cst_constraints_calc_fdxp

! construct right hand side
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        tmp = 0.0d0
        do k=1,NumOfLAtoms
            tmp = tmp + dot_product(CVContextP%CVsDrvs(:,k,ci),VelP(:,k))
        end do
        lambdav(i) = tmp
    end do

    ! DEBUG
    ! write(*,*) 'lambdav=',lambdav

! left side
    call cst_rattlev_calc_jacobian

    ! DEBUG
    ! write(*,*) 'jac=',jac

 ! solve LE
     if( NumOfCONs .gt. 1 ) then
        call dgetrf(NumOfCONs,NumOfCONs,jac,NumOfCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_rattlev_calculate_ma!')
        end if
        call dgetrs('N',NumOfCONs,1,jac,NumOfCONs,indx,lambdav,NumOfCONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] Solution of LE failed in cst_rattlev_calculate_ma!')
        end if
     else
        lambdav(1) = lambdav(1) / jac(1,1)
     end if

! correct velocities
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do k=1,NumOfLAtoms
            VelP(:,k) = VelP(:,k) + lambdav(i)*MassInv(k)*CVContextP%CVsDrvs(:,k,ci)
        end do
    end do

! transform to kcal/mol unit
    lambdav(:) = lambdav(:) * PMF_DT2VDT * PMF_L2CL / fdt ! FIXME

    write(PMF_DEBUG+fmytaskid,*) 'lambdav= ', lambdav(:)

! final check of convergence
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        tmp = 0.0d0
        do k=1,NumOfLAtoms
            tmp = tmp + dot_product(CVContextP%CVsDrvs(:,k,ci),VelP(:,k))
        end do
        if( abs(tmp) .gt. frveltol ) then
            write(PMF_OUT,*) 'RATTLE-V residual velocity: ', tmp
            call pmf_utils_exit(PMF_OUT,1,'[CST] RATTLE-V convergence was not achieved in cst_rattlev_calculate_ma!')
        end if
    end do

! update stats about iterations
    nrupdates = nrupdates + 1.0d0
    invn = 1.0d0 / nrupdates
    dfriter1 = friter - mfriter
    mfriter  = mfriter  + dfriter1 * invn
    dfriter2 = friter - mfriter
    m2friter = m2friter + dfriter1 * dfriter2

end subroutine cst_rattlev_calculate_ma

!===============================================================================
! Subroutine:  cst_rattlev_calc_jacobian
!===============================================================================

subroutine cst_rattlev_calc_jacobian

    use pmf_dat
    use cst_dat
    use cst_constraints

    implicit none
    integer                :: i,ci,j,cj,k
    real(PMFDP)            :: jacv
    ! --------------------------------------------------------------------------

    ! complete Jacobian matrix
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfCONs
            cj = CONList(j)%cvindx
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                jacv = jacv - MassInv(k)*dot_product(CVContextP%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,cj))
            end do
            jac(i,j)=jacv
        end do
    end do

end subroutine cst_rattlev_calc_jacobian

!===============================================================================

end module cst_rattlev

