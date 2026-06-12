!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025-2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_timers_start_timer(PMFLIB_CST_RATTLE_TIMER)

    select case(frattlesolver)
        case(CON_RATTLESOL_MA)
            call cst_rattlev_calculate_ma
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] RATTLE-V solver is not implemented in cst_rattlev_calculate!')
    end select

    call pmf_timers_stop_timer(PMFLIB_CST_RATTLE_TIMER)

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
    integer         :: i,ci,ki,k,m,info
    real(PMFDP)     :: tmp
    real(PMFDP)     :: invn,dfriter1,dfriter2
    ! --------------------------------------------------------------------------

    ! XP and VP

! get derivatives
    call cst_constraints_calc_fdxp

! construct right hand side
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        tmp = 0.0d0
        do ki=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(ki)
            do m=1,3
                tmp = tmp + CVContextP%CVsDrvs(m,k,ci) * VelP(m,k)
            end do
        end do
        lambdav(i) = - tmp
    end do

! left side
    call cst_constraints_calc_zmat_mw(CVContextP%CVsDRvs)

 ! solve LE
     if( NumOfAllCONs .gt. 1 ) then
        indx(:) = 0
        call dgetrf(NumOfAllCONs,NumOfAllCONs,zmat,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_rattlev_calculate_ma!')
        end if
        call dgetrs('N',NumOfAllCONs,1,zmat,NumOfAllCONs,indx,lambdav,NumOfAllCONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] Solution of LE failed in cst_rattlev_calculate_ma!')
        end if
     else
        lambdav(1) = lambdav(1) / zmat(1,1)
     end if

! correct velocities
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do ki=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(ki)
            VelP(:,k) = VelP(:,k) + lambdav(i)*MassInv(k)*CVContextP%CVsDrvs(:,k,ci)
        end do
    end do

! final check of convergence
    do i=1,NumOfAllCONs
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
    friter = 1.0d0
    nrupdates = nrupdates + 1.0d0
    invn = 1.0d0 / nrupdates
    dfriter1 = friter - mfriter
    mfriter  = mfriter  + dfriter1 * invn
    dfriter2 = friter - mfriter
    m2friter = m2friter + dfriter1 * dfriter2

end subroutine cst_rattlev_calculate_ma

!===============================================================================

end module cst_rattlev

