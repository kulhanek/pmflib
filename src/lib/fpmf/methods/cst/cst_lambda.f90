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

    call pmf_timers_start_timer(PMFLIB_CST_LAMBDA_TIMER)

    select case(ftds_lamsol)
        case(CON_LAMSOL_MD)
            ! noting to do
        case(CON_LAMSOL_V1)
            call cst_lambda_calculate_v1()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] LAMBDA solver (ftds_lamsol) is not implemented in cst_lambda_calculate!')
    end select

    call pmf_timers_stop_timer(PMFLIB_CST_LAMBDA_TIMER)

end subroutine cst_lambda_calculate

!===============================================================================
! Subroutine:  cst_lambda_calculate_v1
!===============================================================================

subroutine cst_lambda_calculate_v1

    use pmf_dat
    use pmf_utils
    use cst_dat
    use cst_constraints

    implicit none
    integer             :: i,j,k,m,ci,info
    real(PMFDP)         :: f1,k1,dx,v,lamp
    ! -----------------------------------------------------------------------------

    if( fstep - hist_len .le. 0 ) return

! rhs
    cv(:) = 0.0d0

! force part
    do i = 1,NumOfAllCONs
        ci = CONList(i)%cvindx
        lamp = 0.0d0
        do j=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(j)
            f1 = 0.0d0
            do m=1,3
                f1 = f1 - cvderhist(m,k,ci,hist_len+hist_fidx_tds) * frchist(m,k,hist_len+hist_fidx_tds)
            end do
            lamp = lamp + f1 * MassInv(k)
        end do
        cv(i) = lamp
    end do

! kinetic part
    ! velocity
    do k=1,NumOfLAtoms
        do m=1,3
            v = - 1.0d0 * velhist(m,k,hist_len+hist_fidx_tds-1) + 9.0d0 * velhist(m,k,hist_len+hist_fidx_tds+0) &
                + 9.0d0 * velhist(m,k,hist_len+hist_fidx_tds+1) - 1.0d0 * velhist(m,k,hist_len+hist_fidx_tds+2)
            ! v = v / 16.0d0 <- moved down
            TmpT(m,k) = v
         end do
    end do

    ! cvder time der
    do i = 1,NumOfAllCONs
        ci = CONList(i)%cvindx
        k1 = 0.0d0
        do j=1,CONList(i)%cv%natoms
            k = CONList(i)%cv%lindexes(j)
            do m=1,3
                dx =   1.0d0 * cvderhist(m,k,ci,hist_len+hist_fidx_tds-2) - 8.0d0 * cvderhist(m,k,ci,hist_len+hist_fidx_tds-1) &
                     + 8.0d0 * cvderhist(m,k,ci,hist_len+hist_fidx_tds+1) - 1.0d0 * cvderhist(m,k,ci,hist_len+hist_fidx_tds+2)
                k1 = k1 + dx * TmpT(m,k)
            end do
        end do
        ! one is for cvder time der, the other os for velocity
        cv(i) = cv(i) - k1 * ifdtx / 12.0d0 / 16.0d0
    end do

! zmat
    call cst_lambda_calc_zmat(hist_fidx_tds)

! linear equations
     if( NumOfAllCONs .gt. 1 ) then
        indx(:) = 0
        call dgetrf(NumOfAllCONs,NumOfAllCONs,zmata,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_lambda_calculate_v1!')
        end if
        call dgetrs('N',NumOfAllCONs,1,zmata,NumOfAllCONs,indx,cv,NumOfAllCONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] Solution of LE failed in cst_lambda_calculate_v1!')
        end if
     else
        cv(1) = cv(1) / zmata(1,1)
     end if

    lambdaEhist(:,hist_len+hist_fidx_tds) =  cv(:)

end subroutine cst_lambda_calculate_v1

!===============================================================================
! Subroutine:  cst_lambda_calc_zmat
!===============================================================================

subroutine cst_lambda_calc_zmat(fidx)

    use pmf_dat
    use cst_dat
    use cst_constraints

    implicit none
    integer         :: fidx
    ! --------------------------------------------
    integer         :: i,ci,j,cj,k
    real(PMFDP)     :: jacv
    ! --------------------------------------------------------------------------

    ! complete Jacobian matrix
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfAllCONs
            cj = CONList(j)%cvindx
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                jacv = jacv + MassInv(k)*dot_product(cvderhist(:,k,ci,hist_len+fidx),cvderhist(:,k,cj,hist_len+fidx))
            end do
            zmata(i,j)=jacv
        end do
    end do

end subroutine cst_lambda_calc_zmat

!===============================================================================

end module cst_lambda

