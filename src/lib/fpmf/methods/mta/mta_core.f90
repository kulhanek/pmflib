!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtc_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtc_core_main
!===============================================================================

subroutine mtc_core_main

    use mtc_restart
    use mtc_output
    use mtc_dat
    use mtc_accu
    use pmf_utils
    ! --------------------------------------------------------------------------

    call mtc_core_calc_Zmat
    call mtc_accu_add_data_online
    call mtc_output_write_output
    call mtc_restart_update

end subroutine mtc_core_main

!===============================================================================
! subroutine:  mtc_core_calc_Zmat
!===============================================================================

subroutine mtc_core_calc_Zmat

    use pmf_utils
    use mtc_dat
    use pmf_dat

    implicit none
    integer             :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfMTCCVs
        ci = MTCCVList(i)%cvindx
        do j=1,NumOfMTCCVs
            cj = MTCCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
        end do
    end do

    fzdet = 1.0d0

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfMTCCVs .gt. 1) then
        fzinv(:,:)  = fz(:,:)
        call dgetrf(NumOfMTCCVs,NumOfMTCCVs,fzinv,NumOfMTCCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[MTC] LU decomposition failed in mtc_core_calc_Zmat!')
        end if

         ! and finally determinant
        do i=1,NumOfMTCCVs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * fzinv(i,i)
            else
                fzdet = fzdet * fzinv(i,i)
            end if
        end do

        call dgetri(NumOfMTCCVs,fzinv,NumOfMTCCVs,indx,vv,NumOfMTCCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[MTC] Matrix inversion failed in mtc_core_calc_Zmat!')
        end if
    else
        fzdet       = fz(1,1)
        fzinv(1,1)  = 1.0d0/fz(1,1)
    end if

    return

end subroutine mtc_core_calc_Zmat

!===============================================================================

end module mtc_core
