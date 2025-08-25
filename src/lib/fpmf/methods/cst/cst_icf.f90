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

module cst_icf

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_core_calculate_icf
!===============================================================================

subroutine cst_icf_calculate_icf

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,j,k,m
    real(PMFDP)            :: f1,nv,v1,v2,dh
    ! --------------------------------------------------------------------------

    if( NumOfCONs .ne. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Only 1 CV supported in cst_core_calculate_icf!')
    end if

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

    ! start with dV/dx
    CSTFrc(:,:) = Frc(:,:)

    ! add constraint forces from SHAKE constraints only
    do i=NumOfCONs+1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do k=1,NumOfLAtoms
            CSTFrc(:,k) = CSTFrc(:,k) + lambda(i)*CVContext%CVsDrvs(:,k,ci)
        end do
    end do

! ICF-P
    i = 1   ! CV index
    ci = CONList(i)%cvindx
    f1 = 0.0d0
    nv = 0.0d0
    do j=1,CONList(i)%cv%natoms
        k = CONList(i)%cv%lindexes(j)
        do m=1,3
            ! force part
            nv = nv + CVContext%CVsDrvs(m,k,ci) * CVContext%CVsDrvs(m,k,ci)
            f1 = f1 + CVContext%CVsDrvs(m,k,ci) * CSTFrc(m,k)
        end do
    end do
    icfp(i) = - f1 / nv

    dh = 1e-5

! ICF-K by central differences
    do j=1,CONList(i)%cv%natoms
        k = CONList(i)%cv%lindexes(j)
        do m=1,3
            CSTFrc(:,:) = Crd(:,:)
            CSTFrc(m,k) = CSTFrc(m,k) + dh

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            call CVList(i)%cv%calculate_cv(CSTFrc,CVContextP)
            call calc_icfk_vec

            v1 = icfk_vec(m,k)

            ! write(*,*) 'v1 = ', v1

            CSTFrc(:,:) = Crd(:,:)
            CSTFrc(m,k) = CSTFrc(m,k) - dh

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            call CVList(i)%cv%calculate_cv(CSTFrc,CVContextP)
            call calc_icfk_vec

            v2 = icfk_vec(m,k)

          !  write(7894,*) v1, v2, (v1-v2)/(2.0d0 * dh)

            icfk(i) = icfk(i) + (v1-v2)/(2.0d0 * dh)
      end do
  end do

end subroutine cst_icf_calculate_icf

!===============================================================================
! Subroutine:  calc_icfk_vec
!===============================================================================

subroutine calc_icfk_vec

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,j,k,m
    real(PMFDP)            :: nv
    ! --------------------------------------------------------------------------

    i = 1   ! CV index
    ci = CONList(i)%cvindx
    nv = 0.0d0
    do k=1,NumOfLAtoms
        do m=1,3
            nv = nv + CVContextP%CVsDrvs(m,k,ci) * CVContextP%CVsDrvs(m,k,ci)
        end do
    end do

    ci = CONList(i)%cvindx
    do j=1,CONList(i)%cv%natoms
        k = CONList(i)%cv%lindexes(j)
        do m=1,3
            icfk_vec(m,k) = CVContextP%CVsDrvs(m,k,ci)/nv
        end do
    end do

end subroutine calc_icfk_vec

!===============================================================================

end module cst_icf

