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
! Subroutine:  cst_icf_calculate_icf
!===============================================================================

subroutine cst_icf_calculate_icf

    use pmf_utils
    use cst_dat
    use pmf_timers

    implicit none
    ! --------------------------------------------------------------------------

    call pmf_timers_start_timer(PMFLIB_CST_ICF_TIMER)

    select case(ftds_icfsol)
        case(CON_ICFSOL_V1)
            call cst_icf_calculate_v1()
        case(CON_ICFSOL_V2)
            call cst_icf_calculate_v2()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] ICF solver (ftds_icfsol) is not implemented in cst_icf_calculate_icf!')
    end select

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_TIMER)

end subroutine cst_icf_calculate_icf

!===============================================================================
! Subroutine:  cst_icf_calculate_v1
! numerical divergence
!===============================================================================

subroutine cst_icf_calculate_v1

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,l,cl,k,m
    real(PMFDP)            :: f1,v1,v2
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! ICFP part
    call cst_icf_calculate_zmatinv(CVContext)

    do i=1,NumOfCONs
        call cst_icf_calculate_vi(CVContext,i)
        f1 = 0.0d0
        do k=1,NumOfLAtoms
            do m=1,3
                f1 = f1 + icf_vi(m,k) * Frc(m,k)
            end do
        end do
        icfp(i) = - f1
    end do

! ICFK part

! ICF-K by central differences
    do i=1,NumOfCONs
        do k=1,NumOfLAtoms
            do m=1,3
                icf_he(:,:) = Crd(:,:)
                icf_he(m,k) = icf_he(m,k) + fpmf_div_dh

                CVContextP%CVsValues(:) = 0.0d0
                CVContextP%CVsDrvs(:,:,:) = 0.0d0
                do l=1,NumOfAllCONs
                    cl = CONList(l)%cvindx
                    call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
                end do
                call cst_icf_calculate_zmatinv(CVContextP)
                call cst_icf_calculate_vi(CVContextP,i)

                v1 = icf_vi(m,k)

                ! write(*,*) 'v1 = ', v1

                icf_he(:,:) = Crd(:,:)
                icf_he(m,k) = icf_he(m,k) - fpmf_div_dh

                CVContextP%CVsValues(:) = 0.0d0
                CVContextP%CVsDrvs(:,:,:) = 0.0d0

                do l=1,NumOfAllCONs
                    cl = CONList(l)%cvindx
                    call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
                end do
                call cst_icf_calculate_zmatinv(CVContextP)
                call cst_icf_calculate_vi(CVContextP,i)

                v2 = icf_vi(m,k)

              !  write(7894,*) v1, v2, (v1-v2)/(2.0d0 * dh)

                icfk(i) = icfk(i) + (v1-v2)/(2.0d0 * fpmf_div_dh)
          end do
      end do
  end do

end subroutine cst_icf_calculate_v1

!===============================================================================
! Subroutine:  cst_icf_calculate_v2
! analytical but with numerical/analytical second derivatives
!===============================================================================

subroutine cst_icf_calculate_v2

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: k,l,cl,m,cm,n,cn,o,ol,p
    real(PMFDP)            :: f1,v1
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! update CVs - calculate Values, gradients, and Hessians
    CVContext%CVsValues(:) = 0.0d0
    CVContext%CVsDrvs(:,:,:) = 0.0d0
    CVContext%CVs2ndDrvs(:,:,:,:,:) = 0.0d0

    do l=1,NumOfAllCONs
        call CONList(l)%cv%calculate_cv2ddrvs(Crd,CVContext)
    end do

! get inversion of W
    call cst_icf_calculate_zmatinv(CVContext)

! ICFP part
    do k=1,NumOfCONs
        call cst_icf_calculate_vi(CVContext,k)
        f1 = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                f1 = f1 + icf_vi(m,n) * Frc(m,n)
            end do
        end do
        icfp(k) = - f1
    end do

! ICF-K part
    do k=1,NumOfCONs
        ! simpler part :-)
        do l=1,NumOfAllCONs
            ! get Laplacian
            v1 = 0.0d0
            do ol=1,CONList(l)%cv%natoms
                o = CONList(l)%cv%lindexes(ol)
                do p=1,3
                    v1 = v1 + CVContext%CVs2ndDrvs(p,o,p,o,l)
                end do
            end do
            icfk(k) = icfk(k) + zmata(k,l) * v1
        end do

        ! harder part :-(
        do l=1,NumOfAllCONs
            cl = CONList(l)%cvindx
            do m=1,NumOfAllCONs
                cm = CONList(m)%cvindx
                do n=1,NumOfAllCONs
                    cn = CONList(n)%cvindx
                    v1 = cst_icf_calculate_wxi(cm,cn,cl)
                    icfk(k) = icfk(k) - zmata(k,m) * v1 * zmata(n,l)
                end do
            end do
        end do
  end do

end subroutine cst_icf_calculate_v2

!===============================================================================
! Subroutine:  cst_icf_calculate_vi
!===============================================================================

subroutine cst_icf_calculate_vi(ctx,i)

    use pmf_dat
    use cst_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i
    ! --------------------------------------------
    integer             :: j,cj,k,kj,m
    ! --------------------------------------------------------------------------

    icf_vi(:,:) = 0.0d0

    do j=1,NumOfAllCONs
        cj = CONList(j)%cvindx
        do kj=1,CONList(j)%cv%natoms
            k = CONList(j)%cv%lindexes(kj)
            do m=1,3
                icf_vi(m,k) = icf_vi(m,k) + zmata(i,j) * ctx%CVsDrvs(m,k,cj)
            end do
        end do
    end do

end subroutine cst_icf_calculate_vi

!===============================================================================
! Subroutine:  cst_icf_calculate_wxi
!===============================================================================

function cst_icf_calculate_wxi(cm,cn,cl) result(wxi)

    use pmf_dat
    use cst_dat

    implicit none
    integer             :: cm
    integer             :: cn
    integer             :: cl
    real(PMFDP)         :: wxi
    ! --------------------------------------------
    integer             :: d,o
    ! --------------------------------------------------------------------------

    wxi = 0.0d0

    d = 3*NumOfLAtoms

    call dgemv('N',d,d,1.0d0,CVContext%CVs2ndDrvs(:,:,:,:,cm),d,CVContext%CVsDrvs(:,:,cn),1,0.0d0,icf_vi,1)
    call dgemv('N',d,d,1.0d0,CVContext%CVs2ndDrvs(:,:,:,:,cn),d,CVContext%CVsDrvs(:,:,cm),1,1.0d0,icf_vi,1)

    do o=1,NumOfLAtoms
        wxi = wxi + dot_product(icf_vi(:,o),CVContext%CVsDrvs(:,o,cl))
    end do

end function cst_icf_calculate_wxi

!===============================================================================
! Subroutine:  cst_icf_calculate_zmatinv
!===============================================================================

subroutine cst_icf_calculate_zmatinv(ctx)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    type(CVContextType) :: ctx
    ! --------------------------------------------
    integer             :: i,ci,j,cj,k,info
    real(PMFDP)         :: jacv
    ! --------------------------------------------------------------------------

! this Z matrix is not mass weighted

! get the matrix
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfAllCONs
            cj = CONList(j)%cvindx
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                jacv = jacv + dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
            zmata(i,j) = jacv
        end do
    end do

! invert
    if ( NumOfAllCONs .gt. 1 ) then
        ! LU decomposition
        indx(:) = 0
        call dgetrf(NumOfAllCONs,NumOfAllCONs,zmata,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CST] LU decomposition failed in cst_icf_calculate_zmatinv!')
        end if

        ! invert
        call dgetri(NumOfAllCONs, zmata, NumOfAllCONs, indx, invwork, linvwork, info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                             '[CST] Matrix inversion failed in cst_icf_calculate_zmatinv!')
        end if
    else
        zmata(1,1) = 1.0d0/zmata(1,1)
    end if

end subroutine cst_icf_calculate_zmatinv

!===============================================================================

end module cst_icf

