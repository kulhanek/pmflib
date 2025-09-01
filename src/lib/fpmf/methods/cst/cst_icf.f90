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

    ! if( .not. ( mod(faccustep,ftds_sample) .eq. 0 ) ) return

    call pmf_timers_start_timer(PMFLIB_CST_ICF_TIMER)

    select case(ftds_icfsol)
        case(CON_ICFSOL_V1)
            call cst_icf_calculate_v1()
            icfphist(:,hist_len) = icfp(:)
            icfkhist(:,hist_len) = icfk(:)
        case(CON_ICFSOL_V2)
            call cst_icf_calculate_v2()
            icfphist(:,hist_len) = icfp(:)
            icfkhist(:,hist_len) = icfk(:)
        case(CON_ICFSOL_V3)
            call cst_icf_calculate_v3()
            icfphist(:,hist_len) = icfp(:)
            icfkhist(:,hist_len) = icfk(:)
        case(CON_ICFSOL_V4)
            call cst_icf_calculate_v4()
            icfphist(:,hist_len) = icfp(:)
            icfkhist(:,hist_len) = icfk(:)
        case(CON_ICFSOL_V5)
            call cst_icf_calculate_v5()
            icfphist(:,hist_len) = icfp(:)
            icfkhist(:,hist_len) = icfk(:)
        case(CON_ICFSOL_V6)
            call cst_icf_calculate_v6()
        case(CON_ICFSOL_V7)
            call cst_icf_calculate_v7()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] ICF solver (ftds_icfsol) is not implemented in cst_icf_calculate_icf!')
    end select

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_TIMER)

end subroutine cst_icf_calculate_icf

!===============================================================================
! Subroutine:  cst_icf_calculate_shadow_H
!===============================================================================

subroutine cst_icf_calculate_shadow_H

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: i,k,m,ci
    real(PMFDP)            :: e1,e2,v,df
    ! --------------------------------------------------------------------------

    cfrchist(:,:,hist_len) = frchist(:,:,hist_len)
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        ! FIXME
        cfrchist(:,:,hist_len) = cfrchist(:,:,hist_len) + lambdaMhist(i,hist_len)*cvderhist(:,:,ci,hist_len)
    end do

    e1 = 0.0d0
    e2 = 0.0d0
    do k=1,NumOfLAtoms
        do m=1,3
            v  = - 1.0d0 * velhist(m,k,hist_len+hist_fidx_tds-1) + 9.0d0 * velhist(m,k,hist_len+hist_fidx_tds+0) &
                 + 9.0d0 * velhist(m,k,hist_len+hist_fidx_tds+1) - 1.0d0 * velhist(m,k,hist_len+hist_fidx_tds+2)
            v  = v / 16.0d0
            df = + 1.0d0 * cfrchist(m,k,hist_len+hist_fidx_tds-2) - 8.0d0 * cfrchist(m,k,hist_len+hist_fidx_tds-1) &
                 + 9.0d0 * cfrchist(m,k,hist_len+hist_fidx_tds+1) - 1.0d0 * cfrchist(m,k,hist_len+hist_fidx_tds+2)
            df = df / 12.0d0 * ifdtx
            e1 = e1 + 2.0d0 * v * df - MassInv(k) * (cfrchist(m,k,hist_len+hist_fidx_tds)*cfrchist(m,k,hist_len+hist_fidx_tds))
            e2 = e2 + df**2 * MassInv(k)
        end do
    end do

    shahist(hist_len+hist_fidx_tds) = e1 * fdtx ** 2 / 24.0d0 + e2 * fdtx**4 / 720.0d0

!    write(789,*) e1 * fdtx ** 2 / 24.0d0, e2 * fdtx**4 / 720.0d0

end subroutine cst_icf_calculate_shadow_H

!===============================================================================
! Subroutine:  cst_icf_calculate_v1
! numerical divergence
!===============================================================================

subroutine cst_icf_calculate_v1

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: i,l,cl,k,m
    real(PMFDP)            :: f1,v1,v2
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_icf_calculate_zmatinv(CVContext)

    do i=1,NumOfCONs
        call cst_icf_calculate_vi(CVContext,i,icf_vi1)
        f1 = 0.0d0
        do k=1,NumOfLAtoms
            do m=1,3
                f1 = f1 + icf_vi1(m,k) * Frc(m,k)
            end do
        end do
        icfp(i) = - f1
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
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
                call cst_icf_calculate_vi(CVContextP,i,icf_vi1)

                v1 = icf_vi1(m,k)

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
                call cst_icf_calculate_vi(CVContextP,i,icf_vi1)

                v2 = icf_vi1(m,k)

              !  write(7894,*) v1, v2, (v1-v2)/(2.0d0 * dh)

                icfk(i) = icfk(i) + (v1-v2)/(2.0d0 * fpmf_div_dh)
          end do
      end do
  end do
  call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v1

!===============================================================================
! Subroutine:  cst_icf_calculate_v2
! analytical but with numerical/analytical second derivatives
! optimized looping in ICFK
!===============================================================================

subroutine cst_icf_calculate_v2

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: k,l,cl,m,cm,n,cn,o,ol,p
    real(PMFDP)            :: f1,v1
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! update CVs - calculate Values, gradients, and Hessians
    call pmf_timers_start_timer(PMFLIB_CST_ICF_HESS_TIMER)
    CVContext%CVsValues(:) = 0.0d0
    CVContext%CVsDrvs(:,:,:) = 0.0d0
    CVContext%CVs2ndDrvs(:,:,:,:,:) = 0.0d0

    do l=1,NumOfAllCONs
        call CONList(l)%cv%calculate_cv2ddrvs(Crd,CVContext)
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_HESS_TIMER)

! get inversion of W
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)
    call cst_icf_calculate_zmatinv(CVContext)

! ICFP part
    do k=1,NumOfCONs
        call cst_icf_calculate_vi(CVContext,k,icf_vi1)
        f1 = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                f1 = f1 + icf_vi1(m,n) * Frc(m,n)
            end do
        end do
        icfp(k) = - f1
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICF-K part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
    do k=1,NumOfCONs
        ! simpler part :-)
        do l=1,NumOfAllCONs
            cl = CONList(l)%cvindx
            ! get Laplacian
            v1 = 0.0d0
            do ol=1,CONList(l)%cv%natoms
                o = CONList(l)%cv%lindexes(ol)
                do p=1,3
                    v1 = v1 + CVContext%CVs2ndDrvs(p,o,p,o,cl)
                end do
            end do
            icfk(k) = icfk(k) + zmat(k,l) * v1
        end do

        ! harder part :-(
        ! the loops are reorganized
        ! sparse matrix-vector multiplication in cst_icf_calculate_Higj
        do n=1,NumOfAllCONs
            cn = CONList(n)%cvindx
            icf_vi1(:,:) = 0.0d0
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                icf_vi1(:,:) = icf_vi1(:,:) + zmat(n,l) * CVContext%CVsDrvs(:,:,cl)
            end do
            do m=1,NumOfAllCONs
                cm = CONList(m)%cvindx
                icf_he(:,:) = 0.0d0
                call cst_icf_calculate_Higj(cn,cm)
                call cst_icf_calculate_Higj(cm,cn)
                v1 = 0.0d0
                do o=1,NumOfLAtoms
                    do p=1,3
                        v1 = v1 + icf_he(p,o) * icf_vi1(p,o)
                    end do
                end do
                icfk(k) = icfk(k) - zmat(k,m) * v1
            end do
        end do
  end do
  call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v2

!===============================================================================
! Subroutine:  cst_icf_calculate_v3
! analytical but with numerical/analytical second derivatives
! optimized looping in ICFK
! employ Hessian symmetry
!===============================================================================

subroutine cst_icf_calculate_v3

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: k,l,cl,m,cm,n,cn,o,ol,p,q
    real(PMFDP)            :: f1,v1,v2
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! update CVs - calculate Values, gradients, and Hessians
    call pmf_timers_start_timer(PMFLIB_CST_ICF_HESS_TIMER)
    CVContext%CVsValues(:) = 0.0d0
    CVContext%CVsDrvs(:,:,:) = 0.0d0
    CVContext%CVs2ndDrvs(:,:,:,:,:) = 0.0d0

    do l=1,NumOfAllCONs
        call CONList(l)%cv%calculate_cv2ddrvs(Crd,CVContext)
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_HESS_TIMER)

! get inversion of W
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)
    call cst_icf_calculate_zmatinv(CVContext)

! ICFP part
    do k=1,NumOfCONs
        call cst_icf_calculate_vi(CVContext,k,icf_vi1)
        f1 = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                f1 = f1 + icf_vi1(m,n) * Frc(m,n)
            end do
        end do
        icfp(k) = - f1
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICF-K part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)

    q = 0

    do k=1,NumOfCONs
        ! simpler part :-)
        do l=1,NumOfAllCONs
            cl = CONList(l)%cvindx
            ! get Laplacian
            v1 = 0.0d0
            do ol=1,CONList(l)%cv%natoms
                o = CONList(l)%cv%lindexes(ol)
                do p=1,3
                    v1 = v1 + CVContext%CVs2ndDrvs(p,o,p,o,cl)
                end do
            end do
            icfk(k) = icfk(k) + zmat(k,l) * v1
        end do

        ! harder part :-(
        ! the loops are reorganized
        ! sparse matrix-vector multiplication in cst_icf_calculate_Higj
        icf_vin(:,:,:) = 0.0d0
        do n=1,NumOfAllCONs
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                icf_vin(:,:,n) = icf_vin(:,:,n) + zmat(n,l) * CVContext%CVsDrvs(:,:,cl)
            end do
        end do
        do n=1,NumOfAllCONs
            cn = CONList(n)%cvindx
            do m=1,n                    ! triangular sum
                cm = CONList(m)%cvindx
                icf_he(:,:) = 0.0d0
                call cst_icf_calculate_Higj(cn,cm)
                call cst_icf_calculate_Higj(cm,cn)
                if( n .ne. m ) then
                    v1 = 0.0d0
                    v2 = 0.0d0
                    do o=1,NumOfLAtoms
                        do p=1,3
                            v1 = v1 + icf_he(p,o) * icf_vin(p,o,n)
                            v2 = v2 + icf_he(p,o) * icf_vin(p,o,m)
                        end do
                    end do
                    icfk(k) = icfk(k) - zmat(k,m) * v1 - zmat(k,n) * v2
                else
                    v1 = 0.0d0
                    do o=1,NumOfLAtoms
                        do p=1,3
                            v1 = v1 + icf_he(p,o) * icf_vin(p,o,n)
                        end do
                    end do
                    icfk(k) = icfk(k) - zmat(k,m) * v1
                end if
            end do
        end do
  end do

  call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v3

!===============================================================================
! Subroutine:  cst_icf_calculate_v4
! numerical divergence - stochastic “trace trick” for divergence
!===============================================================================

subroutine cst_icf_calculate_v4

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: k,s,l,cl
    real(PMFDP)            :: v1
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_icf_calculate_zmatll(CVContext%CVsDrvs)

    do k=1,NumOfCONs
        call cst_icf_calculate_vi_ll(CVContext%CVsDrvs,k,icf_vi1)
        icfp(k) = - sum( icf_vi1(:,:) * Frc(:,:) )
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
    do k=1,NumOfCONs

        ! generate z-probes
        call cst_icf_draw_probes_rademacher(sdiv_z)
        if( fpmf_sdiv_qr ) then
            call cst_icf_orthonormalize_probes(sdiv_z)
        end if

        v1 = 0.0d0
        do s=1,fpmf_sdiv_S

            icf_he(:,:) = Crd(:,:) + fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll(CVContextP%CVsDrvs,k,icf_vi1)

            icf_he(:,:) = Crd(:,:) - fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll(CVContextP%CVsDrvs,k,icf_vi2)

            v1 = v1 + sum( sdiv_z(:,:,s)*(icf_vi1(:,:) - icf_vi2(:,:)) )

        end do

        if( fpmf_sdiv_qr ) then
            icfk(k) =  3.0d0 * real(NumOfLAtoms,PMFDP) * v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        else
            icfk(k) =  v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        end if

    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v4

!===============================================================================
! Subroutine:  cst_icf_calculate_v5
! numerical divergence - stochastic “trace trick” for divergence
! mass weighted
!===============================================================================

subroutine cst_icf_calculate_v5

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: k,s,l,cl
    real(PMFDP)            :: v1
    ! --------------------------------------------------------------------------

    icfp(:) = 0.0d0
    icfk(:) = 0.0d0

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_icf_calculate_zmatll_mw(CVContext%CVsDrvs)

    do k=1,NumOfCONs
        call cst_icf_calculate_vi_ll_mw(CVContext%CVsDrvs,k,icf_vi1)
        icfp(k) = - sum( icf_vi1(:,:) * Frc(:,:) )
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
    do k=1,NumOfCONs

        ! generate z-probes
        call cst_icf_draw_probes_rademacher(sdiv_z)
        if( fpmf_sdiv_qr ) then
            call cst_icf_orthonormalize_probes(sdiv_z)
        end if

        v1 = 0.0d0
        do s=1,fpmf_sdiv_S

            icf_he(:,:) = Crd(:,:) + fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll_mw(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll_mw(CVContextP%CVsDrvs,k,icf_vi1)

            icf_he(:,:) = Crd(:,:) - fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll_mw(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll_mw(CVContextP%CVsDrvs,k,icf_vi2)

            v1 = v1 + sum( sdiv_z(:,:,s)*(icf_vi1(:,:) - icf_vi2(:,:)) )

        end do

        if( fpmf_sdiv_qr ) then
            icfk(k) =  3.0d0 * real(NumOfLAtoms,PMFDP) * v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        else
            icfk(k) =  v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        end if

    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v5

!===============================================================================
! Subroutine:  cst_icf_calculate_v6
! numerical divergence - stochastic “trace trick” for divergence
! mass weighted
!===============================================================================

subroutine cst_icf_calculate_v6

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: k,s,l,cl
    real(PMFDP)            :: v1
    ! --------------------------------------------------------------------------

    if( fstep - hist_len .le. 0 ) return

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_icf_calculate_zmatll_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds))

!    ! test orthogonality
!    do k=1,NumOfAllCONs
!        call cst_icf_calculate_vi_ll_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
!        do l=1,NumOfAllCONs
!            cl = CONList(l)%cvindx
!            v1 = sum( icf_vi1(:,:) * cvderhist(:,:,cl,hist_len+hist_fidx_tds))
!            write(78945,*) k,l,v1
!        end do
!    end do
!    stop

    do k=1,NumOfCONs
        call cst_icf_calculate_vi_ll_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
        icf_he(:,:) = ( -1.0d0 * frchist(:,:,hist_len+hist_fidx_tds-2) + 4.0d0 * frchist(:,:,hist_len+hist_fidx_tds-1) &
                        +4.0d0 * frchist(:,:,hist_len+hist_fidx_tds+1) - 1.0d0 * frchist(:,:,hist_len+hist_fidx_tds+2) ) &
                    / 6.0d0
        icfphist(k,hist_len+hist_fidx_tds) = - sum( icf_vi1(:,:) * icf_he(:,:) )
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
    do k=1,NumOfCONs

        ! generate z-probes
        call cst_icf_draw_probes_rademacher(sdiv_z)
        if( fpmf_sdiv_qr ) then
            call cst_icf_orthonormalize_probes(sdiv_z)
        end if

        v1 = 0.0d0
        do s=1,fpmf_sdiv_S

            icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds) + fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll_mw(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll_mw(CVContextP%CVsDrvs,k,icf_vi1)

            icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds) - fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll_mw(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll_mw(CVContextP%CVsDrvs,k,icf_vi2)

            v1 = v1 + sum( sdiv_z(:,:,s)*(icf_vi1(:,:) - icf_vi2(:,:)) )

        end do

        if( fpmf_sdiv_qr ) then
            icfkhist(k,hist_len+hist_fidx_tds) =  3.0d0 * real(NumOfLAtoms,PMFDP) * v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        else
            icfkhist(k,hist_len+hist_fidx_tds) =  v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        end if

    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v6

!===============================================================================
! Subroutine:  cst_icf_calculate_v7
! numerical divergence - stochastic “trace trick” for divergence
!===============================================================================

subroutine cst_icf_calculate_v7

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers

    implicit none
    integer                :: k,s,l,cl
    real(PMFDP)            :: v1
    ! --------------------------------------------------------------------------

    if( fstep - hist_len .le. 0 ) return

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_icf_calculate_zmatll(cvderhist(:,:,:,hist_len+hist_fidx_tds))

!    ! test orthogonality
!    do k=1,NumOfAllCONs
!        call cst_icf_calculate_vi_ll_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
!        do l=1,NumOfAllCONs
!            cl = CONList(l)%cvindx
!            v1 = sum( icf_vi1(:,:) * cvderhist(:,:,cl,hist_len+hist_fidx_tds))
!            write(78945,*) k,l,v1
!        end do
!    end do
!    stop

    do k=1,NumOfCONs
        call cst_icf_calculate_vi_ll(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
        icf_he(:,:) = ( -1.0d0 * frchist(:,:,hist_len+hist_fidx_tds-2) + 4.0d0 * frchist(:,:,hist_len+hist_fidx_tds-1) &
                        +4.0d0 * frchist(:,:,hist_len+hist_fidx_tds+1) - 1.0d0 * frchist(:,:,hist_len+hist_fidx_tds+2) ) &
                    / 6.0d0
        icfphist(k,hist_len+hist_fidx_tds) = - sum( icf_vi1(:,:) * icf_he(:,:) )
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
    do k=1,NumOfCONs

        ! generate z-probes
        call cst_icf_draw_probes_rademacher(sdiv_z)
        if( fpmf_sdiv_qr ) then
            call cst_icf_orthonormalize_probes(sdiv_z)
        end if

        v1 = 0.0d0
        do s=1,fpmf_sdiv_S

            icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds) + fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll(CVContextP%CVsDrvs,k,icf_vi1)

            icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds) - fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_icf_calculate_zmatll(CVContextP%CVsDrvs)
            call cst_icf_calculate_vi_ll(CVContextP%CVsDrvs,k,icf_vi2)

            v1 = v1 + sum( sdiv_z(:,:,s)*(icf_vi1(:,:) - icf_vi2(:,:)) )

        end do

        if( fpmf_sdiv_qr ) then
            icfkhist(k,hist_len+hist_fidx_tds) =  3.0d0 * real(NumOfLAtoms,PMFDP) * v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        else
            icfkhist(k,hist_len+hist_fidx_tds) =  v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        end if

    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v7

!===============================================================================
! Subroutine:  cst_icf_draw_probes_rademacher
!===============================================================================

subroutine cst_icf_draw_probes_rademacher(z)

    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)     :: z(:,:,:)
    ! --------------------------------------------
    integer         :: i,j,k
    real(PMFDP)     :: u
    ! --------------------------------------------------------------------------

    ! dense Rademacher vector: entries are +/-1 with prob 1/2.

    do j = 1, fpmf_sdiv_S
        do i = 1, NumOfLAtoms
            do k=1,3
                call random_number(u)    ! u ~ U(0,1)
                z(k,i,j) = merge( 1.0d0, -1.0d0, u >= 0.5d0)
            end do
        end do
    end do

end subroutine cst_icf_draw_probes_rademacher

!===============================================================================
! Subroutine:  cst_icf_orthonormalize_probes
!===============================================================================

subroutine cst_icf_orthonormalize_probes(z)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)     :: z(:,:,:)
    ! --------------------------------------------
    integer         :: n,k,info
    ! --------------------------------------------

    n = 3*NumOfLAtoms
    k = min(3*NumOfLAtoms,fpmf_sdiv_S)

    call dgeqrf(n,fpmf_sdiv_S,z,n,sdivtau,sdivwork,lsdivwork,info)
    if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CST] QR decomposition failed in cst_icf_orthonormalize_probes!')
    end if

    call dorgqr(n,fpmf_sdiv_S,k,z,n,sdivtau,sdivwork,lsdivwork,info)
    if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CST] QR orthonormalization failed in cst_icf_orthonormalize_probes!')
    end if

end subroutine cst_icf_orthonormalize_probes

!===============================================================================
! Subroutine:  cst_icf_calculate_vi
!===============================================================================

subroutine cst_icf_calculate_vi(ctx,i,icf_vi)

    use pmf_dat
    use cst_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i
    real(PMFDP)         :: icf_vi(:,:)
    ! --------------------------------------------
    integer             :: j,cj,k,kj,m
    ! --------------------------------------------------------------------------

    icf_vi(:,:) = 0.0d0

    do j=1,NumOfAllCONs
        cj = CONList(j)%cvindx
        do kj=1,CONList(j)%cv%natoms
            k = CONList(j)%cv%lindexes(kj)
            do m=1,3
                icf_vi(m,k) = icf_vi(m,k) + zmat(i,j) * ctx%CVsDrvs(m,k,cj)
            end do
        end do
    end do

end subroutine cst_icf_calculate_vi

!===============================================================================
! Subroutine:  cst_icf_calculate_vi_ll
! optimized
!===============================================================================

subroutine cst_icf_calculate_vi_ll(cvsdrvs,i,icf_vi)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)         :: cvsdrvs(:,:,:)
    integer             :: i
    real(PMFDP)         :: icf_vi(:,:)
    ! --------------------------------------------
    integer             :: j,cj,k,kj,m,info
    ! --------------------------------------------------------------------------

    ! form e_i
    cv(:) = 0.0d0
    cv(i) = 1.0d0

    ! solve ZMATA w = e_i
    call dpotrs('L',NumOfAllCONs,1,zmat,NumOfAllCONs,cv,NumOfAllCONs,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] LL linear equation failed in cst_icf_calculate_vi_ll!')
    end if

    icf_vi(:,:) = 0.0d0

    ! sparse version
    do j=1,NumOfAllCONs
        cj = CONList(j)%cvindx
        do kj=1,CONList(j)%cv%natoms
            k = CONList(j)%cv%lindexes(kj)
            do m=1,3
                icf_vi(m,k) = icf_vi(m,k) + cv(j) * cvsdrvs(m,k,cj)
            end do
        end do
    end do

end subroutine cst_icf_calculate_vi_ll

!===============================================================================
! Subroutine:  cst_icf_calculate_vi_ll_mw
! optimized
!===============================================================================

subroutine cst_icf_calculate_vi_ll_mw(cvsdrvs,i,icf_vi)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)         :: cvsdrvs(:,:,:)
    integer             :: i
    real(PMFDP)         :: icf_vi(:,:)
    ! --------------------------------------------
    integer             :: j,cj,k,kj,m,info
    ! --------------------------------------------------------------------------

    ! form e_i
    cv(:) = 0.0d0
    cv(i) = 1.0d0

    ! solve ZMATA w = e_i
    call dpotrs('L',NumOfAllCONs,1,zmat,NumOfAllCONs,cv,NumOfAllCONs,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] LL linear equation failed in cst_icf_calculate_vi_ll!')
    end if

    icf_vi(:,:) = 0.0d0

    ! sparse version
    do j=1,NumOfAllCONs
        cj = CONList(j)%cvindx
        do kj=1,CONList(j)%cv%natoms
            k = CONList(j)%cv%lindexes(kj)
            do m=1,3
                icf_vi(m,k) = icf_vi(m,k) + MassInv(k) * cv(j) * cvsdrvs(m,k,cj)
            end do
        end do
    end do

end subroutine cst_icf_calculate_vi_ll_mw

!===============================================================================
! Subroutine:  cst_icf_calculate_Higj
!===============================================================================

subroutine cst_icf_calculate_Higj(ci,cj)

    use pmf_dat
    use cst_dat

    implicit none
    integer             :: ci
    integer             :: cj
    ! --------------------------------------------
    integer             :: ki,k,o,lj,l,p
    real(PMFDP)         :: v1
    ! --------------------------------------------------------------------------

    do ki=1,CVList(ci)%cv%natoms
        k = CVList(ci)%cv%lindexes(ki)
        do o=1,3
            v1 = 0.0d0
            do lj=1,CVList(cj)%cv%natoms
                l = CVList(cj)%cv%lindexes(lj)
                do p=1,3
                    v1 = v1 + CVContext%CVs2ndDrvs(p,l,o,k,ci) * CVContext%CVsDrvs(p,l,cj)
                end do
            end do
            icf_he(o,k) = icf_he(o,k) + v1
        end do
    end do

end subroutine cst_icf_calculate_Higj

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
            zmat(i,j) = jacv
        end do
    end do

! invert
    if ( NumOfAllCONs .gt. 1 ) then
        ! LU decomposition
        indx(:) = 0
        call dgetrf(NumOfAllCONs,NumOfAllCONs,zmat,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                             '[CST] LU decomposition failed in cst_icf_calculate_zmatinv!')
        end if

        ! invert
        call dgetri(NumOfAllCONs, zmat, NumOfAllCONs, indx, invwork, linvwork, info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1, &
                             '[CST] Matrix inversion failed in cst_icf_calculate_zmatinv!')
        end if
    else
        zmat(1,1) = 1.0d0/zmat(1,1)
    end if

end subroutine cst_icf_calculate_zmatinv

!===============================================================================
! Subroutine:  cst_icf_calculate_zmatll
! optimized
!===============================================================================

subroutine cst_icf_calculate_zmatll(cvsdrvs)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)     :: cvsdrvs(:,:,:)
    ! --------------------------------------------
    integer         :: i,ci,j,cj,info,k,m
    real(PMFDP)     :: jacv
    ! --------------------------------------------------------------------------

! this Z matrix is not mass weighted

! get the matrix
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do j=1,i
            cj = CONList(j)%cvindx
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                do m=1,3
                    jacv = jacv + cvsdrvs(m,k,ci)*cvsdrvs(m,k,cj)
                end do
            end do
            zmat(i,j) = jacv
            zmat(j,i) = jacv
        end do
    end do

! calc LL
    call dpotrf('L',NumOfAllCONs,zmat,NumOfAllCONs,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] LL decomposition failed in cst_icf_calculate_zmatll!')
    end if

end subroutine cst_icf_calculate_zmatll

!===============================================================================
! Subroutine:  cst_icf_calculate_zmatll
! optimized
!===============================================================================

subroutine cst_icf_calculate_zmatll_mw(cvsdrvs)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)             :: cvsdrvs(:,:,:)
    ! --------------------------------------------
    integer                 :: i,ci,j,cj,info,k,m
    real(PMFDP)             :: jacv,v1
    ! --------------------------------------------------------------------------

! get the matrix
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do j=1,i
            cj = CONList(j)%cvindx
            jacv = 0.0d0
            do k=1,NumOfLAtoms
                v1 = 0.0
                do m=1,3
                    v1 = v1 + cvsdrvs(m,k,ci)*cvsdrvs(m,k,cj)
                end do
                jacv = jacv + MassInv(k)*v1
            end do
            zmat(i,j) = jacv
            zmat(j,i) = jacv
        end do
    end do

! calc LL
    call dpotrf('L',NumOfAllCONs,zmat,NumOfAllCONs,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] LL decomposition failed in cst_icf_calculate_zmatll_mw!')
    end if

end subroutine cst_icf_calculate_zmatll_mw

!===============================================================================

end module cst_icf

