!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025-2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    integer     :: faccustep_test
    ! --------------------------------------------------------------------------

    ! this is not optimal but better than nothing
    faccustep_test = faccustep
    if( .not. (fintcalc .or. ftdscalc) ) return
    if( enevalidhist(hist_len+hist_fidx_tds) ) faccustep_test = faccustep_test + 1
    if( .not. ( (mod(faccustep_test,ftds_sample) .eq. 0) .and. enevalidhist(hist_len+hist_fidx_tds) ) ) return

    call pmf_timers_start_timer(PMFLIB_CST_ICF_TIMER)

    select case(ftds_icfsol)
        case(CON_ICFSOL_V1)
            call cst_icf_calculate_v1()
        case(CON_ICFSOL_V1MW)
            call cst_icf_calculate_v1_mw()
        case(CON_ICFSOL_V2)
            call cst_icf_calculate_v2()
        case(CON_ICFSOL_V2MW)
            call cst_icf_calculate_v2_mw()
        case(CON_ICFSOL_V3)
            call cst_icf_calculate_v3()
        case(CON_ICFSOL_V3MW)
            call cst_icf_calculate_v3_mw()
        case default
            call pmf_utils_exit(PMF_OUT,1,'[CST] ICF solver (ftds_icfsol) is not implemented in cst_icf_calculate_icf!')
    end select

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_TIMER)

end subroutine cst_icf_calculate_icf

!===============================================================================
! Subroutine:  cst_icf_calculate_v1
! numerical divergence, no mass-weighted
!===============================================================================

subroutine cst_icf_calculate_v1

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers
    use cst_constraints

    implicit none
    integer                :: k,l,cl,n,m
    real(PMFDP)            :: icfp,icfk
    ! --------------------------------------------------------------------------

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_constraints_calc_zmat(cvderhist(:,:,:,hist_len+hist_fidx_tds))
    call cst_icf_calculate_zmatll

    do k=1,NumOfCONs
        call cst_icf_calculate_vill(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
        icfp = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icfp = icfp + icf_vi1(m,n) * frchist(m,n,hist_len+hist_fidx_tds)
            end do
        end do
        icfphist(k,hist_len+hist_fidx_tds) = - icfp   ! FRC is force not gradient
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
! ICF-K by central differences
    do k=1,NumOfCONs
        icfk = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds)
                icf_he(m,n) = icf_he(m,n) + fpmf_div_dh

                CVContextP%CVsValues(:) = 0.0d0
                CVContextP%CVsDrvs(:,:,:) = 0.0d0
                do l=1,NumOfAllCONs
                    cl = CONList(l)%cvindx
                    call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
                end do
                call cst_constraints_calc_zmat(CVContextP%CVsDrvs)
                call cst_icf_calculate_zmatll
                call cst_icf_calculate_vill(CVContextP%CVsDrvs,k,icf_vi1)

                icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds)
                icf_he(m,n) = icf_he(m,n) - fpmf_div_dh

                CVContextP%CVsValues(:) = 0.0d0
                CVContextP%CVsDrvs(:,:,:) = 0.0d0

                do l=1,NumOfAllCONs
                    cl = CONList(l)%cvindx
                    call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
                end do
                call cst_constraints_calc_zmat(CVContextP%CVsDrvs)
                call cst_icf_calculate_zmatll
                call cst_icf_calculate_vill(CVContextP%CVsDrvs,k,icf_vi2)

                icfk = icfk + (icf_vi1(m,n)-icf_vi2(m,n))/(2.0d0 * fpmf_div_dh)
            end do
        end do
        icfkhist(k,hist_len+hist_fidx_tds) = icfk
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v1

!===============================================================================
! Subroutine:  cst_icf_calculate_v1_mw
! numerical divergence, mass-weighted
!===============================================================================

subroutine cst_icf_calculate_v1_mw

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers
    use cst_constraints

    implicit none
    integer                :: k,l,cl,n,m
    real(PMFDP)            :: icfp,icfk
    ! --------------------------------------------------------------------------

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_constraints_calc_zmat_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds))
    call cst_icf_calculate_zmatll

    do k=1,NumOfCONs
        call cst_icf_calculate_vill_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
        icfp = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icfp = icfp + icf_vi1(m,n) * frchist(m,n,hist_len+hist_fidx_tds)
            end do
        end do
        icfphist(k,hist_len+hist_fidx_tds) = - icfp   ! FRC is force not gradient
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICFK part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)
! ICF-K by central differences
    do k=1,NumOfCONs
        icfk = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds)
                icf_he(m,n) = icf_he(m,n) + fpmf_div_dh

                CVContextP%CVsValues(:) = 0.0d0
                CVContextP%CVsDrvs(:,:,:) = 0.0d0
                do l=1,NumOfAllCONs
                    cl = CONList(l)%cvindx
                    call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
                end do
                call cst_constraints_calc_zmat_mw(CVContextP%CVsDrvs)
                call cst_icf_calculate_zmatll
                call cst_icf_calculate_vill_mw(CVContextP%CVsDrvs,k,icf_vi1)

                icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds)
                icf_he(m,n) = icf_he(m,n) - fpmf_div_dh

                CVContextP%CVsValues(:) = 0.0d0
                CVContextP%CVsDrvs(:,:,:) = 0.0d0

                do l=1,NumOfAllCONs
                    cl = CONList(l)%cvindx
                    call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
                end do
                call cst_constraints_calc_zmat_mw(CVContextP%CVsDrvs)
                call cst_icf_calculate_zmatll
                call cst_icf_calculate_vill_mw(CVContextP%CVsDrvs,k,icf_vi2)

                icfk = icfk + (icf_vi1(m,n)-icf_vi2(m,n))/(2.0d0 * fpmf_div_dh)
            end do
        end do
        icfkhist(k,hist_len+hist_fidx_tds) = icfk
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v1_mw


!===============================================================================
! Subroutine:  cst_icf_calculate_v2
! analytical but with numerical/analytical second derivatives
! optimized looping in ICFK, employ Hessian symmetry
!===============================================================================

subroutine cst_icf_calculate_v2

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers
    use cst_constraints

    implicit none
    integer                :: k,l,cl,m,cm,n,cn,o,ol,p
    real(PMFDP)            :: v1,v2,icfp,icfk
    ! --------------------------------------------------------------------------

! update CVs - calculate Values, gradients, and Hessians
    call pmf_timers_start_timer(PMFLIB_CST_ICF_HESS_TIMER)
    CVContextP%CVsValues(:) = 0.0d0
    CVContextP%CVsDrvs(:,:,:) = 0.0d0
    CVContextP%CVs2ndDrvs(:,:,:,:,:) = 0.0d0

    do l=1,NumOfAllCONs
        call CONList(l)%cv%calculate_cv2ddrvs(crdhist(:,:,hist_len+hist_fidx_tds),CVContextP)
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_HESS_TIMER)

! get inversion of W
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_constraints_calc_zmat(CVContextP%CVsDrvs)
    call cst_icf_calculate_zmatinv

! ICFP part
    do k=1,NumOfCONs
        call cst_icf_calculate_vi(CVContextP%CVsDrvs,k,icf_vi1)
        icfp = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icfp = icfp + icf_vi1(m,n) * frchist(m,n,hist_len+hist_fidx_tds)
            end do
        end do
        icfphist(k,hist_len+hist_fidx_tds) = - icfp   ! FRC is force not gradient
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICF-K part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)

    do k=1,NumOfCONs
        icfk = 0.0d0

        ! simpler part :-)
        do l=1,NumOfAllCONs
            cl = CONList(l)%cvindx
            ! get Laplacian
            v1 = 0.0d0
            do ol=1,CONList(l)%cv%natoms
                o = CONList(l)%cv%lindexes(ol)
                do p=1,3
                    v1 = v1 + CVContextP%CVs2ndDrvs(p,o,p,o,cl)
                end do
            end do
            icfk = icfk + zmat(k,l) * v1
        end do

        ! harder part :-(
        ! the loops are reorganized
        ! sparse matrix-vector multiplication in cst_icf_calculate_Higj
        icf_vin(:,:,:) = 0.0d0
        do n=1,NumOfAllCONs
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                icf_vin(:,:,n) = icf_vin(:,:,n) + zmat(n,l) * CVContextP%CVsDrvs(:,:,cl)
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
                    icfk = icfk - zmat(k,m) * v1 - zmat(k,n) * v2
                else
                    v1 = 0.0d0
                    do o=1,NumOfLAtoms
                        do p=1,3
                            v1 = v1 + icf_he(p,o) * icf_vin(p,o,n)
                        end do
                    end do
                    icfk = icfk - zmat(k,m) * v1
                end if
            end do
        end do
        icfkhist(k,hist_len+hist_fidx_tds) = icfk
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v2

!===============================================================================
! Subroutine:  cst_icf_calculate_v2_mw
! analytical but with numerical/analytical second derivatives
! optimized looping in ICFK, employ Hessian symmetry
! mass-weighted
!===============================================================================

subroutine cst_icf_calculate_v2_mw

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers
    use cst_constraints

    implicit none
    integer                :: k,l,cl,m,cm,n,cn,o,ol,p
    real(PMFDP)            :: v1,v2,icfp,icfk
    ! --------------------------------------------------------------------------

! update CVs - calculate Values, gradients, and Hessians
    call pmf_timers_start_timer(PMFLIB_CST_ICF_HESS_TIMER)
    CVContextP%CVsValues(:) = 0.0d0
    CVContextP%CVsDrvs(:,:,:) = 0.0d0
    CVContextP%CVs2ndDrvs(:,:,:,:,:) = 0.0d0

    do l=1,NumOfAllCONs
        call CONList(l)%cv%calculate_cv2ddrvs(crdhist(:,:,hist_len+hist_fidx_tds),CVContextP)
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_HESS_TIMER)

! get inversion of mass-weighted W
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_constraints_calc_zmat_mw(CVContextP%CVsDrvs)
    call cst_icf_calculate_zmatinv

! ICFP part
    do k=1,NumOfCONs
        call cst_icf_calculate_vi_mw(CVContextP%CVsDrvs,k,icf_vi1)
        icfp = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icfp = icfp + icf_vi1(m,n) * frchist(m,n,hist_len+hist_fidx_tds)
            end do
        end do
        icfphist(k,hist_len+hist_fidx_tds) = - icfp   ! FRC is force not gradient
    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFP_TIMER)

! ICF-K part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFK_TIMER)

    do k=1,NumOfCONs
        icfk = 0.0d0

        ! simpler part: div(M^-1 grad xi_l)
        do l=1,NumOfAllCONs
            cl = CONList(l)%cvindx
            ! get mass-weighted Laplacian
            v1 = 0.0d0
            do ol=1,CONList(l)%cv%natoms
                o = CONList(l)%cv%lindexes(ol)
                do p=1,3
                    v1 = v1 + MassInv(o) * CVContextP%CVs2ndDrvs(p,o,p,o,cl)
                end do
            end do
            icfk = icfk + zmat(k,l) * v1
        end do

        ! harder part: derivative of inverse mass-weighted W
        ! sparse matrix-vector multiplication in cst_icf_calculate_Higj_mw
        icf_vin(:,:,:) = 0.0d0
        do n=1,NumOfAllCONs
            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                do o=1,NumOfLAtoms
                    icf_vin(:,o,n) = icf_vin(:,o,n) + MassInv(o) * zmat(n,l) * CVContextP%CVsDrvs(:,o,cl)
                end do
            end do
        end do
        do n=1,NumOfAllCONs
            cn = CONList(n)%cvindx
            do m=1,n                    ! triangular sum
                cm = CONList(m)%cvindx
                icf_he(:,:) = 0.0d0
                call cst_icf_calculate_Higj_mw(cn,cm)
                call cst_icf_calculate_Higj_mw(cm,cn)
                if( n .ne. m ) then
                    v1 = 0.0d0
                    v2 = 0.0d0
                    do o=1,NumOfLAtoms
                        do p=1,3
                            v1 = v1 + icf_he(p,o) * icf_vin(p,o,n)
                            v2 = v2 + icf_he(p,o) * icf_vin(p,o,m)
                        end do
                    end do
                    icfk = icfk - zmat(k,m) * v1 - zmat(k,n) * v2
                else
                    v1 = 0.0d0
                    do o=1,NumOfLAtoms
                        do p=1,3
                            v1 = v1 + icf_he(p,o) * icf_vin(p,o,n)
                        end do
                    end do
                    icfk = icfk - zmat(k,m) * v1
                end if
            end do
        end do
        icfkhist(k,hist_len+hist_fidx_tds) = icfk
    end do

    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v2_mw


!===============================================================================
! Subroutine:  cst_icf_calculate_v3
! numerical divergence - stochastic “trace trick” for divergence
!===============================================================================

subroutine cst_icf_calculate_v3

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers
    use cst_constraints

    implicit none
    integer                :: k,s,l,cl,n,m
    real(PMFDP)            :: v1,icfp
    ! --------------------------------------------------------------------------

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_constraints_calc_zmat(cvderhist(:,:,:,hist_len+hist_fidx_tds))
    call cst_icf_calculate_zmatll

    do k=1,NumOfCONs
        call cst_icf_calculate_vill(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
        icfp = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icfp = icfp + icf_vi1(m,n) * frchist(m,n,hist_len+hist_fidx_tds)
            end do
        end do
        icfphist(k,hist_len+hist_fidx_tds) = - icfp   ! FRC is force not gradient
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
            call cst_constraints_calc_zmat(CVContextP%CVsDrvs)
            call cst_icf_calculate_zmatll
            call cst_icf_calculate_vill(CVContextP%CVsDrvs,k,icf_vi1)

            icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds) - fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_constraints_calc_zmat(CVContextP%CVsDrvs)
            call cst_icf_calculate_zmatll
            call cst_icf_calculate_vill(CVContextP%CVsDrvs,k,icf_vi2)

            v1 = v1 + sum( sdiv_z(:,:,s)*(icf_vi1(:,:) - icf_vi2(:,:)))

        end do

        if( fpmf_sdiv_qr ) then
            icfkhist(k,hist_len+hist_fidx_tds) =  3.0d0 * real(NumOfLAtoms,PMFDP) * v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        else
            icfkhist(k,hist_len+hist_fidx_tds) =  v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        end if

    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v3

!===============================================================================
! Subroutine:  cst_icf_calculate_v3_mw
! numerical divergence - stochastic “trace trick” for divergence
! mass-weighted
!===============================================================================

subroutine cst_icf_calculate_v3_mw

    use pmf_utils
    use pmf_dat
    use cst_dat
    use pmf_timers
    use cst_constraints

    implicit none
    integer                :: k,s,l,cl,n,m
    real(PMFDP)            :: v1,icfp
    ! --------------------------------------------------------------------------

! ICFP part
    call pmf_timers_start_timer(PMFLIB_CST_ICF_ICFP_TIMER)

    call cst_constraints_calc_zmat_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds))
    call cst_icf_calculate_zmatll

    do k=1,NumOfCONs
        call cst_icf_calculate_vill_mw(cvderhist(:,:,:,hist_len+hist_fidx_tds),k,icf_vi1)
        icfp = 0.0d0
        do n=1,NumOfLAtoms
            do m=1,3
                icfp = icfp + icf_vi1(m,n) * frchist(m,n,hist_len+hist_fidx_tds)
            end do
        end do
        icfphist(k,hist_len+hist_fidx_tds) = - icfp   ! FRC is force not gradient
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
            call cst_constraints_calc_zmat_mw(CVContextP%CVsDrvs)
            call cst_icf_calculate_zmatll
            call cst_icf_calculate_vill_mw(CVContextP%CVsDrvs,k,icf_vi1)

            icf_he(:,:) = crdhist(:,:,hist_len+hist_fidx_tds) - fpmf_sdiv_dh * sdiv_z(:,:,s)

            CVContextP%CVsValues(:) = 0.0d0
            CVContextP%CVsDrvs(:,:,:) = 0.0d0

            do l=1,NumOfAllCONs
                cl = CONList(l)%cvindx
                call CVList(cl)%cv%calculate_cv(icf_he,CVContextP)
            end do
            call cst_constraints_calc_zmat_mw(CVContextP%CVsDrvs)
            call cst_icf_calculate_zmatll
            call cst_icf_calculate_vill_mw(CVContextP%CVsDrvs,k,icf_vi2)

            v1 = v1 + sum( sdiv_z(:,:,s)*(icf_vi1(:,:) - icf_vi2(:,:)) )

        end do

        if( fpmf_sdiv_qr ) then
            icfkhist(k,hist_len+hist_fidx_tds) =  3.0d0 * real(NumOfLAtoms,PMFDP) * v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        else
            icfkhist(k,hist_len+hist_fidx_tds) =  v1 / (2.0d0 * fpmf_sdiv_dh * fpmf_sdiv_S)
        end if

    end do
    call pmf_timers_stop_timer(PMFLIB_CST_ICF_ICFK_TIMER)

end subroutine cst_icf_calculate_v3_mw

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
                    v1 = v1 + CVContextP%CVs2ndDrvs(p,l,o,k,ci) * CVContextP%CVsDrvs(p,l,cj)
                end do
            end do
            icf_he(o,k) = icf_he(o,k) + v1
        end do
    end do

end subroutine cst_icf_calculate_Higj

!===============================================================================
! Subroutine:  cst_icf_calculate_Higj_mw
! Calculate H_i M^-1 g_j contribution for the derivative of mass-weighted W.
!===============================================================================

subroutine cst_icf_calculate_Higj_mw(ci,cj)

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
                    v1 = v1 + CVContextP%CVs2ndDrvs(p,l,o,k,ci) * MassInv(l) * CVContextP%CVsDrvs(p,l,cj)
                end do
            end do
            icf_he(o,k) = icf_he(o,k) + v1
        end do
    end do

end subroutine cst_icf_calculate_Higj_mw


!===============================================================================
! Subroutine:  cst_icf_calculate_zmatll
! LL factorization of ZMAT
!===============================================================================

subroutine cst_icf_calculate_zmatll

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                 :: i,info
    ! --------------------------------------------------------------------------

! zmat diagonal regularization
    do i=1,NumOfAllCONs
        zmat(i,i) = zmat(i,i) + ficf_fdamp
    end do

! calc LL
    call dpotrf('L',NumOfAllCONs,zmat,NumOfAllCONs,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] LL decomposition failed in cst_icf_calculate_zmatll!')
    end if

end subroutine cst_icf_calculate_zmatll

!===============================================================================
! Subroutine:  cst_icf_calculate_vill
! via LL, optimized
!===============================================================================

subroutine cst_icf_calculate_vill(cvsdrvs,i,icf_vi)

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
                         '[CST] LL linear equation failed in cst_icf_calculate_vill!')
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

end subroutine cst_icf_calculate_vill

!===============================================================================
! Subroutine:  cst_icf_calculate_vill_mw
! via LL, optimized, mass-weighted
!===============================================================================

subroutine cst_icf_calculate_vill_mw(cvsdrvs,i,icf_vi)

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
                         '[CST] LL linear equation failed in cst_icf_calculate_vill!')
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

end subroutine cst_icf_calculate_vill_mw

!===============================================================================
! Subroutine:  cst_icf_calculate_zmatinv
!===============================================================================

subroutine cst_icf_calculate_zmatinv

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                 :: i,info
    ! --------------------------------------------------------------------------

! zmat diagonal regularization
    do i=1,NumOfAllCONs
        zmat(i,i) = zmat(i,i) + ficf_fdamp
    end do

! invert
    if ( NumOfAllCONs .gt. 1 ) then
        ! LU decomposition
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
        if( abs(zmat(1,1)) .lt. eps_zmat ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] Singular one-dimensional Z matrix in cst_icf_calculate_zmatinv!')
        end if
        zmat(1,1) = 1.0d0/zmat(1,1)
    end if

end subroutine cst_icf_calculate_zmatinv

!===============================================================================
! Subroutine:  cst_icf_calculate_vi
! zmat is inverted
!===============================================================================

subroutine cst_icf_calculate_vi(cvsdrvs,i,icf_vi)

    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)         :: cvsdrvs(:,:,:)
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
                icf_vi(m,k) = icf_vi(m,k) + zmat(i,j) * cvsdrvs(m,k,cj)
            end do
        end do
    end do

end subroutine cst_icf_calculate_vi

!===============================================================================
! Subroutine:  cst_icf_calculate_vi_mw
! zmat is inverted, mass-weighted
!===============================================================================

subroutine cst_icf_calculate_vi_mw(cvsdrvs,i,icf_vi)

    use pmf_dat
    use cst_dat

    implicit none
    real(PMFDP)         :: cvsdrvs(:,:,:)
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
                icf_vi(m,k) = icf_vi(m,k) + MassInv(k) * zmat(i,j) * cvsdrvs(m,k,cj)
            end do
        end do
    end do

end subroutine cst_icf_calculate_vi_mw


!===============================================================================

end module cst_icf

