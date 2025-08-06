!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module cst_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_core_main_lf
!===============================================================================

subroutine cst_core_main_lf

    use cst_constraints
    use cst_shake
    use cst_output
    use cst_restart
    use cst_trajectory
    use pmf_utils

    implicit none
    ! --------------------------------------------------------------------------

    select case(fintalg)
        case(IA_LEAP_FROG)
            call cst_constraints_increment
            call cst_core_shift_histbuffs
            lambda(:) = 0.0d0
        case(IA_LF_MIDDLE)
            ! nothing to be here
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_core_main_lf!')
    end select

    call cst_core_calculate_zdet
    call cst_shake_calculate

    lambda(:) = lambda(:) + lambdax(:) * isfdts

    select case(fintalg)
        case(IA_LEAP_FROG)
            lambdahist(:,hist_len) = lambda(:)
            epothist(hist_len) = PotEne - fepotaverage
            ersthist(hist_len) = PMFEne
            if( fenthalpy_der ) then
                call cst_core_calculate_icfp
                icfphist(:,hist_len) = icfp(:)
            end if
            call cst_core_analyze
            call cst_output_write
            call cst_restart_update
            call cst_trajectory_write_snapshot
        case(IA_LF_MIDDLE)
            ! nothing to be here
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_core_main_lf!')
    end select

end subroutine cst_core_main_lf

!===============================================================================
! Subroutine:  cst_core_rattlev_lf
!===============================================================================

subroutine cst_core_rattlev_lf(cid)

    use cst_constraints
    use cst_rattlev
    use pmf_dat
    use pmf_utils
    use cst_output
    use cst_restart
    use cst_trajectory

    implicit none
    integer :: cid      ! call id from MD engine
                        ! in the LF-middle, there are two rattle-v calls
    ! --------------------------------------------------------------------------

    select case(fintalg)
        case(IA_LEAP_FROG)
            ! nothing to be here
        case(IA_LF_MIDDLE)
            if( cid .eq. 1 ) then
                call cst_constraints_increment
                call cst_core_shift_histbuffs
                lambda(:) = 0.0d0
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_core_main_lf!')
    end select

    call cst_rattlev_calculate

    lambda(:) = lambda(:) + lambdav(:) * isfdtr

    select case(fintalg)
        case(IA_LEAP_FROG)
            ! nothing to be here
        case(IA_LF_MIDDLE)
            if( cid .eq. 2 ) then
                lambdahist(:,hist_len) = lambda(:)
                epothist(hist_len) = PotEne - fepotaverage
                ersthist(hist_len) = PMFEne
                if( fenthalpy_der ) then
                    call cst_core_calculate_icfp
                    icfphist(:,hist_len) = icfp(:)
                end if
                call cst_core_analyze
                call cst_output_write
                call cst_restart_update
                call cst_trajectory_write_snapshot
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported integration algorithm in cst_core_main_lf!')
    end select

end subroutine cst_core_rattlev_lf

!===============================================================================
! Subroutine:  cst_core_register_ekin
!===============================================================================

subroutine cst_core_register_ekin_lf

    use pmf_dat
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    enevalidhist(hist_len)  = KinEne%Valid
    ekinhist(hist_len)      = KinEne%KinEneVV - fekinaverage

end subroutine cst_core_register_ekin_lf

!===============================================================================
! Subroutine:  cst_core_calculate_zdet
!===============================================================================

subroutine cst_core_calculate_zdet

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,j,cj,k,info
    real(PMFDP)            :: jacv,isrz
    ! --------------------------------------------------------------------------

! calculate Z matrix at Crd (in t)
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfCONs
            cj = CONList(j)%cvindx
            jacv = 0.0
            do k=1,NumOfLAtoms
                jacv = jacv + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            jac(i,j) = jacv
        end do
    end do

! calculate Z determinant ------------------------------------
    if( NumOfCONs .gt. 1 ) then
        ! LU decomposition
        call dgetrf(NumOfCONs,NumOfCONs,jac,NumOfCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_core_calculate_zdet!')
        end if
        fzdet = 1.0d0
        ! and finally determinant
        do i=1,NumOfCONs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * jac(i,i)
            else
                fzdet = fzdet * jac(i,i)
            end if
        end do
    else
        fzdet = jac(1,1)
    end if

! record data
    isrz    = 1.0d0/sqrt(fzdet)
    isrzhist(hist_len) = isrz

end subroutine cst_core_calculate_zdet

!===============================================================================
! Subroutine:  cst_core_calculate_icfp
!===============================================================================

subroutine cst_core_calculate_icfp

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,k,m
    real(PMFDP)            :: f1,nv
    ! --------------------------------------------------------------------------

    ! start with dV/dx
    CSTFrc(:,:) = Frc(:,:)

    ! add constraint forces
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        do k=1,NumOfLAtoms
            CSTFrc(:,k) = CSTFrc(:,k) + lambda(i)*CVContext%CVsDrvs(:,k,ci)
        end do
    end do

    ! project to CVs
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        f1 = 0.0d0
        nv = 0.0d0
        do k=1,NumOfLAtoms
            do m=1,3
                ! force part
                nv = nv + CVContext%CVsDrvs(m,k,ci) * CVContext%CVsDrvs(m,k,ci)
                f1 = f1 + CVContext%CVsDrvs(m,k,ci) * CSTFrc(m,k)
            end do
        end do
        icfp(i) = f1 / nv
    end do

end subroutine cst_core_calculate_icfp

!===============================================================================
! Subroutine:  cst_core_shift_histbuffs
!===============================================================================

subroutine cst_core_shift_histbuffs

    use cst_dat

    implicit none
    integer :: i
    ! --------------------------------------------------------------------------

    do i=1,hist_len-1
        lambdahist(:,i) = lambdahist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        isrzhist(i)     = isrzhist(i+1)
        enevalidhist(i) = enevalidhist(i+1)
    end do

end subroutine cst_core_shift_histbuffs

!===============================================================================
! Subroutine:  cst_core_analyze
!===============================================================================

subroutine cst_core_analyze

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer         :: i,ci
    ! --------------------------------------------------------------------------

! reset accumulators
    if ( faccurst .eq. 0 ) then
        faccurst = -1

        ! free energy calculation
        nsamples    = 0.0d0
        mlambda(:)  = 0.0d0
        m2lambda(:) = 0.0d0
        misrz       = 0.0d0
        m2isrz      = 0.0d0

        ! accumulator setup for entropy and enthalpy
        fene_step = 0

        if( fenthalpy .or. fentropy ) then
            ntds = 0.0d0
        end if

        if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
            meint       = 0.0d0
            m2eint      = 0.0d0
            mepot       = 0.0d0
            m2epot      = 0.0d0
            merst       = 0.0d0
            m2erst      = 0.0d0
            mekin       = 0.0d0
            m2ekin      = 0.0d0
        end if

        if( fenthalpy .and. fenthalpy_der ) then
            micfp(:)    = 0.0d0
            m2icfp(:)   = 0.0d0
            c11pp(:)    = 0.0d0
        end if

        if( fentropy ) then
            metot       = 0.0d0
            m2etot      = 0.0d0
            mpp(:)      = 0.0d0
            m2pp(:)     = 0.0d0
            mpn(:)      = 0.0d0
            m2pn(:)     = 0.0d0
            mhicf(:)    = 0.0d0
            m2hicf(:)   = 0.0d0
        end if

        if( fentropy .and. fentdecomp ) then
            c11hp(:)    = 0.0d0
            c11hr(:)    = 0.0d0
            c11hk(:)    = 0.0d0
        end if

        CONList(:)%sdevtot = 0.0d0

        write(CST_OUT,'(A)') '#-------------------------------------------------------------------------------'
        write(CST_OUT,'(A)') '# INFO: ALL ACCUMULATORS WERE RESETED                                           '
        write(CST_OUT,'(A)') '#       PRODUCTION STAGE OF ACCUMULATION IS STARTED                             '
        write(CST_OUT,'(A)') '#-------------------------------------------------------------------------------'
    end if

    ! accumulate results -----------------------------------------------------
    if( faccurst .gt. 0 ) then
        faccurst = faccurst - 1
    end if

! calculate final constraint deviations
    do i=1,NumOfCONs
        ci = CONList(i)%cvindx
        CONList(i)%deviation = CONList(i)%cv%get_deviation(CVContextP%CVsValues(ci),CONList(i)%value)   ! t+dt
        CONList(i)%sdevtot = CONList(i)%sdevtot + CONList(i)%deviation**2                               ! t+dt
    end do

! do we have enough samples?
    if( fstep .le. 2 ) return

! record data
    call cst_core_analyze_lam
    call cst_core_analyze_dhTds

end subroutine cst_core_analyze

!===============================================================================
! Subroutine:  cst_core_analyze_lam
! free energy
!===============================================================================

subroutine cst_core_analyze_lam

    use pmf_dat
    use cst_dat

    implicit none
    integer         :: i
    real(PMFDP)     :: isrz,lam,dval1,dval2,invn
    ! --------------------------------------------------------------------------

    if( mod(fstep,flamsample) .ne. 0 ) return

! values
    lambda(:)   = lambdahist(:,hist_len+hist_fidx)
    isrz        = isrzhist(hist_len+hist_fidx)     ! t-dt

    nsamples = nsamples + 1
    if( nsamples .le. 0 ) return
    invn = 1.0d0/nsamples

    do i=1,NumOfCONs
        ! lambda
        lam             = lambda(i)
        dval1           = lam - mlambda(i)
        mlambda(i)      = mlambda(i)  + dval1 * invn
        dval2           = lam - mlambda(i)
        m2lambda(i)     = m2lambda(i) + dval1 * dval2
    end do

! isrz
    dval1   = isrz - misrz
    misrz   = misrz  + dval1 * invn
    dval2   = isrz - misrz
    m2isrz  = m2isrz + dval1*dval2

end subroutine cst_core_analyze_lam

!===============================================================================
! Subroutine:  cst_core_analyze_dhTds
! enthalpy and entropy
!===============================================================================

subroutine cst_core_analyze_dhTds

    use pmf_dat
    use cst_dat

    implicit none
    integer         :: i
    real(PMFDP)     :: lam,dval1,dval2,invn
    real(PMFDP)     :: etot,epot,erst,eint,ekin
    real(PMFDP)     :: detot1, detot2
    real(PMFDP)     :: depot1, depot2
    real(PMFDP)     :: derst1, derst2
    real(PMFDP)     :: deint1, deint2
    real(PMFDP)     :: dekin1, dekin2
    real(PMFDP)     :: dpp, dpp1, dpp2
    real(PMFDP)     :: dpn, dpn1, dpn2
    real(PMFDP)     :: dicf1, dicf2
    ! --------------------------------------------------------------------------

    if( enevalidhist(hist_len+hist_fidx) ) fene_step = fene_step + 1
    if( .not. ( (mod(fene_step,fenesample) .eq. 0) .and. enevalidhist(hist_len+hist_fidx) ) ) return

    ntds = ntds + 1.0d0
    invn = 1.0d0/ntds

    lambda(:)   = lambdahist(:,hist_len+hist_fidx)
    epot        = epothist(hist_len+hist_fidx)     ! t-dt
    erst        = ersthist(hist_len+hist_fidx)     ! t-dt
    ekin        = ekinhist(hist_len+hist_fidx)     ! t-dt
    etot        = epot + erst + ekin               ! t-dt
    eint        = epot + erst

    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
        ! internal energy
        deint1 = eint - meint
        meint  = meint + deint1 * invn
        deint2 = eint - meint
        m2eint = m2eint + deint1 * deint2

        ! potential energy
        depot1 = epot - mepot
        mepot  = mepot + depot1 * invn
        depot2 = epot - mepot
        m2epot = m2epot + depot1 * depot2

        ! restraint energy
        derst1 = erst - merst
        merst  = merst + derst1 * invn
        derst2 = erst - merst
        m2erst = m2erst + derst1 * derst2

        ! kinetic energy
        dekin1 = ekin - mekin
        mekin  = mekin + dekin1 * invn
        dekin2 = ekin - mekin
        m2ekin = m2ekin + dekin1 * dekin2
    end if

    if( fenthalpy .and. fenthalpy_der ) then
        do i=1,NumOfCONs
            dicf1     = icfp(i) - micfp(i)
            micfp(i)  = micfp(i) + dicf1 * invn
            dicf2     = icfp(i) - micfp(i)
            m2icfp(i) = m2icfp(i) + dicf1 * dicf2

            c11pp(i)  = c11pp(i) + dicf1 * depot2
        end do
    end if

    if( fentropy ) then
        ! total energy
        detot1 = etot - metot
        metot  = metot + detot1 * invn
        detot2 = etot - metot
        m2etot = m2etot + detot1 * detot2
    end if

! lambda and entropy
    if( fentropy ) then
        do i=1,NumOfCONs

            ! lambda
            lam             = lambda(i)
            dval1           = lam - mhicf(i)
            mhicf(i)        = mhicf(i) + dval1 * invn
            dval2           = lam - mhicf(i)
            m2hicf(i)       = m2hicf(i) + dval1 * dval2

            dpp     = lam + etot
            dpp1    = dpp - mpp(i)
            mpp(i)  = mpp(i) + dpp1 * invn
            dpp2    = dpp - mpp(i)
            m2pp(i) = m2pp(i) + dpp1 * dpp2

            dpn     = lam - etot
            dpn1    = dpn - mpn(i)
            mpn(i)  = mpn(i) + dpn1 * invn
            dpn2    = dpn - mpn(i)
            m2pn(i) = m2pn(i) + dpn1 * dpn2

            if( fentdecomp ) then
                c11hp(i)   = c11hp(i) + dval1 * depot2
                c11hr(i)   = c11hr(i) + dval1 * derst2
                c11hk(i)   = c11hk(i) + dval1 * dekin2
            end if
        end do
    end if

end subroutine cst_core_analyze_dhTds

!===============================================================================

end module cst_core

