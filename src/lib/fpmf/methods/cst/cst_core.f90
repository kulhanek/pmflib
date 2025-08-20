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

    call cst_core_calculate_fw
    call cst_shake_calculate

    lambda(:) = lambda(:) + lambdax(:) * isfdts

    if( fintene .and. fintene_der ) then
        call cst_core_calculate_icf
        icfphist(:,hist_len) = icfp(:)
        icfkhist(:,hist_len) = icfk(:)
    end if

    select case(fintalg)
        case(IA_LEAP_FROG)
            lambdahist(:,hist_len) = lambda(:)
            epothist(hist_len) = PotEne - fepotaverage
            ersthist(hist_len) = PMFEne
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
!                if( fenthalpy_der ) then
!                    call cst_core_calculate_icfp
!                    icfphist(:,hist_len) = icfp(:)
!                end if
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
! Subroutine:  cst_core_calculate_fw
! it uses CVContext
!===============================================================================

subroutine cst_core_calculate_fw

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                :: i,ci,j,cj,k,info
    real(PMFDP)            :: jacv,fzdeta,fzdets,fzdet
    ! --------------------------------------------------------------------------

! ALL constraints ================================

! calculate Z matrix at Crd (in t)
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        do j=1,NumOfAllCONs
            cj = CONList(j)%cvindx
            jacv = 0.0
            do k=1,NumOfLAtoms
                jacv = jacv + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            zmata(i,j) = jacv
        end do
    end do

! calculate Z determinant ------------------------------------
    if( NumOfAllCONs .gt. 1 ) then
        ! LU decomposition
        call dgetrf(NumOfAllCONs,NumOfAllCONs,zmata,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_core_calculate_fw!')
        end if
        fzdeta = 1.0d0
        ! and finally determinant
        do i=1,NumOfAllCONs
            if( indx(i) .ne. i ) then
                fzdeta = - fzdeta * zmata(i,i)
            else
                fzdeta = fzdeta * zmata(i,i)
            end if
        end do
    else
        fzdeta = zmata(1,1)
    end if

! SHAKE constraints ==============================

    fzdets = 1.0d0

    if( frmshake_zdet ) then
    ! calculate Z matrix at Crd (in t)
        do i=1,NumOfSHAKECONs
            ci = CONList(i+NumOfCONs)%cvindx
            do j=1,NumOfSHAKECONs
                cj = CONList(j+NumOfCONs)%cvindx
                jacv = 0.0
                do k=1,NumOfLAtoms
                    jacv = jacv + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
                end do
                zmats(i,j) = jacv
            end do
        end do

    ! calculate Z determinant ------------------------------------
        if( NumOfSHAKECONs .gt. 1 ) then
            ! LU decomposition
            call dgetrf(NumOfSHAKECONs,NumOfSHAKECONs,zmats,NumOfSHAKECONs,indx,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_core_calculate_fw!')
            end if
            fzdets = 1.0d0
            ! and finally determinant
            do i=1,NumOfSHAKECONs
                if( indx(i) .ne. i ) then
                    fzdets = - fzdets * zmats(i,i)
                else
                    fzdets = fzdets * zmats(i,i)
                end if
            end do
        else if( NumOfSHAKECONs .eq. 1 ) then
            fzdets = zmats(1,1)
        else
            fzdets = 1.0d0
        end if
    end if

! record data
! DOI: 10.1080/00268970310001592746 - eq. 6
    fzdet   = fzdeta / fzdets
    fwfac   = 1.0d0/sqrt(fzdet)
    fwhist(hist_len) = fwfac

end subroutine cst_core_calculate_fw


!===============================================================================
! Subroutine:  cst_core_calculate_icf
!===============================================================================

subroutine cst_core_calculate_icf

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

end subroutine cst_core_calculate_icf

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
        fwhist(i)       = fwhist(i+1)
        icfphist(:,i)   = icfphist(:,i+1)
        icfkhist(:,i)   = icfkhist(:,i+1)
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
    use cst_accu

    implicit none
    integer         :: i,ci
    ! --------------------------------------------------------------------------

! reset accumulators
    if ( faccurst .eq. 0 ) then
        faccurst = -1

        call cst_accu_clear

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
    do i=1,NumOfAllCONs
        ci = CONList(i)%cvindx
        CONList(i)%deviation = CONList(i)%cv%get_deviation(CVContextP%CVsValues(ci),CONList(i)%value)   ! t+dt
        CONList(i)%sdevtot = CONList(i)%sdevtot + CONList(i)%deviation**2                               ! t+dt
    end do

! do we have enough samples?
    if( fstep .le. 2 ) return

! record data
    call cst_accu_add_lam
    call cst_accu_add_dhTds

end subroutine cst_core_analyze

!===============================================================================

end module cst_core

