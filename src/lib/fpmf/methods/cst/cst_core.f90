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
    use cst_lambda
    use cst_icf
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

    epothist(hist_len)          = PotEne - fepotaverage
    ersthist(hist_len)          = PMFEne

    crdhist(:,:,hist_len)       = Crd(:,:)
    cvderhist(:,:,:,hist_len)   = CVContext%CVsDrvs(:,:,:)
    frchist(:,:,hist_len)       = Frc(:,:)
    velhist(:,:,hist_len)       = Vel(:,:)

    if( fintcalc .and. fint_der ) then
        call cst_icf_calculate_icf
    end if

    if( ftdscalc ) then
        call cst_lambda_calculate
    end if

    select case(fintalg)
        case(IA_LEAP_FROG)
            lambdaMhist(:,hist_len)      = lambda(:)
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
    use cst_icf

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
                lambdaMhist(:,hist_len) = lambda(:)
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
    integer                :: i,ci,j,cj,k,m,info
    real(PMFDP)            :: z1,v1,fzdet,fzdets
    ! --------------------------------------------------------------------------

! ALL constraints ================================

! calculate Z matrix at Crd (in t)
    call cst_constraints_calc_zmat_mw(CVContext%CVsDRvs)

! calculate Z determinant ------------------------------------
    if( NumOfAllCONs .gt. 1 ) then
        ! LU decomposition
        call dgetrf(NumOfAllCONs,NumOfAllCONs,zmat,NumOfAllCONs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_core_calculate_fw!')
        end if
        fzdet = 1.0d0
        ! and finally determinant
        do i=1,NumOfAllCONs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * zmat(i,i)
            else
                fzdet = fzdet * zmat(i,i)
            end if
        end do
    else
        fzdet = zmat(1,1)
    end if

! SHAKE constraints ==============================

    fzdets = 1.0d0

    if( frmmdcon_zdet ) then
    ! calculate Z matrix at Crd (in t)
        do i=1,NumOfMDCONs
            ci = CONList(i+NumOfCONs)%cvindx
            do j=1,i
                cj = CONList(j+NumOfCONs)%cvindx
                z1 = 0.0d0
                do k=1,NumOfLAtoms
                    v1 = 0.0
                    do m=1,3
                        v1 = v1 + CVContext%CVsDrvs(m,k,ci)*CVContext%CVsDrvs(m,k,cj)
                    end do
                    z1 = z1 + MassInv(k)*v1
                end do
                zmats(i,j)=z1
                zmats(j,i)=z1
            end do
        end do

    ! calculate Z determinant ------------------------------------
        if( NumOfMDCONs .gt. 1 ) then
            ! LU decomposition
            call dgetrf(NumOfMDCONs,NumOfMDCONs,zmats,NumOfMDCONs,indx,info)
            if( info .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_core_calculate_fw!')
            end if
            fzdets = 1.0d0
            ! and finally determinant
            do i=1,NumOfMDCONs
                if( indx(i) .ne. i ) then
                    fzdets = - fzdets * zmats(i,i)
                else
                    fzdets = fzdets * zmats(i,i)
                end if
            end do
        else if( NumOfMDCONs .eq. 1 ) then
            fzdets = zmats(1,1)
        end if
    end if

    fzdet   = fzdet / fzdets
    fwfac   = 1.0d0/sqrt(fzdet)
    fwhist(hist_len) = fwfac

end subroutine cst_core_calculate_fw

!===============================================================================
! Subroutine:  cst_core_shift_histbuffs
!===============================================================================

subroutine cst_core_shift_histbuffs

    use cst_dat

    implicit none
    integer :: i
    ! --------------------------------------------------------------------------

    do i=1,hist_len-1
        lambdaMhist(:,i)    = lambdaMhist(:,i+1)
        lambdaEhist(:,i)    = lambdaEhist(:,i+1)
        fwhist(i)           = fwhist(i+1)

        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        enevalidhist(i)     = enevalidhist(i+1)

        icfphist(:,i)       = icfphist(:,i+1)
        icfkhist(:,i)       = icfkhist(:,i+1)

        crdhist(:,:,i)      = crdhist(:,:,i+1)
        cvderhist(:,:,:,i)  = cvderhist(:,:,:,i+1)
        frchist(:,:,i)      = frchist(:,:,i+1)
        velhist(:,:,i)      = velhist(:,:,i+1)
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

