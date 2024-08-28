!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more detajls.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module abf_core

use pmf_sizes
use pmf_constants

implicit none
contains


!===============================================================================
! Subroutine:  abf_core_get_us_bias
! get bias from US restraints
!===============================================================================

subroutine abf_core_get_us_bias(values,gfx,bene)

    use abf_dat

    implicit none
    real(PMFDP)     :: values(:)
    real(PMFDP)     :: gfx(:)
    real(PMFDP)     :: bene
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    bene = 0.0

    do i=1,NumOfABFCVs

        ABFCVList(i)%deviation = ABFCVList(i)%cv%get_deviation(values(i),ABFCVList(i)%target_value)

        ABFCVList(i)%energy = 0.5d0*ABFCVList(i)%force_constant*ABFCVList(i)%deviation**2
        bene = bene + ABFCVList(i)%energy

        gfx(i) = - ABFCVList(i)%force_constant*ABFCVList(i)%deviation
    end do

end subroutine abf_core_get_us_bias

!===============================================================================
! Subroutine:  abf_core_update_history_force
! apply ABF force and update history buffers
!===============================================================================

subroutine abf_core_update_history_force()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,ci
    real(PMFDP)            :: eus
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        xhist(:,:,i)        = xhist(:,:,i+1)
        micfhist(:,i)       = micfhist(:,i+1)
        icfhist(:,i)        = icfhist(:,i+1)
        fdetzhist(i)        = fdetzhist(i+1)
        fziihist(:,i)       = fziihist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
    end do

    xhist(:,:,hist_len)     = Crd(:,:)

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)
    fdetzhist(hist_len) = fzdet
    do i=1,NumOfABFCVs
        fziihist(i,hist_len) = fz(i,i)
    end do

! apply force filters
    la(:) = 0.0d0
    if( fapply_abf ) then
        if( fusmode ) then
            ! US-ABF
            call abf_core_get_us_bias(cvhist(:,hist_len),la,eus)
            PMFEne = PMFEne + eus
        else
            ! regular ABF
            ! calculate abf force to be applied
            select case(feimode)
                case(0)
                    call abf_accu_get_data(cvhist(:,hist_len),la)
                case(1)
                    call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
                case(2)
                    call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                        call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
                    call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
                case(3)
                    call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
                case default
                    call pmf_utils_exit(PMF_OUT,1, &
                         '[ABF] Not implemented extrapolation/interpolation mode in abf_core_update_history!')
            end select
        end if

        ! project abf force along coordinate
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if
    micfhist(:,hist_len) = la(:)

    return

end subroutine abf_core_update_history_force

!===============================================================================
! Subroutine:  abf_core_update_cvder
!===============================================================================

subroutine abf_core_update_cvder()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,m,ki
    ! --------------------------------------------------------------------------

    do i=1,NumOfABFCVs
        ki = ABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            do m=1,3
                cvderhist(m,j,i,hist_len) = CVContext%CVsDrvs(m,j,ki)
            end do
        end do
    end do

end subroutine abf_core_update_cvder

!===============================================================================
! subroutine:  abf_core_calc_Zmat
!===============================================================================

subroutine abf_core_calc_Zmat(ctx)

    use pmf_utils
    use abf_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        do j=1,NumOfABFCVs
            cj = ABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(ctx%CVsDrvs(:,k,ci),ctx%CVsDrvs(:,k,cj))
            end do
        end do
    end do

    fzdet = 1.0d0

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then
        fzinv(:,:)  = fz(:,:)
        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzinv,NumOfABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if

         ! and finally determinant
        do i=1,NumOfABFCVs
            if( indx(i) .ne. i ) then
                fzdet = - fzdet * fzinv(i,i)
            else
                fzdet = fzdet * fzinv(i,i)
            end if
        end do

        call dgetri(NumOfABFCVs,fzinv,NumOfABFCVs,indx,vv,NumOfABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_core_calc_Zmat!')
        end if
    else
        fzdet       = fz(1,1)
        fzinv(1,1)  = 1.0d0/fz(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat

!===============================================================================
! Subroutine:  abf_core_update_zdhist
!===============================================================================

subroutine abf_core_update_zdhist()

    use pmf_dat
    use pmf_cvs
    use abf_dat

    implicit none
    integer                :: i,j,k,m,ki
    real(PMFDP)            :: v
    ! --------------------------------------------------------------------------

    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zdhist(m,j,i,hist_len) = v
            end do
        end do
    end do

end subroutine abf_core_update_zdhist

!===============================================================================

end module abf_core
