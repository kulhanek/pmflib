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
! Subroutine:  abf_core_update_history
! apply ABF force and update history buffers
!===============================================================================

subroutine abf_core_update_history()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,ci
    real(PMFDP)            :: epot,erst,ekin,eus
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        epotrwhist(i)       = epotrwhist(i+1)
        erstrwhist(i)       = erstrwhist(i+1)
        ekinlfhist(i)       = ekinlfhist(i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
    end do

! raw data
    epotrwhist(hist_len)    = PotEne - fepotaverage
    erstrwhist(hist_len)    = PMFEne
    ekinlfhist(hist_len)    = KinEne%KinEneLF - fekinaverage   ! shifted by -1/2dt

! process EPOT
    select case(ftds_epot_src)
        case(1)
            epot = epotrwhist(hist_len)
            erst = erstrwhist(hist_len)
            epothist(hist_len) = epot
            ersthist(hist_len) = erst
        case(2)
            epot = 0.5d0*(epotrwhist(hist_len-0)+epotrwhist(hist_len-2))
            erst = 0.5d0*(erstrwhist(hist_len-0)+erstrwhist(hist_len-2))
            epothist(hist_len-1) = epot
            ersthist(hist_len-1) = erst
        case(3)
            epot = (1.0d0/3.0d0)*(epotrwhist(hist_len-0)+epotrwhist(hist_len-1)+epotrwhist(hist_len-2))
            erst = (1.0d0/3.0d0)*(erstrwhist(hist_len-0)+erstrwhist(hist_len-1)+erstrwhist(hist_len-2))
            epothist(hist_len-1) = epot
            ersthist(hist_len-1) = erst
        case(4)
            epot = (1.0d0/35.0d0)*( -3.0d0*epotrwhist(hist_len-0)+12.0d0*epotrwhist(hist_len-1) &
                                   +17.0d0*epotrwhist(hist_len-2)+12.0d0*epotrwhist(hist_len-3) &
                                    -3.0d0*epotrwhist(hist_len-4))
            erst = (1.0d0/35.0d0)*( -3.0d0*erstrwhist(hist_len-0)+12.0d0*erstrwhist(hist_len-1) &
                                   +17.0d0*erstrwhist(hist_len-2)+12.0d0*erstrwhist(hist_len-3) &
                                    -3.0d0*erstrwhist(hist_len-4))
            epothist(hist_len-2) = epot
            ersthist(hist_len-2) = erst
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_epot_src in abf_core_update_history!')
    end select

! process EKIN
    select case(ftds_ekin_src)
        ! KE from velocities at full-step interpolated from velocities at half-step
        case(1)
            ekinhist(hist_len-1)    = KinEne%KinEneVV - fekinaverage   ! shifted by -dt
            ekin                    = ekinhist(hist_len-1)
        ! KE from interpolated KE at half-step
        case(4)
            ekinhist(hist_len-1)    = 0.5d0*(ekinlfhist(hist_len-0) + ekinlfhist(hist_len-1))
            ekin                    = ekinhist(hist_len-1)
        case(5)
            ekinhist(hist_len-2)    = (1.0d0/16.0d0)*(      -ekinlfhist(hist_len-0)+9.0d0*ekinlfhist(hist_len-1) &
                                                      +9.0d0*ekinlfhist(hist_len-2)      -ekinlfhist(hist_len-3))
            ekin                    = ekinhist(hist_len-2)
        case(6)
            ekinhist(hist_len-3)    = (1.0d0/256.0d0)*(  +3.0d0*ekinlfhist(hist_len-0)  -25.0d0*ekinlfhist(hist_len-1) &
                                                       +150.0d0*ekinlfhist(hist_len-2) +150.0d0*ekinlfhist(hist_len-3) &
                                                        -25.0d0*ekinlfhist(hist_len-4)   +3.0d0*ekinlfhist(hist_len-5))
            ekin                    = ekinhist(hist_len-3)
        case(9)
            ekinhist(hist_len-2)    = KinEne%KinEneVV - fekinaverage
            ekin                    = ekinhist(hist_len-2)
        case(10)
            ekinhist(hist_len-0)    = KinEne%KinEneVV - fekinaverage
            ekin                    = ekinhist(hist_len-0)
        case(11)
            ekinhist(hist_len-1)    = KinEne%KinEneHA - fekinaverage
            ekin                    = ekinhist(hist_len-1)
    case default
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented ftds_ekin_src mode in abf_core_update_history!')
    end select

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

! apply force filters
    la(:) = 0.0d0
    if( fapply_abf ) then
        if( fusmode ) then
            ! US-ABF
            call abf_core_get_us_bias(cvhist(:,hist_len),la,eus)
            ! since ersthist(hist_len) is already updated
            ! we can add the bias into energy PMFEne
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
        ! FIXME
        do i=1,abfaccu%tot_cvs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if
    micfhist(:,hist_len) = la(:)

    if( frecord ) then
        ! record time progress of data
        call abf_accu_add_data_record_lf(cvhist(:,hist_len),micfhist(:,hist_len), &
                                         epot,erst,ekin)
    end if

    return

end subroutine abf_core_update_history

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
! Subroutine:  abf_core_register_rawdata
!===============================================================================

subroutine abf_core_register_rawdata(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin)

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    real(PMFDP),intent(in)  :: cvs(:)
    real(PMFDP),intent(in)  :: ficf(:)
    real(PMFDP),intent(in)  :: sicf(:)
    real(PMFDP),intent(in)  :: vicf(:)
    real(PMFDP),intent(in)  :: bicf(:)
    real(PMFDP),intent(in)  :: epot
    real(PMFDP),intent(in)  :: erst
    real(PMFDP),intent(in)  :: ekin
    ! --------------------------------------------
    real(PMFDP)             :: etot, ekin_scaled
    ! --------------------------------------------------------------------------

    ! total ABF force
    pxia(:) = ficf(:) + sicf(:) + vicf(:) ! adaptive correction
    ! bicf  ! current bias

    ! scale ekin
    ekin_scaled = ekin * ftds_ekin_scale

    ! total energy
    etot = epot + erst + ekin_scaled

    ! add data to accumulator
    ! subroutine abf_accu_add_data_online(cvs,gfx,epot,erst,ekin,etot)
    call abf_accu_add_data_online(cvs,pxia,bicf,epot,erst,ekin_scaled,etot)

    if( fentropy .and. fentdecomp ) then
        ! subroutine abf_accu_add_data_entropy_decompose(cvs,fx,sx,vx,lx,bx,epot,erst,ekin)
        call abf_accu_add_data_entropy_decompose(cvs,ficf,sicf,vicf,bicf,epot,erst,ekin_scaled)
    end if

end subroutine abf_core_register_rawdata

!===============================================================================
! subroutine:  abf_core_calc_Zmat
!===============================================================================

subroutine abf_core_calc_Zmat(ctx)

    use pmf_utils
    use abf_dat

    implicit none
    type(CVContextType) :: ctx
    integer             :: i,ci,j,cj,k,info
    real(PMFDP)         :: v,logdet
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

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then
        fzinv(:,:)  = fz(:,:)
        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzinv,NumOfABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if

        call dgetri(NumOfABFCVs,fzinv,NumOfABFCVs,indx,vv,NumOfABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_core_calc_Zmat!')
        end if
    else
        fzinv(1,1)  = 1.0d0/fz(1,1)
    end if

    return

end subroutine abf_core_calc_Zmat

!===============================================================================

end module abf_core
