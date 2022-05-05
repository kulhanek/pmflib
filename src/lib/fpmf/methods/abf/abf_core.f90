!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

integer,parameter ::  DEBUG_ABF_FMODE1    = 4750
integer,parameter ::  DEBUG_ABF_FMODE2    = 4751
integer,parameter ::  DEBUG_ABF_FMODE3    = 4752
integer,parameter ::  DEBUG_ABF_FMODE4    = 4753
integer,parameter ::  DEBUG_ABF_GPR       = 4789

contains

!===============================================================================
! Subroutine:  abf_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine abf_core_main

    use abf_trajectory
    use abf_restart
    use abf_output
    use abf_client
    use abf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    select case(fmode)
        ! standard algorithms
        case(1)
            !call abf_core_force_2p
            call abf_core_force_3pB
        case(2)
            call abf_core_force_gpr
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    call abf_client_exchange_data(.false.)

end subroutine abf_core_main

!===============================================================================
! Subroutine:  abf_core_force_3pB
! this is leap-frog ABF version, simplified algorithm
! forces from velocities
!===============================================================================

subroutine abf_core_force_3pB()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v1,v2,f1,epot,erst,ekin
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)

! shift epot ene
    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v1 = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v1 = v1 + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zdhist(m,j,i,hist_len) = v1
            end do
        end do
    end do

! apply force filters
    la(:) = 0.0d0
    if( fapply_abf ) then
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
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if
    micfhist(:,hist_len) = la(:)

    if( frecord ) then
        ! record time progress of data
        call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
                                         epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
     end if

! ABF part
    if( fstep .ge. hist_len ) then
        do i=1,NumOfABFCVs
            f1  = 0.0d0
            v1 = 0.0d0
            v2 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f1 = f1 + zdhist(m,j,i,hist_len-1) * ( &
                         (vhist(m,j,hist_len-0)-vhist(m,j,hist_len-1) / LNG_c_explic) / LNG_c_implic )
                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                end do
            end do
            pxi0(i) = f1*ifdtx
            pxi1(i) = 0.5d0*(v1+v2)*ifdtx
        end do

        ! total ABF force
        pxip(:) = pxi0(:) + pxi1(:) - micfhist(:,hist_len-1)  ! unbiased estimate

        epot = epothist(hist_len-1)
        erst = ersthist(hist_len-1)
        ekin = ekinhist(hist_len-1)

        ! debug
        ! write(1225,*) epot,erst,ekin

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxip,epot,erst,ekin)

        if( fentropy .and. fentdecomp ) then
            call abf_accu_add_data_entropy_decompose(cvhist(:,hist_len-1),pxi0,pxi1,micfhist(:,hist_len-1),epot,erst,ekin)
        end if
    end if

    return

end subroutine abf_core_force_3pB

!===============================================================================
! Subroutine:  abf_core_force_2p
! this is leap-frog ABF version, simplified algorithm
!===============================================================================

subroutine abf_core_force_2p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use abf_output

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,e
    ! --------------------------------------------------------------------------

    ! shift accuvalue history
    cvaluehist0(:) = cvaluehist1(:)

    ! save coordinate value to history
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvaluehist1(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    if( fenthalpy ) then
        epothist1 = PotEne - fepotaverage
    else
        epothist1 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    if( fentropy ) then
        ekinhist1 = KinEne - fekinaverage
    else
        ekinhist1 = 0.0d0
    end if

    ! shift erst ene
    ersthist0 = ersthist1
    if( fentropy ) then
        ersthist1 = PMFEne
    else
        ersthist1 = 0.0d0
    end if

    ! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)

    do i=1,NumOfABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfABFCVs
        v = 0.0d0
        e = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! zd0 in t-dt
                ! Vel in t-1/2dt
                ! v0 (OldVel) in t-3/2dt
                ! a <- Vel(m,j)-v0(m,j))/fdtx in t-dt, do not use forces (Frc) because of SHAKE
                v = v + zd0(m,j,i)*(Vel(m,j)-v0(m,j))
                ! zd1 in t
                ! zd0 in t-dt
                ! vel in t-1/2dt
                e = e + (zd1(m,j,i)-zd0(m,j,i))* Vel(m,j)
            end do
        end do
        pxi0(i) = v / fdtx ! in t-dt
        pxip(i) = e / fdtx ! in t-1/2dt
    end do

    write(1225,*) fstep-1,v,e

    if( fstep .ge. 4 ) then

        ! complete ICF in t-dt
        ! pxi0 in t-dt
        ! pxi1 - old ABF forces in t-dt
        ! pxip in t-1/2dt
        ! pxim in t-3/2dt
        pxi0(:) = pxi0(:) - pxi1(:)
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))

        ! total ABF force
        pxi0(:) = pxi0(:) + pxim(:)

        ! add data to accumulator
        call abf_accu_add_data_online(cvaluehist0,pxi0(:),epothist0,ekinhist1,ersthist0)
    end if

    ! backup to the next step
    zd0  = zd1
    pxim = pxip
    v0   = Vel

    ! apply ABF bias
    la(:) = 0.0d0

    ! apply force filters
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvaluehist1(:),la)
            case(1)
                call abf_accu_get_data_lramp(cvaluehist1(:),la)
            case(2)
                ! call abf_accu_get_data_gks(cvaluehist1(:),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if

    ! keep ABF forces to subtract them in the next step
    pxi1 = la

    return

end subroutine abf_core_force_2p

!===============================================================================
! Subroutine:  abf_core_force_gpr
! this is leap-frog ABF version, GPR algorithm
! forces GPR process
!===============================================================================

subroutine abf_core_force_gpr()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers
    use abf_init

    implicit none
    integer     :: i,j,k,l1,l2,p1,p2
    integer     :: ci,gpr_mid
    real(PMFDP) :: v,epot,erst,ekin,mean
    real(PMFDP) :: invn,depot1,depot2
    ! --------------------------------------------------------------------------

    ! in this algorithm, we use circular buffer
    cbuff_pos = cbuff_pos + 1
    if( cbuff_pos .gt. gpr_len ) then
        cbuff_pos = 1
    end if

! update values
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,cbuff_pos)  = CVContext%CVsValues(ci)
    end do

    epothist(cbuff_pos)     = PotEne - fepotaverage
    ersthist(cbuff_pos)     = PMFEne
    if( cbuff_pos-1 .le. 0 ) then
        ekinhist(gpr_len)    = KinEne - fekinaverage    ! shifted by -dt
    else
        ekinhist(cbuff_pos-1) = KinEne - fekinaverage    ! shifted by -dt
    end if

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)
    fzinvhist(:,:,cbuff_pos) = fzinv(:,:)

! apply ABF force
    la(:) = 0.0d0
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call abf_accu_get_data(cvhist(:,cbuff_pos),la)
            case(1)
                call abf_accu_get_data_lramp(cvhist(:,cbuff_pos),la)
            case(2)
                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                    call abf_accu_get_data_ksmooth(cvhist(:,cbuff_pos),la)
                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
            case(3)
                call abf_accu_get_data_lsmooth(cvhist(:,cbuff_pos),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfABFCVs
            ci = ABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if
    micfhist(:,cbuff_pos)    = la(:)

    if( frecord ) then
        ! record time progress of data
        if( cbuff_pos-1 .le. 0 ) then
            ekin = ekinhist(gpr_len)
        else
            ekin = ekinhist(cbuff_pos-1)
        end if
        call abf_accu_add_data_record_lf(cvhist(:,cbuff_pos),fzinvhist(:,:,cbuff_pos),micfhist(:,cbuff_pos), &
                                         epothist(cbuff_pos),ersthist(cbuff_pos),ekin)
     end if

    if( fstep .le. gpr_len )  return

    ! predict cv velocity at gpr_len/2+1 centered around cbuff_pos - gpr_len/2
    gpr_mid = cbuff_pos - gpr_len/2
    if( gpr_mid .le. 0 ) then
        gpr_mid = gpr_mid + gpr_len
    end if

    ! put CV velocities in pxi1
    do i=1,NumOfABFCVs
        ! calculate mean value
        mean = 0.0d0

        if( .not. gpr_cvs_nomean ) then
            do k=1,gpr_len
                mean = mean + cvhist(i,k)
            end do
            mean = mean / real(gpr_len,PMFDP)
        end if

        ! shift data
        l1 = cbuff_pos
        do k=gpr_len,1,-1
            gpr_data(k) = cvhist(i,l1) - mean
            l1 = l1 - 1
            if( l1 .le. 0 ) l1 = l1 + gpr_len
        end do

        pxi1(i) = dot_product(gpr_data,gpr_kfd_cvs)

        if( gpr_calc_logxx ) then
            ! calculate logML
            ! solve GPR
            call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_cvs,gpr_len,gpr_data,1,0.0d0,gpr_model,1)

                ! calculate CV derivative in time - derivative is shift invariant
            gpr_logml(i) = -0.5d0*dot_product(gpr_data,gpr_model)   &
                        - 0.5d0*gpr_K_cvs_logdet                    &
                        - 0.5d0*real(gpr_len,PMFDP)*log(2.0*PMF_PI)

            gpr_logpl(i) = - 0.5d0*real(gpr_len,PMFDP)*log(2.0*PMF_PI)
            do k=1,gpr_len
                gpr_logpl(i) = gpr_logpl(i) + 0.5d0*log(gpr_K_cvs(k,k)) - 0.5d0*gpr_model(k)**2/gpr_K_cvs(k,k)
            end do

            if ( fdebug ) then
                write(789,*) fstep, gpr_logml(1), gpr_logpl(1)
            end if

            ! increase number of samples
            gpr_nlogxx = gpr_nlogxx + 1.0d0
            invn = 1.0d0 / gpr_nlogxx

            ! statistics
            depot1 = gpr_logml(i) - gpr_mlogml(i)
            gpr_mlogml(i)  = gpr_mlogml(i)  + depot1 * invn
            depot2 = gpr_logml(i) - gpr_mlogml(i)
            gpr_m2logml(i) = gpr_m2logml(i) + depot1 * depot2

            ! statistics
            depot1 = gpr_logpl(i) - gpr_mlogpl(i)
            gpr_mlogpl(i)  = gpr_mlogpl(i)  + depot1 * invn
            depot2 = gpr_logpl(i) - gpr_mlogpl(i)
            gpr_m2logpl(i) = gpr_m2logpl(i) + depot1 * depot2
        end if
    end do

! calculate momenta
    do i=1,NumOfABFCVs
        v = 0.0d0
        do j=1,NumOfABFCVs
            v = v + fzinvhist(i,j,gpr_mid) * pxi1(j)
        end do
        xphist(i,gpr_mid) = v
    end do

    ! write(148569,*) fstep, gpr_mid, cvhist(1,gpr_mid), pxi1(1), xphist(1,gpr_mid)

! calculate derivatives of CV momenta
    select case(gpr_icf_cdf)
        case(0)
            ! move gpr_mid by one element left so we have three valid points
            gpr_mid = gpr_mid - 1
            if( gpr_mid .le. 0 ) gpr_mid = gpr_mid + gpr_len

            ! central differences - 2p+1 - indexes
            l1 = gpr_mid - 1
            if( l1 .le. 0 ) l1 = l1 + gpr_len
            p1 = gpr_mid + 1
            if( p1 .gt. gpr_len ) p1 = p1 - gpr_len

            ! central differences - 2p+1
            do i=1,NumOfABFCVs
                icfhist(i,gpr_mid) = 0.5d0*(xphist(i,p1)-xphist(i,l1))*ifdtx
            end do
        case(1)
            ! move gpr_mid by two elements left so we have five valid points
            gpr_mid = gpr_mid - 2
            if( gpr_mid .le. 0 ) gpr_mid = gpr_mid + gpr_len

            ! central differences - 4p+1 - indexes
            l2 = gpr_mid - 2
            if( l2 .le. 0 ) l2 = l2 + gpr_len
            l1 = gpr_mid - 1
            if( l1 .le. 0 ) l1 = l1 + gpr_len
            p1 = gpr_mid + 1
            if( p1 .gt. gpr_len ) p1 = p1 - gpr_len
            p2 = gpr_mid + 2
            if( p2 .gt. gpr_len ) p2 = p2 - gpr_len

            ! central differences - 4p+1
            do i=1,NumOfABFCVs
                icfhist(i,gpr_mid) = (1.0d0/12.0d0)*(-xphist(i,p2) +8.0d0*xphist(i,p1) -8.0d0*xphist(i,l1) +xphist(i,l2))*ifdtx
            end do
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented gpr_icf_cdf mode!')
    end select

! record new data
    pxi0(:) = icfhist(:,gpr_mid) - micfhist(:,gpr_mid)

    epot = epothist(gpr_mid)
    erst = ersthist(gpr_mid)
    ekin = ekinhist(gpr_mid)

    if( fdebug ) then
        write(DEBUG_ABF_FMODE4,*) fstep-hist_len+k, cvhist(:,k), pxi0, epot, erst, ekin
    end if

    ! add data to accumulator
    call abf_accu_add_data_online(cvhist(:,gpr_mid),pxi0,epot,erst,ekin)

end subroutine abf_core_force_gpr

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

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFCVs .gt. 1) then

        fzinv(:,:)  = fz(:,:)

        call dgetrf(NumOfABFCVs,NumOfABFCVs,fzinv,NumOfABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if

        fzdet = 1.0d0
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

end module abf_core
