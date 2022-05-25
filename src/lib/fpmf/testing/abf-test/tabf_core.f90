!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module tabf_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  tabf_core_main
! this is leap-frog and velocity-verlet ABF version
!===============================================================================

subroutine tabf_core_main

    use tabf_trajectory
    use tabf_restart
    use tabf_output
    use tabf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    select case(fmode)
        case(1)
            ! numerical differentiation
            call tabf_core_force_2p
        case(2)
            ! standard ABF algorithm
            call tabf_core_force_4p
        case(3)
            ! analytical/numerical differentiation
            call tabf_core_force_2p_frc  ! SHAKE must be off
        case default
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Not implemented fmode in tabf_core_main!')
    end select

    call tabf_output_write
    call tabf_trajectory_write_snapshot
    call tabf_restart_update

end subroutine tabf_core_main

!!===============================================================================
!! Subroutine:  abf_core_shake
!! correct for forces from SHAKE
!!===============================================================================
!
!subroutine abf_core_shake
!
!    use abf_dat
!    use pmf_utils
!
!    ! --------------------------------------------------------------------------
!
!    select case(fmode)
!        case(10)
!            ! fix total forces from SHAKE
!            fhist(:,:,hist_len) =  fhist(:,:,hist_len) + SHAKEFrc(:,:)
!        case(1,2,11)
!            ! nothing to be done
!        case default
!            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main_shake!')
!    end select
!
!end subroutine abf_core_shake

!===============================================================================
! Subroutine:  tabf_core_force_4p
! this is leap-frog ABF version, original implementation
!===============================================================================

subroutine tabf_core_force_4p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use tabf_dat
    use tabf_accu
    use tabf_output

    implicit none
    integer     :: i,j,k,m
    integer     :: ci,ck
    real(PMFDP) :: v,avg_epot,avg_ekin,avg_erst
    ! --------------------------------------------------------------------------

    ! calculate acceleration in time t for all pmf atoms
    do i=1,NumOfLAtoms
        a1(:,i) = MassInv(i)*Frc(:,i)
    end do

    ! shift accuvalue history
    cvaluehist0(:) = cvaluehist1(:)
    cvaluehist1(:) = cvaluehist2(:)
    cvaluehist2(:) = cvaluehist3(:)

    ! save coordinate value to history
    do i=1,NumOfTABFCVs
        ci = TABFCVList(i)%cvindx
        cvaluehist3(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    epothist1 = epothist2
    epothist2 = epothist3
    if( fenthalpy .or. fentropy ) then
        epothist3 = PotEne + fepotoffset
    else
        epothist3 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    ekinhist1 = ekinhist2
    ekinhist2 = ekinhist3
    if( fentropy ) then
        ekinhist3 = KinEneVV + fekinoffset
    else
        ekinhist3 = 0.0d0
    end if

    ! shift erst ene
    ersthist0 = ersthist1
    ersthist1 = ersthist2
    ersthist2 = ersthist3
    if( fentropy ) then
        ersthist3 = PMFEne
    else
        ersthist3 = 0.0d0
    end if

    ! calculate abf force to be applied -------------
    select case(feimode)
        case(0)
            call tabf_accu_get_data(cvaluehist3(:),la)
        case(1)
            call tabf_accu_get_data_lramp(cvaluehist3(:),la)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Not implemented extrapolation/interpolation mode!')
    end select

    ! apply force filters
    if( .not. fapply_abf ) then
        ! we do not want to apply ABF force - set the forces to zero
        la(:) = 0.0d0
    end if

    ! project abf force along coordinate ------------
    do i=1,NumOfTABFCVs
        ci = TABFCVList(i)%cvindx
        do j=1,NumOfLAtoms
            a1(:,j) = a1(:,j) + la(i) * MassInv(j) * CVContext%CVsDrvs(:,j,ci)
        end do
    end do

    ! rest of ABF stuff -----------------------------

    ! calculate Z matrix and its inverse
    call tabf_core_calc_Zmat

    ! pxip = zd0(t-dt)*[v(t-dt/2)/2 - dt*a1(t)/12]
    do i=1,NumOfTABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) - fdtx*a1(m,j)/12.0)
            end do
        end do
        pxip(i) = v
    end do

    ! ZD0(3, NumOfLAtoms, NumOfTABFCVs)   <-- 1/dt * m_ksi grad ksi(r0)
    do i=1,NumOfTABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfTABFCVs
                    ck = TABFCVList(k)%cvindx
                    v = v + fzinv(i,k) * CVContext%CVsDrvs(m,j,ck)
                end do
                zd0(m,j,i) = v / fdtx
            end do
        end do
    end do

    ! pxim = zd0(t)*[v(t-dt/2)/2 + dt*a0(t-dt)/12]
    do i=1,NumOfTABFCVs
        v = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                v = v + zd0(m,j,i) * (0.5d0*Vel(m,j) + fdtx*a0(m,j)/12.0)
            end do
        end do
        pxim(i) = v
    end do

    !a0 <-- a1
    a0(:,:) = a1(:,:)

    ! pxi0 <--- pxi0 + pxip
    pxi0(:) = pxi0(:) + pxip(:)

    ! update accumulator - we need at least four samples
    if( fstep .ge. 4 ) then
        ! calculate coordinate values at time t-3/2dt
        do i=1,NumOfTABFCVs
            avg_values(i) = TABFCVList(i)%cv%get_average_value(cvaluehist1(i),cvaluehist2(i))
        end do

        avg_epot = 0.5d0*(epothist1 + epothist2) ! t - 3/2*dt
        avg_ekin = 0.5d0*(ekinhist2 + ekinhist3) ! t - 1/2*dt; ekin already shifted by -dt
        avg_erst = 0.5d0*(ersthist1 + ersthist2) ! t - 3/2*dt

        pdum(:) = 0.0d0

        ! add data to accumulator
        if( fblock_size .gt. 0 ) then
            call tabf_accu_add_data_blocking(cvaluehist0,pxi0(:),pdum(:),avg_epot,avg_ekin,avg_erst)
        else
            call tabf_accu_add_data_online(cvaluehist0,pxi0(:),pdum(:),avg_epot,avg_ekin,avg_erst)
        end if

        ! write icf
        call tabf_output_write_icf(avg_values,pxi0(:))
    end if

    ! pxi0 <--- -pxip + pxim + pxi1 - la/2
    pxi0(:) = -pxip(:) + pxim(:) + pxi1(:) - 0.5d0*la(:)

    ! pxi1 <--- -pxim - la/2
    pxi1(:) = -pxim(:) - 0.5d0*la(:)

    ! project updated acceleration back
    do i=1,NumOfLAtoms
        Frc(:,i) = a1(:,i)*Mass(i)
    end do

    return

end subroutine tabf_core_force_4p

!===============================================================================
! Subroutine:  tabf_core_force_2p
! this is leap-frog ABF version
!===============================================================================

subroutine tabf_core_force_2p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use tabf_dat
    use tabf_accu
    use tabf_output

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,e
    ! --------------------------------------------------------------------------

    ! shift accuvalue history
    cvaluehist0(:) = cvaluehist1(:)

    ! save coordinate value to history
    do i=1,NumOfTABFCVs
        ci = TABFCVList(i)%cvindx
        cvaluehist1(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    if( fenthalpy .or. fentropy ) then
        epothist1 = PotEne + fepotoffset
    else
        epothist1 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    if( fentropy ) then
        ekinhist1 = KinEneVV + fekinoffset
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
    call tabf_core_calc_Zmat

    do i=1,NumOfTABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfTABFCVs
                    ki = TABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfTABFCVs
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

    if( fstep .ge. 4 ) then

        ! complete ICF in t-dt
        ! pxi0 in t-dt
        ! pxi1 - old ABF forces in t-dt
        ! pxip in t-1/2dt
        ! pxim in t-3/2dt
        pxi0(:) = pxi0(:) - pxi1(:)
        pxim(:) = 0.5d0*(pxim(:)+pxip(:))

        ! add data to accumulator
        if( fblock_size .gt. 0 ) then
            call tabf_accu_add_data_blocking(cvaluehist0,pxi0(:),pxim(:),epothist0,ekinhist1,ersthist0)
        else
            call tabf_accu_add_data_online(cvaluehist0,pxi0(:),pxim(:),epothist0,ekinhist1,ersthist0)
        end if
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
                call tabf_accu_get_data(cvaluehist1(:),la)
            case(1)
                call tabf_accu_get_data_lramp(cvaluehist1(:),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[TABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfTABFCVs
            ci = TABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if

    ! keep ABF forces to subtract them in the next step
    pxi1 = la

    return

end subroutine tabf_core_force_2p

!===============================================================================
! Subroutine:  tabf_core_force_2p_frc
! this is leap-frog ABF version
!===============================================================================

subroutine tabf_core_force_2p_frc()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use tabf_dat
    use tabf_accu
    use tabf_output

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,e
    ! --------------------------------------------------------------------------

    ! shift accuvalue history
    cvaluehist0(:) = cvaluehist1(:)

    ! save coordinate value to history
    do i=1,NumOfTABFCVs
        ci = TABFCVList(i)%cvindx
        cvaluehist1(i) = CVContext%CVsValues(ci)
    end do

    ! shift epot ene
    epothist0 = epothist1
    if( fenthalpy .or. fentropy ) then
        epothist1 = PotEne + fepotoffset
    else
        epothist1 = 0.0d0
    end if

    ! shift ekin ene
    ekinhist0 = ekinhist1
    if( fentropy ) then
        ekinhist1 = KinEneVV + fekinoffset
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
    call tabf_core_calc_Zmat

    do i=1,NumOfTABFCVs
        do j=1,NumOfLAtoms
            do m=1,3
                v = 0.0d0
                do k=1,NumOfTABFCVs
                    ki = TABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zd1(m,j,i) = v
            end do
        end do
    end do

    do i=1,NumOfTABFCVs
        v = 0.0d0
        e = 0.0d0
        do j=1,NumOfLAtoms
            do m=1,3
                ! zd1 in t
                ! Frc in t
                v = v + zd1(m,j,i)*MassInv(j)*Frc(m,j)
                ! zd1 in t
                ! zd0 in t-dt
                ! vel in t-1/2dt
                e = e + (zd1(m,j,i)-zd0(m,j,i))* Vel(m,j)
            end do
        end do
        pxi1(i) = v        ! in t
        pxip(i) = e / fdtx ! in t-1/2dt
    end do

    if( fstep .ge. 4 ) then

        ! complete ICF in t-dt
        ! pxi0 in t-dt         - potene contribution
        ! pxip in t-1/2dt      - kinetic contributions
        ! pxim in t-3/2dt
        pxim(:) =  0.5d0*(pxim(:)+pxip(:))

        ! add data to accumulator
        if( fblock_size .gt. 0 ) then
            call tabf_accu_add_data_blocking(cvaluehist0,pxi0(:),pxim(:),epothist0,ekinhist1,ersthist0)
        else
            call tabf_accu_add_data_online(cvaluehist0,pxi0(:),pxim(:),epothist0,ekinhist1,ersthist0)
        end if
    end if

    ! backup to the next step
    zd0  = zd1
    pxi0 = pxi1
    pxim = pxip

    ! apply ABF bias
    la(:) = 0.0d0

    ! apply force filters
    if( fapply_abf ) then
        ! calculate abf force to be applied
        select case(feimode)
            case(0)
                call tabf_accu_get_data(cvaluehist1(:),la)
            case(1)
                call tabf_accu_get_data_lramp(cvaluehist1(:),la)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[TABF] Not implemented extrapolation/interpolation mode!')
        end select

        ! project abf force along coordinate
        do i=1,NumOfTABFCVs
            ci = TABFCVList(i)%cvindx
            do j=1,NumOfLAtoms
                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
            end do
        end do
    end if

    return

end subroutine tabf_core_force_2p_frc

!!===============================================================================
!! Subroutine:  abf_core_force_3pA
!! this is leap-frog ABF version, simplified algorithm
!! employing forces from potential and SHAKE
!!===============================================================================
!
!subroutine abf_core_force_3pA()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!
!    implicit none
!    integer                :: i,j,k,m,ifac
!    integer                :: ci,ki
!    real(PMFDP)            :: v,v1,v2,f,etot,epot,erst,ekin
!    ! --------------------------------------------------------------------------
!
!    if( fdebug .and. (mod(fstep,5000) .eq. 0) ) then
!        open(unit=4789,file='abf-fmode1.bicf',status='UNKNOWN')
!        ifac = 10
!        do i=1,ABFCVList(1)%nbins*ifac
!            cvave(1) = ABFCVList(1)%min_value &
!                     + real(i,PMFDP) * (ABFCVList(1)%max_value-ABFCVList(1)%min_value)/(ABFCVList(1)%nbins * ifac)
!            la(:) = 0.0d0
!            pxi0(:) = 0.0d0
!            call abf_accu_get_data(cvave,pxi0)
!            select case(feimode)
!                case(0)
!                    call abf_accu_get_data(cvave,la)
!                case(1)
!                    call abf_accu_get_data_lramp(cvave,la)
!                case(2)
!                    call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                        call abf_accu_get_data_ksmooth(cvave,la)
!                    call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!                case(3)
!                    call abf_accu_get_data_lsmooth(cvave,la)
!                case default
!                    call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
!            end select
!            write(4789,*) cvave, pxi0, la, sfac
!        end do
!        close(4789)
!    end if
!
!! shift values
!    do i=1,hist_len-1
!        cvhist(:,i)     = cvhist(:,i+1)
!        epothist(i)     = epothist(i+1)
!        ersthist(i)     = ersthist(i+1)
!        ekinhist(i)     = ekinhist(i+1)
!        vhist(:,:,i)    = vhist(:,:,i+1)
!        fhist(:,:,i)    = fhist(:,:,i+1)
!        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
!        micfhist(:,i)   = micfhist(:,i+1)
!    end do
!
!! update values
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len) = CVContext%CVsValues(ci)
!    end do
!    vhist(:,:,hist_len)     = Vel(:,:)   ! in t-dt/2
!
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
!
!! apply ABF force
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len)    = la(:)
!    fhist(:,:,hist_len)     = Frc(:,:)   ! in t, SHAKE forces are added later
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!
!    do i=1,NumOfABFCVs
!        do j=1,NumOfLAtoms
!            do m=1,3
!                v = 0.0d0
!                do k=1,NumOfABFCVs
!                    ki = ABFCVList(k)%cvindx
!                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
!                end do
!                zdhist(m,j,i,hist_len) = v
!            end do
!        end do
!    end do
!
!    ! record time progress of data
!    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
!                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!
!! ABF part
!    if( fstep .ge. hist_len ) then
!
!        do i=1,NumOfABFCVs
!            f  = 0.0d0
!            v1 = 0.0d0
!            v2 = 0.0d0
!            do j=1,NumOfLAtoms
!                do m=1,3
!                    ! force part
!                    f = f + zdhist(m,j,i,hist_len-1)*fhist(m,j,hist_len-1)  * MassInv(j)
!                    ! velocity part
!                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
!                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
!                end do
!            end do
!            pxi0(i) = f + 0.5d0*(v1+v2)*ifdtx
!        end do
!
!        ! total ABF force
!        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)
!
!        epot = epothist(hist_len-1)
!        erst = ersthist(hist_len-1)
!        ekin = ekinhist(hist_len-1)
!        etot = epot + erst + ekin
!
!        if( fdebug ) then
!            write(DEBUG_ABF_FMODE1,*) fstep-1, cvhist(:,hist_len-1), pxi0, epot, erst, ekin, etot
!        end if
!
!        ! add data to accumulator
!        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,ekin,etot)
!    end if
!
!    return
!
!end subroutine abf_core_force_3pA
!
!!===============================================================================
!! Subroutine:  abf_core_force_3pC
!! this is leap-frog ABF version, simplified algorithm
!! forces from positions
!!===============================================================================
!
!subroutine abf_core_force_3pC()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!
!    implicit none
!    integer                :: i,j,k,m
!    integer                :: ci,ki
!    real(PMFDP)            :: v,f,etot,epot,erst,ekin
!    ! --------------------------------------------------------------------------
!
!! shift accuvalue history
!    do i=1,hist_len-1
!        cvhist(:,i)     = cvhist(:,i+1)
!        epothist(i)     = epothist(i+1)
!        ersthist(i)     = ersthist(i+1)
!        ekinhist(i)     = ekinhist(i+1)
!        xhist(:,:,i)    = xhist(:,:,i+1)
!        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
!        micfhist(:,i)   = micfhist(:,i+1)
!    end do
!
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len) = CVContext%CVsValues(ci)
!    end do
!    xhist(:,:,hist_len)    = Crd(:,:)
!
!! shift epot ene
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!
!    do i=1,NumOfABFCVs
!        do j=1,NumOfLAtoms
!            do m=1,3
!                v = 0.0d0
!                do k=1,NumOfABFCVs
!                    ki = ABFCVList(k)%cvindx
!                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
!                end do
!                zdhist(m,j,i,hist_len) = v
!            end do
!        end do
!    end do
!
!! apply force filters
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len) = la(:)
!
!    ! record time progress of data
!    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
!                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!
!
!! ABF part
!    if( fstep .ge. hist_len ) then
!        do i=1,NumOfABFCVs
!            f  = 0.0d0
!            v  = 0.0d0
!            do j=1,NumOfLAtoms
!                do m=1,3
!                    ! force part
!                    f = f + zdhist(m,j,i,hist_len-1) &
!                      * (xhist(m,j,hist_len-0)-2.0d0*xhist(m,j,hist_len-1)+xhist(m,j,hist_len-2))
!                    ! velocity part
!                    v = v + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-2)) &
!                       * (xhist(m,j,hist_len-0)-xhist(m,j,hist_len-2))
!                end do
!            end do
!            pxi0(i) = (f + 0.25d0*v) * ifdtx * ifdtx
!        end do
!
!        ! total ABF force
!        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)  ! unbiased estimate
!
!        epot = epothist(hist_len-1)
!        erst = ersthist(hist_len-1)
!        ekin = ekinhist(hist_len-1)
!        etot = epot + erst + ekin
!
!        ! add data to accumulator
!        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,ekin,etot)
!    end if
!
!    return
!
!end subroutine abf_core_force_3pC
!
!!===============================================================================
!! Subroutine:  abf_core_force_gpr
!! this is leap-frog ABF version, simplified algorithm
!! forces GPR process
!!===============================================================================
!
!subroutine abf_core_force_gpr()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!    use abf_init
!
!    implicit none
!    integer                     :: i,j,k
!    integer                     :: ci,skip
!    real(PMFDP)                 :: v,etot,epot,erst,ekin,mean
!    character(len=PMF_MAX_TYPE) :: cvid
!    ! --------------------------------------------------------------------------
!
!! shift values
!    do i=1,hist_len-1
!        cvhist(:,i)         = cvhist(:,i+1)
!        epothist(i)         = epothist(i+1)
!        ersthist(i)         = ersthist(i+1)
!        ekinhist(i)         = ekinhist(i+1)
!        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
!        micfhist(:,i)       = micfhist(:,i+1)
!    end do
!
!! update values
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
!    end do
!
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!    fzinvhist(:,:,hist_len) = fzinv(:,:)
!
!! apply ABF force
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len)    = la(:)
!
!    ! record time progress of data
!    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
!                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!
!
!    ! use hist_len not gpr_len here, hist_len = gpr_len + 1 due to ekin delay
!    if( mod(fstep,hist_len) .ne. 0 )  return
!
!! calculate CV velocities
!    do i=1,NumOfABFCVs
!        ! calculate mean value
!        mean = 0.0d0
!        do k=1,gpr_len
!            mean = mean + cvhist(i,k)
!        end do
!        mean = mean / real(gpr_len,PMFDP)
!
!        ! shift data
!        do k=1,gpr_len
!            gpr_data(k) = cvhist(i,k) - mean
!        end do
!
!        ! solve GPR
!        call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_cvs,gpr_len,gpr_data,1,0.0d0,gpr_model,1)
!
!        do k=1,gpr_len
!            ! calculate CV derivative in time - derivative is shift invariant
!            xvelhist(i,k) = dot_product(gpr_model,gpr_kfd_cvs(:,k))
!        end do
!
!        if( fdebug ) then
!            write(cvid,'(I0.3)') i
!            open(unit=4789,file='abf-gpr.cvhist_'//trim(cvid),status='UNKNOWN')
!            do k=1,gpr_len
!                write(4789,*) k, dot_product(gpr_model,gpr_kff_cvs(:,k))+mean, cvhist(i,k)
!            end do
!            close(4789)
!            open(unit=4789,file='abf-gpr.xvelhist_'//trim(cvid),status='UNKNOWN')
!            write(4789,*) 1, xvelhist(i,1)
!            do k=2,gpr_len-1
!                write(4789,*) k, xvelhist(i,k), 0.5d0*(cvhist(i,k+1)-cvhist(i,k-1))*ifdtx
!            end do
!            write(4789,*) gpr_len, xvelhist(i,gpr_len)
!            close(4789)
!        end if
!    end do
!
!! calculate momenta
!    do k=1,gpr_len
!        do i=1,NumOfABFCVs
!            v = 0.0d0
!            do j=1,NumOfABFCVs
!                v = v + fzinvhist(i,j,k) * xvelhist(j,k)
!            end do
!            xphist(i,k) = v
!        end do
!    end do
!
!! calculate derivatives of CV momenta
!    if( gpr_icf_cdf ) then
!        do i=1,NumOfABFCVs
!            do k=2,gpr_len-1
!                icfhist(i,k) = 0.5d0*(xphist(i,k+1)-xphist(i,k-1))*ifdtx
!            end do
!        end do
!    else
!        do i=1,NumOfABFCVs
!
!            ! setup data
!            do k=1,gpr_len
!                gpr_data(k) = xphist(i,k)
!            end do
!
!            ! solve GPR
!            call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_icf,gpr_len,gpr_data,1,0.0d0,gpr_model,1)
!
!            do k=1,gpr_len
!                ! calculate momenta derivative in time - derivative is shift invariant
!                icfhist(i,k) = dot_product(gpr_model,gpr_kfd_icf(:,k))
!
!            end do
!
!            if( fdebug ) then
!                write(cvid,'(I0.3)') i
!                open(unit=4789,file='abf-gpr.xphist_'//trim(cvid),status='UNKNOWN')
!                do k=1,gpr_len
!                    write(4789,*) k, dot_product(gpr_model,gpr_kff_icf(:,k)), xphist(i,k)
!                end do
!                close(4789)
!                open(unit=4789,file='abf-gpr.icfhist_'//trim(cvid),status='UNKNOWN')
!                write(4789,*) 1, icfhist(i,1)
!                do k=2,gpr_len-1
!                    write(4789,*) k, icfhist(i,k), 0.5d0*(xphist(i,k+1)-xphist(i,k-1))*ifdtx
!                end do
!                write(4789,*) gpr_len, icfhist(i,gpr_len)
!                close(4789)
!            end if
!        end do
!    end if
!! record data
!    if( gpr_ene_smooth .gt. 0 ) then
!        select case(gpr_ene_smooth)
!            case(1)
!                mean = 0.0d0
!                do k=1,gpr_len
!                    mean = mean + epothist(k) + ersthist(k) + ekinhist(k)
!                end do
!                mean = mean / real(gpr_len,PMFDP)
!
!                ! shift data
!                do k=1,gpr_len
!                    gpr_data(k) = epothist(k) + ersthist(k) + ekinhist(k) - mean
!                end do
!            case(2)
!                mean = 0.0d0
!                do k=1,gpr_len
!                    mean = mean + ekinhist(k)
!                end do
!                mean = mean / real(gpr_len,PMFDP)
!
!                ! shift data
!                do k=1,gpr_len
!                    gpr_data(k) = ekinhist(k) - mean
!                end do
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported gpr_ene_smooth in abf_core_force_gpr!')
!        end select
!
!        ! solve GPR
!        call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_ene,gpr_len,gpr_data,1,0.0d0,gpr_model,1)
!
!        select case(gpr_ene_smooth)
!            case(1)
!                open(unit=4789,file='abf-gpr.etot',status='UNKNOWN')
!                do k=1,gpr_len
!                    write(4789,*) k, dot_product(gpr_model,gpr_kff_ene(:,k))+mean, gpr_data(k)+mean
!                end do
!                close(4789)
!            case(2)
!                open(unit=4789,file='abf-gpr.ekin',status='UNKNOWN')
!                do k=1,gpr_len
!                    write(4789,*) k, dot_product(gpr_model,gpr_kff_ene(:,k))+mean, gpr_data(k)+mean
!                end do
!                close(4789)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported gpr_ene_smooth in abf_core_force_gpr!')
!        end select
!    end if
!
!    skip = 0
!    if( gpr_icf_cdf ) then
!        skip = 1
!    end if
!
!    do k=1+gpr_boundary+skip,gpr_len-gpr_boundary-skip
!
!        ! total ABF force
!        pxi0(:) = icfhist(:,k) - micfhist(:,k)
!
!        epot = epothist(k)
!        erst = ersthist(k)
!
!        select case(gpr_ene_smooth)
!            case(0)
!                ekin = ekinhist(k)
!                etot = epot + erst + ekin
!            case(1)
!                etot = dot_product(gpr_model,gpr_kff_ene(:,k)) + mean
!            case(2)
!                ekin = dot_product(gpr_model,gpr_kff_ene(:,k)) + mean
!                etot = epot + erst + ekin
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported gpr_ene_smooth in abf_core_force_gpr!')
!        end select
!
!        if( fdebug ) then
!            write(DEBUG_ABF_FMODE4,*) fstep-hist_len+k, cvhist(:,k), pxi0, etot
!        end if
!
!        ! add data to accumulator
!        call abf_accu_add_data_online(cvhist(:,k),pxi0,epot,erst,ekin,etot)
!    end do
!
!end subroutine abf_core_force_gpr
!
!!===============================================================================
!! Subroutine:  abf_core_force_sg
!! this is leap-frog ABF version
!! SG filter for differentiation
!!===============================================================================
!
!subroutine abf_core_force_sg()
!
!    use pmf_utils
!    use pmf_dat
!    use pmf_cvs
!    use abf_dat
!    use abf_accu
!    use pmf_timers
!
!    implicit none
!    integer                :: i, j, ci
!    real(PMFDP)            :: vp, vm, epot, erst, ekin, etot
!    ! --------------------------------------------------------------------------
!
!! shift values
!    do i=1,hist_len-1
!        cvhist(:,i)         = cvhist(:,i+1)
!        epothist(i)         = epothist(i+1)
!        ersthist(i)         = ersthist(i+1)
!        ekinhist(i)         = ekinhist(i+1)
!        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
!        micfhist(:,i)       = micfhist(:,i+1)
!    end do
!
!! update values
!    do i=1,NumOfABFCVs
!        ci = ABFCVList(i)%cvindx
!        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
!    end do
!
!    epothist(hist_len)      = PotEne - fepotaverage
!    ersthist(hist_len)      = PMFEne
!    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt
!
!! calculate Z matrix and its inverse
!    call abf_core_calc_Zmat(CVContext)
!    fzinvhist(:,:,hist_len) = fzinv(:,:)
!
!! apply ABF force
!    la(:) = 0.0d0
!    if( fapply_abf ) then
!        ! calculate abf force to be applied
!        select case(feimode)
!            case(0)
!                call abf_accu_get_data(cvhist(:,hist_len),la)
!            case(1)
!                call abf_accu_get_data_lramp(cvhist(:,hist_len),la)
!            case(2)
!                call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
!                    call abf_accu_get_data_ksmooth(cvhist(:,hist_len),la)
!                call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
!            case(3)
!                call abf_accu_get_data_lsmooth(cvhist(:,hist_len),la)
!            case default
!                call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
!        end select
!
!        ! project abf force along coordinate
!        do i=1,NumOfABFCVs
!            ci = ABFCVList(i)%cvindx
!            do j=1,NumOfLAtoms
!                Frc(:,j) = Frc(:,j) + la(i) * CVContext%CVsDrvs(:,j,ci)
!            end do
!        end do
!    end if
!    micfhist(:,hist_len)    = la(:)
!
!    ! record time progress of data
!    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
!                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))
!
!    if( fstep .lt. hist_len ) return
!
!    do i=1,NumOfABFCVs
!        pxi0(i) = dot_product(sg_c1(:),cvhist(i,1:hist_len-2))
!        pxi1(i) = dot_product(sg_c1(:),cvhist(i,3:hist_len-0))
!    end do
!
!! calculate momenta
!    do i=1,NumOfABFCVs
!        vm = 0.0d0
!        vp = 0.0d0
!        do j=1,NumOfABFCVs
!            vm = vm + fzinvhist(i,j,hist_len/2+0) * pxi0(j)
!            vp = vp + fzinvhist(i,j,hist_len/2+2) * pxi1(j)
!        end do
!        pxip(i) = vp
!        pxim(i) = vm
!    end do
!
!! calculate derivatives of CV momenta
!    do i=1,NumOfABFCVs
!        pxi0(i) = 0.5d0*(pxip(i)-pxim(i))*ifdtx
!    end do
!
!    epot = epothist(hist_len/2+1)   ! hist_len is +2 bigger than fsgframelen
!    erst = ersthist(hist_len/2+1)
!    ekin = ekinhist(hist_len/2+1)
!    etot = epot + erst + ekin
!
!    pxi0(:)     = pxi0(:) - micfhist(:,hist_len/2+1)
!
!    ! add data to accumulator
!    call abf_accu_add_data_online(cvhist(:,hist_len/2+1),pxi0,epot,erst,ekin,etot)
!
!end subroutine abf_core_force_sg

!===============================================================================
! subroutine:  tabf_core_calc_Zmat
!===============================================================================

subroutine tabf_core_calc_Zmat()

    use pmf_utils
    use tabf_dat

    implicit none
    integer         :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfTABFCVs
        ci = TABFCVList(i)%cvindx
        do j=1,NumOfTABFCVs
            cj = TABFCVList(j)%cvindx
            fz(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                fz(i,j) = fz(i,j) + MassInv(k)*dot_product(CVContext%CVsDrvs(:,k,ci),CVContext%CVsDrvs(:,k,cj))
            end do
            fzinv(i,j) = fz(i,j)            ! we need this for LAPACK
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfTABFCVs .gt. 1) then
        call dgetrf(NumOfTABFCVs,NumOfTABFCVs,fzinv,NumOfTABFCVs,indx,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[TABF] LU decomposition failed in tabf_calc_Zmat!')
        end if

        call dgetri(NumOfTABFCVs,fzinv,NumOfTABFCVs,indx,vv,NumOfTABFCVs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[TABF] Matrix inversion failed in tabf_calc_Zmat!')
        end if
    else
        fzinv(1,1)=1.0d0/fz(1,1)
    end if

    return

end subroutine tabf_core_calc_Zmat

!===============================================================================

end module tabf_core
