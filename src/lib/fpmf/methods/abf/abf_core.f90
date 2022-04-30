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
            call abf_core_force_3pB
        case(2)
            call abf_core_force_gpr
        case(3)
            call abf_core_force_sg
        ! experimental algorithms
        case(10)
            call abf_core_force_3pA
        case(11)
            call abf_core_force_3pC
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main!')
    end select

    call abf_output_write
    call abf_trajectory_write_snapshot
    call abf_restart_update
    call abf_client_exchange_data(.false.)

end subroutine abf_core_main

!===============================================================================
! Subroutine:  abf_core_shake
! correct for forces from SHAKE
!===============================================================================

subroutine abf_core_shake

    use abf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    select case(fmode)
        case(10)
            ! fix total forces from SHAKE
            fhist(:,:,hist_len) =  fhist(:,:,hist_len) + SHAKEFrc(:,:)
        case(1,2,11)
            ! nothing to be done
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented fmode in abf_core_main_shake!')
    end select

end subroutine abf_core_shake

!===============================================================================
! Subroutine:  abf_core_force_3pA
! this is leap-frog ABF version, simplified algorithm
! employing forces from potential and SHAKE
!===============================================================================

subroutine abf_core_force_3pA()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m,ifac
    integer                :: ci,ki
    real(PMFDP)            :: v,v1,v2,f,etot,epot,erst,ekin
    ! --------------------------------------------------------------------------

    if( fdebug .and. (mod(fstep,5000) .eq. 0) ) then
        open(unit=4789,file='abf-fmode1.bicf',status='UNKNOWN')
        ifac = 10
        do i=1,ABFCVList(1)%nbins*ifac
            cvave(1) = ABFCVList(1)%min_value &
                     + real(i,PMFDP) * (ABFCVList(1)%max_value-ABFCVList(1)%min_value)/(ABFCVList(1)%nbins * ifac)
            la(:) = 0.0d0
            pxi0(:) = 0.0d0
            call abf_accu_get_data(cvave,pxi0)
            select case(feimode)
                case(0)
                    call abf_accu_get_data(cvave,la)
                case(1)
                    call abf_accu_get_data_lramp(cvave,la)
                case(2)
                    call pmf_timers_start_timer(PMFLIB_ABF_KS_TIMER)
                        call abf_accu_get_data_ksmooth(cvave,la)
                    call pmf_timers_stop_timer(PMFLIB_ABF_KS_TIMER)
                case(3)
                    call abf_accu_get_data_lsmooth(cvave,la)
                case default
                    call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented extrapolation/interpolation mode!')
            end select
            write(4789,*) cvave, pxi0, la, sfac
        end do
        close(4789)
    end if

! shift values
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        fhist(:,:,i)    = fhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

! update values
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)     = Vel(:,:)   ! in t-dt/2

    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt

! apply ABF force
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
    micfhist(:,hist_len)    = la(:)
    fhist(:,:,hist_len)     = Frc(:,:)   ! in t, SHAKE forces are added later

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
                zdhist(m,j,i,hist_len) = v
            end do
        end do
    end do

    ! record time progress of data
    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))

! ABF part
    if( fstep .ge. hist_len ) then

        do i=1,NumOfABFCVs
            f  = 0.0d0
            v1 = 0.0d0
            v2 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f = f + zdhist(m,j,i,hist_len-1)*fhist(m,j,hist_len-1)  * MassInv(j)
                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                end do
            end do
            pxi0(i) = f + 0.5d0*(v1+v2)*ifdtx
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)

        epot = epothist(hist_len-1)
        erst = ersthist(hist_len-1)
        ekin = ekinhist(hist_len-1)
        etot = epot + erst + ekin

        if( fdebug ) then
            write(DEBUG_ABF_FMODE1,*) fstep-1, cvhist(:,hist_len-1), pxi0, epot, erst, ekin, etot
        end if

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,ekin,etot)
    end if

    return

end subroutine abf_core_force_3pA

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
    real(PMFDP)            :: v,v1,v2,f,etot,epot,erst,ekin
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
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
                v = 0.0d0
                do k=1,NumOfABFCVs
                    ki = ABFCVList(k)%cvindx
                    v = v + fzinv(i,k)*CVContext%CVsDrvs(m,j,ki)
                end do
                zdhist(m,j,i,hist_len) = v
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

    ! record time progress of data
    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))

! ABF part
    if( fstep .ge. hist_len ) then
        do i=1,NumOfABFCVs
            f  = 0.0d0
            v1 = 0.0d0
            v2 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f = f + zdhist(m,j,i,hist_len-1)*(vhist(m,j,hist_len-0)-vhist(m,j,hist_len-1))
                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                    v2 = v2 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                end do
            end do
            pxi0(i) = (f + 0.5d0*(v1+v2)) * ifdtx
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)  ! unbiased estimate

        epot = epothist(hist_len-1)
        erst = ersthist(hist_len-1)
        ekin = ekinhist(hist_len-1)
        etot = epot + erst + ekin

        ! debug
        ! write(1225,*) epot,erst,ekin,etot

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,ekin,etot)
    end if

    return

end subroutine abf_core_force_3pB

!===============================================================================
! Subroutine:  abf_core_force_3pC
! this is leap-frog ABF version, simplified algorithm
! forces from positions
!===============================================================================

subroutine abf_core_force_3pC()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,f,etot,epot,erst,ekin
    ! --------------------------------------------------------------------------

! shift accuvalue history
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        xhist(:,:,i)    = xhist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    xhist(:,:,hist_len)    = Crd(:,:)

! shift epot ene
    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt

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
                zdhist(m,j,i,hist_len) = v
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

    ! record time progress of data
    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))


! ABF part
    if( fstep .ge. hist_len ) then
        do i=1,NumOfABFCVs
            f  = 0.0d0
            v  = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    f = f + zdhist(m,j,i,hist_len-1) &
                      * (xhist(m,j,hist_len-0)-2.0d0*xhist(m,j,hist_len-1)+xhist(m,j,hist_len-2))
                    ! velocity part
                    v = v + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-2)) &
                       * (xhist(m,j,hist_len-0)-xhist(m,j,hist_len-2))
                end do
            end do
            pxi0(i) = (f + 0.25d0*v) * ifdtx * ifdtx
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-1)  ! unbiased estimate

        epot = epothist(hist_len-1)
        erst = ersthist(hist_len-1)
        ekin = ekinhist(hist_len-1)
        etot = epot + erst + ekin

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,ekin,etot)
    end if

    return

end subroutine abf_core_force_3pC

!===============================================================================
! Subroutine:  abf_core_force_gpr
! this is leap-frog ABF version, simplified algorithm
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
    integer                     :: i,j,k
    integer                     :: ci,skip
    real(PMFDP)                 :: v,etot,epot,erst,ekin,mean
    character(len=PMF_MAX_TYPE) :: cvid
    ! --------------------------------------------------------------------------

! shift values
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

! update values
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
    end do

    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)
    fzinvhist(:,:,hist_len) = fzinv(:,:)

! apply ABF force
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
    micfhist(:,hist_len)    = la(:)

    ! record time progress of data
    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))


    ! use hist_len not gpr_len here, hist_len = gpr_len + 1 due to ekin delay
    if( mod(fstep,hist_len) .ne. 0 )  return

! calculate CV velocities
    do i=1,NumOfABFCVs
        ! calculate mean value
        mean = 0.0d0
        do k=1,gpr_len
            mean = mean + cvhist(i,k)
        end do
        mean = mean / real(gpr_len,PMFDP)

        ! shift data
        do k=1,gpr_len
            gpr_data(k) = cvhist(i,k) - mean
        end do

        ! solve GPR
        call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_cvs,gpr_len,gpr_data,1,0.0d0,gpr_model,1)

        do k=1,gpr_len
            ! calculate CV derivative in time - derivative is shift invariant
            xvelhist(i,k) = dot_product(gpr_model,gpr_kfd_cvs(:,k))
        end do

        if( fdebug ) then
            write(cvid,'(I0.3)') i
            open(unit=4789,file='abf-gpr.cvhist_'//trim(cvid),status='UNKNOWN')
            do k=1,gpr_len
                write(4789,*) k, dot_product(gpr_model,gpr_kff_cvs(:,k))+mean, cvhist(i,k)
            end do
            close(4789)
            open(unit=4789,file='abf-gpr.xvelhist_'//trim(cvid),status='UNKNOWN')
            write(4789,*) 1, xvelhist(i,1)
            do k=2,gpr_len-1
                write(4789,*) k, xvelhist(i,k), 0.5d0*(cvhist(i,k+1)-cvhist(i,k-1))*ifdtx
            end do
            write(4789,*) gpr_len, xvelhist(i,gpr_len)
            close(4789)
        end if
    end do

! calculate momenta
    do k=1,gpr_len
        do i=1,NumOfABFCVs
            v = 0.0d0
            do j=1,NumOfABFCVs
                v = v + fzinvhist(i,j,k) * xvelhist(j,k)
            end do
            xphist(i,k) = v
        end do
    end do

! calculate derivatives of CV momenta
    if( gpr_icf_cdf ) then
        do i=1,NumOfABFCVs
            do k=2,gpr_len-1
                icfhist(i,k) = 0.5d0*(xphist(i,k+1)-xphist(i,k-1))*ifdtx
            end do
        end do
    else
        do i=1,NumOfABFCVs

            ! setup data
            do k=1,gpr_len
                gpr_data(k) = xphist(i,k)
            end do

            ! solve GPR
            call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_icf,gpr_len,gpr_data,1,0.0d0,gpr_model,1)

            do k=1,gpr_len
                ! calculate momenta derivative in time - derivative is shift invariant
                icfhist(i,k) = dot_product(gpr_model,gpr_kfd_icf(:,k))

            end do

            if( fdebug ) then
                write(cvid,'(I0.3)') i
                open(unit=4789,file='abf-gpr.xphist_'//trim(cvid),status='UNKNOWN')
                do k=1,gpr_len
                    write(4789,*) k, dot_product(gpr_model,gpr_kff_icf(:,k)), xphist(i,k)
                end do
                close(4789)
                open(unit=4789,file='abf-gpr.icfhist_'//trim(cvid),status='UNKNOWN')
                write(4789,*) 1, icfhist(i,1)
                do k=2,gpr_len-1
                    write(4789,*) k, icfhist(i,k), 0.5d0*(xphist(i,k+1)-xphist(i,k-1))*ifdtx
                end do
                write(4789,*) gpr_len, icfhist(i,gpr_len)
                close(4789)
            end if
        end do
    end if
! record data
    if( gpr_ene_smooth .gt. 0 ) then
        select case(gpr_ene_smooth)
            case(1)
                mean = 0.0d0
                do k=1,gpr_len
                    mean = mean + epothist(k) + ersthist(k) + ekinhist(k)
                end do
                mean = mean / real(gpr_len,PMFDP)

                ! shift data
                do k=1,gpr_len
                    gpr_data(k) = epothist(k) + ersthist(k) + ekinhist(k) - mean
                end do
            case(2)
                mean = 0.0d0
                do k=1,gpr_len
                    mean = mean + ekinhist(k)
                end do
                mean = mean / real(gpr_len,PMFDP)

                ! shift data
                do k=1,gpr_len
                    gpr_data(k) = ekinhist(k) - mean
                end do
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported gpr_ene_smooth in abf_core_force_gpr!')
        end select

        ! solve GPR
        call dgemv('N',gpr_len,gpr_len,1.0d0,gpr_K_ene,gpr_len,gpr_data,1,0.0d0,gpr_model,1)

        select case(gpr_ene_smooth)
            case(1)
                open(unit=4789,file='abf-gpr.etot',status='UNKNOWN')
                do k=1,gpr_len
                    write(4789,*) k, dot_product(gpr_model,gpr_kff_ene(:,k))+mean, gpr_data(k)+mean
                end do
                close(4789)
            case(2)
                open(unit=4789,file='abf-gpr.ekin',status='UNKNOWN')
                do k=1,gpr_len
                    write(4789,*) k, dot_product(gpr_model,gpr_kff_ene(:,k))+mean, gpr_data(k)+mean
                end do
                close(4789)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported gpr_ene_smooth in abf_core_force_gpr!')
        end select
    end if

    skip = 0
    if( gpr_icf_cdf ) then
        skip = 1
    end if

    do k=1+gpr_boundary+skip,gpr_len-gpr_boundary-skip

        ! total ABF force
        pxi0(:) = icfhist(:,k) - micfhist(:,k)

        epot = epothist(k)
        erst = ersthist(k)

        select case(gpr_ene_smooth)
            case(0)
                ekin = ekinhist(k)
                etot = epot + erst + ekin
            case(1)
                etot = dot_product(gpr_model,gpr_kff_ene(:,k)) + mean
            case(2)
                ekin = dot_product(gpr_model,gpr_kff_ene(:,k)) + mean
                etot = epot + erst + ekin
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported gpr_ene_smooth in abf_core_force_gpr!')
        end select

        if( fdebug ) then
            write(DEBUG_ABF_FMODE4,*) fstep-hist_len+k, cvhist(:,k), pxi0, etot
        end if

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,k),pxi0,epot,erst,ekin,etot)
    end do

end subroutine abf_core_force_gpr

!===============================================================================
! Subroutine:  abf_core_force_sg
! this is leap-frog ABF version
! SG filter for differentiation
!===============================================================================

subroutine abf_core_force_sg()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i, j, ci
    real(PMFDP)            :: vp, vm, epot, erst, ekin, etot
    ! --------------------------------------------------------------------------

! shift values
    do i=1,hist_len-1
        cvhist(:,i)         = cvhist(:,i+1)
        epothist(i)         = epothist(i+1)
        ersthist(i)         = ersthist(i+1)
        ekinhist(i)         = ekinhist(i+1)
        fzinvhist(:,:,i)    = fzinvhist(:,:,i+1)
        micfhist(:,i)       = micfhist(:,i+1)
    end do

! update values
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
    end do

    epothist(hist_len)      = PotEne - fepotaverage
    ersthist(hist_len)      = PMFEne
    ekinhist(hist_len-1)    = KinEne - fekinaverage    ! shifted by -dt

! calculate Z matrix and its inverse
    call abf_core_calc_Zmat(CVContext)
    fzinvhist(:,:,hist_len) = fzinv(:,:)

! apply ABF force
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
    micfhist(:,hist_len)    = la(:)

    ! record time progress of data
    call abf_accu_add_data_record_lf(cvhist(:,hist_len),fzinv,la, &
                                     epothist(hist_len),ersthist(hist_len),ekinhist(hist_len-1))

    if( fstep .lt. hist_len ) return

    do i=1,NumOfABFCVs
        pxi0(i) = dot_product(sg_c1(:),cvhist(i,1:hist_len-2))
        pxi1(i) = dot_product(sg_c1(:),cvhist(i,3:hist_len-0))
    end do

! calculate momenta
    do i=1,NumOfABFCVs
        vm = 0.0d0
        vp = 0.0d0
        do j=1,NumOfABFCVs
            vm = vm + fzinvhist(i,j,hist_len/2+0) * pxi0(j)
            vp = vp + fzinvhist(i,j,hist_len/2+2) * pxi1(j)
        end do
        pxip(i) = vp
        pxim(i) = vm
    end do

! calculate derivatives of CV momenta
    do i=1,NumOfABFCVs
        pxi0(i) = 0.5d0*(pxip(i)-pxim(i))*ifdtx
    end do

    epot = epothist(hist_len/2+1)   ! hist_len is +2 bigger than fsgframelen
    erst = ersthist(hist_len/2+1)
    ekin = ekinhist(hist_len/2+1)
    etot = epot + erst + ekin

    pxi0(:)     = pxi0(:) - micfhist(:,hist_len/2+1)

    ! add data to accumulator
    call abf_accu_add_data_online(cvhist(hist_len/2+1,:),pxi0,epot,erst,ekin,etot)

end subroutine abf_core_force_sg

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
