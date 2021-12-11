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
        case(1)
            call abf_core_force_3pA
        case(2)
            call abf_core_force_3pB
        case(3)
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
! Subroutine:  abf_core_shake
! correct for forces from SHAKE
!===============================================================================

subroutine abf_core_shake

    use abf_dat
    use pmf_utils

    ! --------------------------------------------------------------------------

    select case(fmode)
        case(1)
            ! fix total forces from SHAKE
            fhist(:,:,hist_len) =  fhist(:,:,hist_len) + SHAKEFrc(:,:)
        case(2)
            ! nothing to be done
        case(3)
            ! record SHAKE forces
            fshist(:,:,hist_len) = SHAKEFrc(:,:)
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
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,v1,v2,f,etot,epot,erst,ekin
    ! --------------------------------------------------------------------------

! shift values
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        fhist(:,:,i)    = fhist(:,:,i+1)
        fshist(:,:,i)   = fshist(:,:,i+1)
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

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,etot)
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

! ABF part
    if( fstep .ge. hist_len ) then
        do i=1,NumOfABFCVs
            f  = 0.0d0
            v  = 0.0d0
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

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,hist_len-1),pxi0,epot,erst,etot)
    end if

    return

end subroutine abf_core_force_3pB

!===============================================================================
! Subroutine:  abf_core_force_gpr
!===============================================================================

subroutine abf_core_force_gpr()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki,gi
    real(PMFDP)            :: v,v1,fp,fs,etot,epot,erst,ekin
    ! --------------------------------------------------------------------------

! shift values
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        fhist(:,:,i)    = fhist(:,:,i+1)
        fshist(:,:,i)   = fshist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

! update values
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len)  = CVContext%CVsValues(ci)
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

! ABF part
    if( fstep .ge. 2 ) then
        ! move history
        do i=1,gpr_len-1
            fpgprhist(:,i) = fpgprhist(:,i+1)
            fsgprhist(:,i) = fsgprhist(:,i+1)
            v1gprhist(:,i) = v1gprhist(:,i+1)
        end do

        do i=1,NumOfABFCVs
            fp = 0.0d0
            fs = 0.0d0
            v1 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    fp = fp + zdhist(m,j,i,hist_len-1)*fhist(m,j,hist_len-1)  * MassInv(j)
                    fs = fs + zdhist(m,j,i,hist_len-1)*fshist(m,j,hist_len-1) * MassInv(j)
                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-0)-zdhist(m,j,i,hist_len-1)) * vhist(m,j,hist_len-0)
                end do
            end do
            fpgprhist(i,gpr_len) = fp
            fsgprhist(i,gpr_len) = fs
            v1gprhist(i,gpr_len) = v1 * ifdtx
        end do
    end if

! GPR part
    if( fstep .ge. hist_len ) then

        gi = gpr_len/2 + 1

        pxi0(:) = fpgprhist(:,gi) + fsgprhist(:,gi) ! + (v1gprhist(:,gi)+v1gprhist(:,gi-1))*0.5d0

        do i=1,NumOfABFCVs
            ! shift data
            gpr_model(:) = v1gprhist(i,:)

            ! solve GPR
            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)

            if( gpr_info .ne. 0 ) then
                ! throw error
                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model in abf_core_force_gpr!')
            end if

            ! predict
            v1 = dot_product(gpr_model,gpr_kffhm)
            pxi0(i) = pxi0(i) + v1

            ! write(458,*) v1, (v1gprhist(:,gi)+v1gprhist(:,gi-1))*0.5d0

!            ! shift data
!            gpr_model(:) = fsgprhist(i,:)
!
!            ! solve GPR
!            call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)
!
!            if( gpr_info .ne. 0 ) then
!                ! throw error
!                call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model in abf_core_force_gpr!')
!            end if
!
!            ! predict
!            v1 = dot_product(gpr_model,gpr_kff)
!            pxi0(i) = pxi0(i) + v1

            ! write(459,*) v1,  fsgprhist(:,gi)
        end do

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,gi)

        epot = epothist(gi)
        erst = ersthist(gi)
        ekin = ekinhist(gi)

!        gpr_model(:) = ekinhist(1:gpr_len)
!
!        ! solve GPR
!        call dgetrs('N',gpr_len,1,gpr_K,gpr_len,gpr_indx,gpr_model,gpr_len,gpr_info)
!
!        if( gpr_info .ne. 0 ) then
!            ! throw error
!            call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to solve GPR model in abf_core_force_gpr!')
!        end if
!
!        ! predict
!        ekin = dot_product(gpr_model,gpr_kff)

        etot = epot + erst + 0.93d0*ekin

        ! add data to accumulator
        call abf_accu_add_data_online(cvhist(:,gi),pxi0,epot,erst,etot)
    end if

    return

end subroutine abf_core_force_gpr

!===============================================================================
! Subroutine:  abf_core_force_6p
! this is leap-frog ABF version, simplified algorithm
! employing forces from potential and SHAKE
!===============================================================================

subroutine abf_core_force_6p()

    use pmf_utils
    use pmf_dat
    use pmf_cvs
    use abf_dat
    use abf_accu
    use pmf_timers

    implicit none
    integer                :: i,j,k,m
    integer                :: ci,ki
    real(PMFDP)            :: v,v1,v2,v3,v4,fp,fs,asf,etot,epot,erst,ekin
    ! --------------------------------------------------------------------------

! shift values
    do i=1,hist_len-1
        cvhist(:,i)     = cvhist(:,i+1)
        epothist(i)     = epothist(i+1)
        ersthist(i)     = ersthist(i+1)
        ekinhist(i)     = ekinhist(i+1)
        vhist(:,:,i)    = vhist(:,:,i+1)
        fhist(:,:,i)    = fhist(:,:,i+1)
        fshist(:,:,i)   = fshist(:,:,i+1)
        zdhist(:,:,:,i) = zdhist(:,:,:,i+1)
        micfhist(:,i)   = micfhist(:,i+1)
    end do

! update values
    do i=1,NumOfABFCVs
        ci = ABFCVList(i)%cvindx
        cvhist(i,hist_len) = CVContext%CVsValues(ci)
    end do
    vhist(:,:,hist_len)    = Vel(:,:)   ! in t-dt/2

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

! ABF part
    if( fstep .ge. hist_len ) then

        do i=1,NumOfABFCVs
            fp = 0.0d0
            fs = 0.0d0
            v1 = 0.0d0
            v2 = 0.0d0
            v3 = 0.0d0
            v4 = 0.0d0
            do j=1,NumOfLAtoms
                do m=1,3
                    ! force part
                    fp = fp + zdhist(m,j,i,hist_len-3) * fhist(m,j,hist_len-3)  * MassInv(j)

!                    asf = ( -3.0d0*fshist(m,j,hist_len-1) &
!                           +12.0d0*fshist(m,j,hist_len-2) &
!                           +17.0d0*fshist(m,j,hist_len-3) &
!                           +12.0d0*fshist(m,j,hist_len-4) &
!                            -3.0d0*fshist(m,j,hist_len-5))/35.0d0

                    asf = ( fshist(m,j,hist_len-2) &
                           +fshist(m,j,hist_len-3) &
                           +fshist(m,j,hist_len-4))/3.0d0

                    fs = fs + zdhist(m,j,i,hist_len-3) * asf * MassInv(j)

                    ! velocity part
                    v1 = v1 + (zdhist(m,j,i,hist_len-1)-zdhist(m,j,i,hist_len-2)) * vhist(m,j,hist_len-1)
                    v2 = v2 + (zdhist(m,j,i,hist_len-2)-zdhist(m,j,i,hist_len-3)) * vhist(m,j,hist_len-2)
                    v3 = v3 + (zdhist(m,j,i,hist_len-3)-zdhist(m,j,i,hist_len-4)) * vhist(m,j,hist_len-3)
                    v4 = v4 + (zdhist(m,j,i,hist_len-4)-zdhist(m,j,i,hist_len-5)) * vhist(m,j,hist_len-4)
                end do
            end do
            pxi0(i) = fp + fs + 0.25d0*(v1+v2+v3+v4)*ifdtx
        end do

        ! write(7894,*) fp, fs

        ! total ABF force
        pxi0(:) = pxi0(:) - micfhist(:,hist_len-3)  ! unbiased estimate

        epot = epothist(hist_len-3)
        erst = ersthist(hist_len-3)
        ekin = ekinhist(hist_len-3)
        etot = epot + erst + ekin

        ! add data to accumulator
 ! FIXME
 !       call abf_accu_add_data_online(cvhist(:,hist_len-3),pxi0,epot,erst,etot)
    end if

    return

end subroutine abf_core_force_6p

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
