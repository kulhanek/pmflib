!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module abf_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_accu_init
!===============================================================================

subroutine abf_accu_init()

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    allocate(abfaccu%sizes(NumOfABFCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator!')
    endif

    abfaccu%tot_cvs     = NumOfABFCVs
    abfaccu%tot_nbins   = 1
    do i=1, abfaccu%tot_cvs
        abfaccu%sizes(i)%min_value  = ABFCVList(i)%min_value
        abfaccu%sizes(i)%max_value  = ABFCVList(i)%max_value
        abfaccu%sizes(i)%nbins      = ABFCVList(i)%nbins
        abfaccu%sizes(i)%width      = abs(abfaccu%sizes(i)%max_value - abfaccu%sizes(i)%min_value)
        abfaccu%sizes(i)%bin_width  = abfaccu%sizes(i)%width / abfaccu%sizes(i)%nbins
        abfaccu%sizes(i)%cv         => ABFCVList(i)%cv
        abfaccu%tot_nbins           = abfaccu%tot_nbins * abfaccu%sizes(i)%nbins
    end do

    ! ABF force arrays
    allocate(   abfaccu%weights(abfaccu%tot_nbins), &
                abfaccu%binpos(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                abfaccu%nsamples(abfaccu%tot_nbins), &
                abfaccu%micf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                abfaccu%m2icf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                abfaccu%bnsamples(abfaccu%tot_nbins), &
                abfaccu%bmicf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (abfforce)!')
    endif

    if( fenthalpy ) then
        allocate(   abfaccu%mepot(abfaccu%tot_nbins), &
                    abfaccu%m2epot(abfaccu%tot_nbins), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (enthalpy)!')
        endif
    end if

    if( fentropy ) then
        allocate(   abfaccu%cvs(abfaccu%tot_cvs,fnstlim), &
                    abfaccu%zinv(abfaccu%tot_cvs,abfaccu%tot_cvs,fnstlim), &
                    abfaccu%icf(abfaccu%tot_cvs,fnstlim), &
                    abfaccu%abf(abfaccu%tot_cvs,fnstlim), &
                    abfaccu%epot(fnstlim), &
                    abfaccu%erst(fnstlim), &
                    abfaccu%ekin(fnstlim), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (entropy)!')
        endif
    end if

    if( fserver_enabled ) then
        ! ABF incremental force arrays - allocate all to simplify F90/C++ interface
        allocate(   abfaccu%inc_nsamples(abfaccu%tot_nbins), &
                    abfaccu%inc_micf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                    abfaccu%inc_m2icf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                    abfaccu%inc_mepot(abfaccu%tot_nbins), &
                    abfaccu%inc_m2epot(abfaccu%tot_nbins), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (inc_micf)!')
        endif
    end if

    call abf_accu_clear()

    return

end subroutine abf_accu_init

!===============================================================================
! Subroutine:  abf_accu_clear
!===============================================================================

subroutine abf_accu_clear()

    use abf_dat
    use pmf_dat

    implicit none
    integer :: i
    ! --------------------------------------------------------------------------

    ! init binpos
    do i=1,abfaccu%tot_nbins
        ! get CV values on a grid point
        call pmf_accu_get_point(abfaccu%PMFAccuType,i,abfaccu%binpos(:,i))
    end do

    abfaccu%weights(:)      = 1.0d0
    abfaccu%nsamples(:)     = 0
    abfaccu%micf(:,:)       = 0.0d0
    abfaccu%m2icf(:,:)      = 0.0d0

    abfaccu%bnsamples(:)    = 0
    abfaccu%bmicf(:,:)      = 0.0d0

    if( fenthalpy ) then
        abfaccu%mepot(:)        = 0.0d0
        abfaccu%m2epot(:)       = 0.0d0
    end if

    if( fentropy ) then
        abfaccu%cvs(:,:)        = 0.0d0
        abfaccu%zinv(:,:,:)     = 0.0d0
        abfaccu%icf(:,:)        = 0.0d0
        abfaccu%abf(:,:)        = 0.0d0
        abfaccu%epot(:)         = 0.0d0
        abfaccu%erst(:)         = 0.0d0
        abfaccu%ekin(:)         = 0.0d0
    end if

    if( fserver_enabled ) then
        abfaccu%inc_nsamples(:) = 0
        abfaccu%inc_micf(:,:)   = 0.0d0
        abfaccu%inc_m2icf(:,:)  = 0.0d0
        abfaccu%inc_mepot(:)    = 0.0d0
        abfaccu%inc_m2epot(:)   = 0.0d0
    end if

end subroutine abf_accu_clear

!===============================================================================
! Subroutine:  abf_accu_read
!===============================================================================

subroutine abf_accu_read(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    character(len=PMF_KEYLINE)      :: keyline
    ! --------------------------------------------------------------------------

    do while(.true.)

        ! read keyline
        read(iounit,5,end=500,err=300) keyline

        ! process keyline
        if( pmf_accu_is_header_key(keyline) ) then
            call pmf_accu_read_header(abfaccu%PMFAccuType,iounit,'ABF',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%nsamples)
            ! ------------------------------------
                case('MICF')
                    call pmf_accu_read_rbuf_M(abfaccu%PMFAccuType,iounit,keyline,abfaccu%micf)
            ! ------------------------------------
                case('M2ICF')
                    call pmf_accu_read_rbuf_M(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2icf)
            ! ------------------------------------
                case('MEPOT')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%mepot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2EPOT')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2epot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from the accumulator - keyline!')

end subroutine abf_accu_read

!===============================================================================
! Subroutine:  abf_accu_read_mask
!===============================================================================

subroutine abf_accu_read_mask(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    character(len=PMF_KEYLINE)      :: keyline
    ! --------------------------------------------------------------------------

    do while(.true.)

        ! read keyline
        read(iounit,5,end=500,err=300) keyline

        ! process keyline
        if( pmf_accu_is_header_key(keyline) ) then
            call pmf_accu_read_header(abfaccu%PMFAccuType,iounit,keyline,'ABF')
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('WEIGHTS')
                    call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%weights)
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[ABFMASK] Unable to read from the accumulator - keyline!')

end subroutine abf_accu_read_mask

!===============================================================================
! Subroutine:  abf_accu_write
!===============================================================================

subroutine abf_accu_write(iounit,full)

    use abf_dat

    implicit none
    integer  :: iounit
    logical  :: full
    !---------------------------------------------------------------------------

    abfaccu%method = 'ABF'
    call pmf_accu_write_header(abfaccu%PMFAccuType,iounit)
    call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'NSAMPLES',   'AD',abfaccu%nsamples)
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MICF',       'WA',abfaccu%micf,  'NSAMPLES')
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2ICF',      'M2',abfaccu%m2icf, 'NSAMPLES','MICF')

    if( fenthalpy ) then
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MEPOT',  'WA',abfaccu%mepot, 'NSAMPLES')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2EPOT', 'M2',abfaccu%m2epot,'NSAMPLES','MEPOT')
    end if

    call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'BNSAMPLES',  'IG',abfaccu%bnsamples)
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'BMICF',      'IG',abfaccu%bmicf, 'BNSAMPLES')

    if( full ) then
        if( fentropy ) then
            call pmf_accu_write_rbuf_TC(abfaccu%PMFAccuType,iounit,     'TCVS', abfaccu%cvs)
            call pmf_accu_write_rbuf_TCC(abfaccu%PMFAccuType,iounit,    'TZINV',abfaccu%zinv)
            call pmf_accu_write_rbuf_TC(abfaccu%PMFAccuType,iounit,     'TICF', abfaccu%icf)
            call pmf_accu_write_rbuf_TC(abfaccu%PMFAccuType,iounit,     'TMICF',abfaccu%abf)
            call pmf_accu_write_rbuf_T(iounit,                          'TEPOT',abfaccu%epot)
            call pmf_accu_write_rbuf_T(iounit,                          'TERST',abfaccu%erst)
            call pmf_accu_write_rbuf_T(iounit,                          'TEKIN',abfaccu%ekin)
        end if
    end if

end subroutine abf_accu_write

!===============================================================================
! Subroutine:  abf_accu_add_data_online
!===============================================================================

subroutine abf_accu_add_data_online(cvs,gfx,epot)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)    ! ICF
    real(PMFDP)    :: epot
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, icf
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dicf1, dicf2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    abfaccu%nsamples(gi0) = abfaccu%nsamples(gi0) + 1.0d0
    invn = 1.0d0 / abfaccu%nsamples(gi0)

    if( fenthalpy ) then
        ! potential energy
        depot1 = epot - abfaccu%mepot(gi0)
        abfaccu%mepot(gi0)  = abfaccu%mepot(gi0)  + depot1 * invn
        depot2 = epot - abfaccu%mepot(gi0)
        abfaccu%m2epot(gi0) = abfaccu%m2epot(gi0) + depot1 * depot2
    end if

    do i=1,NumOfABFCVs
        icf = - gfx(i)
        dicf1 = icf - abfaccu%micf(i,gi0)
        abfaccu%micf(i,gi0)  = abfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = icf - abfaccu%micf(i,gi0)
        abfaccu%m2icf(i,gi0) = abfaccu%m2icf(i,gi0) + dicf1 * dicf2
    end do

    if( fserver_enabled ) then
        abfaccu%inc_nsamples(gi0) = abfaccu%inc_nsamples(gi0) + 1.0d0
        invn = 1.0d0 / abfaccu%inc_nsamples(gi0)

        if( fenthalpy ) then
            depot1 = epot - abfaccu%inc_mepot(gi0)
            abfaccu%inc_mepot(gi0)  = abfaccu%inc_mepot(gi0)  + depot1 * invn
            depot2 = epot - abfaccu%inc_mepot(gi0)
            abfaccu%inc_m2epot(gi0) = abfaccu%inc_m2epot(gi0) + depot1 * depot2
        end if

        do i=1,NumOfABFCVs
            icf = - gfx(i)
            dicf1 = icf - abfaccu%inc_micf(i,gi0)
            abfaccu%inc_micf(i,gi0)  = abfaccu%inc_micf(i,gi0)  + dicf1 * invn
            dicf2 = icf -  abfaccu%inc_micf(i,gi0)
            abfaccu%inc_m2icf(i,gi0) = abfaccu%inc_m2icf(i,gi0) + dicf1 * dicf2
        end do

    end if

    ! increase number of samples for applied bias - this uses different counter for number of samples
    abfaccu%bnsamples(gi0) = abfaccu%bnsamples(gi0) + 1.0d0
    invn = 1.0d0 / abfaccu%bnsamples(gi0)

    do i=1,NumOfABFCVs
        icf = - gfx(i)
        dicf1 = icf - abfaccu%bmicf(i,gi0)
        abfaccu%bmicf(i,gi0)  = abfaccu%bmicf(i,gi0)  + dicf1 * invn
    end do

end subroutine abf_accu_add_data_online

!===============================================================================
! Function:  abf_get_skernel
!===============================================================================

function abf_get_skernel(cvs1,cvs2) result(kval)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs1(:)
    real(PMFDP)    :: cvs2(:)
    real(PMFDP)    :: kval
    ! -----------------------------------------------
    integer        :: i
    real(PMFDP)    :: dx,u2
    ! --------------------------------------------------------------------------

! calculate length between two points
    u2 = 0.0d0
    do i=1,NumOfABFCVs
        dx = (cvs1(i) - cvs2(i))/(ABFCVList(i)%wfac*abfaccu%PMFAccuType%sizes(i)%bin_width)
        u2 = u2 + dx**2
    end do

! get value of kernel
! https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
    kval = 0.0d0
    select case(fsmooth_kernel)
    ! Epanechnikov (parabolic)
        case(0)
            if( u2 .lt.  1.0d0 ) then
                kval = 3.0d0/4.0d0*(1-u2)
            end if
     ! Triweight
        case(1)
            if( u2 .lt.  1.0d0 ) then
                kval = 35.0d0/32.0d0*(1-u2)**3
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented smoothing kernel in abf_get_skernel!')
    end select

end function abf_get_skernel

!===============================================================================
! Subroutine:  abf_accu_add_data_ksmooth
!===============================================================================

subroutine abf_accu_add_data_ksmooth(cvs,gfx,epot)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)    ! ICF
    real(PMFDP)    :: epot
    ! -----------------------------------------------
    integer        :: gi0, si0, i, ni
    real(PMFDP)    :: w, stot_weight, invn, icf
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dicf1, dicf2
    ! --------------------------------------------------------------------------

! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

! first calculate kernel values
    stot_weight = 0.0d0
    sweights(:) = 0.0d0
    do ni=1,max_snb_size
        if( snb_list(ni,gi0) .le. 0 ) cycle
        w = abf_get_skernel(abfaccu%binpos(:,snb_list(ni,gi0)),cvs(:))
        stot_weight  = stot_weight + w
        sweights(ni) = w
    end do
    if( stot_weight .eq. 0.0d0 ) return ! out of valid area

! normalize weights
    do ni=1,max_snb_size
        if( snb_list(ni,gi0) .le. 0 ) cycle
        sweights(ni) = sweights(ni) / stot_weight
    end do

! apply distributed sample
    do ni=1,max_snb_size
        si0 = snb_list(ni,gi0)
        if( si0 .le. 0 ) cycle
        w = sweights(ni)
        if( sweights(ni) .le. 0.0d0 ) cycle

        ! increase number of samples
        abfaccu%nsamples(si0) = abfaccu%nsamples(si0) + w
        invn = w / abfaccu%nsamples(si0)

        if( fenthalpy ) then
            ! potential energy
            depot1 = epot - abfaccu%mepot(si0)
            abfaccu%mepot(si0)  = abfaccu%mepot(si0)  + depot1 * invn
            depot2 = epot - abfaccu%mepot(si0)
            abfaccu%m2epot(si0) = abfaccu%m2epot(si0) + w * depot1 * depot2
        end if

        do i=1,NumOfABFCVs
            icf = - gfx(i)
            dicf1 = icf - abfaccu%micf(i,si0)
            abfaccu%micf(i,si0)  = abfaccu%micf(i,si0)  + dicf1 * invn
            dicf2 = icf - abfaccu%micf(i,si0)
            abfaccu%m2icf(i,si0) = abfaccu%m2icf(i,si0) + w * dicf1 * dicf2
        end do

        if( fserver_enabled ) then
            abfaccu%inc_nsamples(si0) = abfaccu%inc_nsamples(si0) + w
            invn = w / abfaccu%inc_nsamples(si0)

            if( fenthalpy ) then
                depot1 = epot - abfaccu%inc_mepot(si0)
                abfaccu%inc_mepot(si0)  = abfaccu%inc_mepot(si0)  + depot1 * invn
                depot2 = epot - abfaccu%inc_mepot(si0)
                abfaccu%inc_m2epot(si0) = abfaccu%inc_m2epot(si0) + w * depot1 * depot2
            end if

            do i=1,NumOfABFCVs
                icf = - gfx(i)
                dicf1 = icf - abfaccu%inc_micf(i,si0)
                abfaccu%inc_micf(i,si0)  = abfaccu%inc_micf(i,si0)  + dicf1 * invn
                dicf2 = icf -  abfaccu%inc_micf(i,si0)
                abfaccu%inc_m2icf(i,si0) = abfaccu%inc_m2icf(i,si0) + w * dicf1 * dicf2
            end do

        end if

        ! increase number of samples for applied bias - this uses different counter for number of samples
        abfaccu%bnsamples(si0) = abfaccu%bnsamples(si0) + w
        invn = w / abfaccu%bnsamples(si0)

        do i=1,NumOfABFCVs
            icf = - gfx(i)
            dicf1 = icf - abfaccu%bmicf(i,si0)
            abfaccu%bmicf(i,si0)  = abfaccu%bmicf(i,si0)  + dicf1 * invn
        end do
    end do

end subroutine abf_accu_add_data_ksmooth

!===============================================================================
! Subroutine:  abf_accu_add_data_entropy
!===============================================================================

subroutine abf_accu_add_data_entropy(cvs,zinv,gfx,abf,epot,erst,ekin)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: zinv(:,:)
    real(PMFDP)    :: gfx(:)    ! ICF
    real(PMFDP)    :: abf(:)    ! MICF
    real(PMFDP)    :: epot
    real(PMFDP)    :: erst
    real(PMFDP)    :: ekin
    ! --------------------------------------------------------------------------

    if( fentropy ) then
        abfaccu%cvs(:,fstep)    = cvs(:)
        abfaccu%zinv(:,:,fstep) = zinv(:,:)
        abfaccu%icf(:,fstep)    = gfx(:)
        abfaccu%abf(:,fstep)    = abf(:)
        abfaccu%epot(fstep)     = epot
        abfaccu%erst(fstep)     = erst
        abfaccu%ekin(fstep)     = ekin
    end if

end subroutine abf_accu_add_data_entropy

!===============================================================================
! Subroutine:  abf_accu_get_data
!===============================================================================

subroutine abf_accu_get_data(cvs,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0
    real(PMFDP)    :: w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out of valid area

    w      = abfaccu%weights(gi0)
    gfx(:) = w * abfaccu%bmicf(:,gi0)

end subroutine abf_accu_get_data

!===============================================================================
! Subroutine:  abf_accu_get_data_lramp
!===============================================================================

subroutine abf_accu_get_data_lramp(cvs,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0
    real(PMFDP)    :: n, sc_ramp, w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get number of samples
    n = abfaccu%nsamples(gi0)
    if( n .gt. 0 ) then
        sc_ramp = 1.0d0
        if( n .le. fhramp_max ) then
            sc_ramp = 0.0d0
            if( n .gt. fhramp_min ) then
                sc_ramp = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        w      = abfaccu%weights(gi0)
        gfx(:) = w * sc_ramp * abfaccu%bmicf(:,gi0)
    end if

end subroutine abf_accu_get_data_lramp

!===============================================================================
! Subroutine:  abf_accu_get_data_ksmooth
!===============================================================================

subroutine abf_accu_get_data_ksmooth(cvs,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0,si0,ni
    real(PMFDP)    :: stot_weight,w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out of valid area

! first calculate kernel values
    sweights(:) = 0.0d0
    stot_weight = 0.0d0
    do ni=1,max_snb_size
        if( snb_list(ni,gi0) .le. 0 ) cycle
        w = abf_get_skernel(abfaccu%binpos(:,snb_list(ni,gi0)),cvs(:))
        stot_weight  = stot_weight + w
        sweights(ni) = w
    end do
    if( stot_weight .eq. 0.0d0 ) return ! out of valid area

! normalize weights
    do ni=1,max_snb_size
        if( snb_list(ni,gi0) .le. 0 ) cycle
        sweights(ni) = sweights(ni) / stot_weight
    end do

    write(84621,*) sweights

! get smoothed mean forces
    do ni=1,max_snb_size
        si0 = snb_list(ni,gi0)
        if( si0 .le. 0 ) cycle
        w = sweights(ni)
        if( sweights(ni) .le. 0.0d0 ) cycle

        w = w * abfaccu%weights(si0)
        gfx(:) = gfx(:) + w * abfaccu%bmicf(:,si0)
    end do

end subroutine abf_accu_get_data_ksmooth

!===============================================================================

end module abf_accu

