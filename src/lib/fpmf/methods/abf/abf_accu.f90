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
                abfaccu%nsamples(abfaccu%tot_nbins), &
                abfaccu%micf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                abfaccu%m2icf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                abfaccu%mepot(abfaccu%tot_nbins), &
                abfaccu%m2epot(abfaccu%tot_nbins), &
                stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (abfforce)!')
    endif

    if( fentropy ) then
        allocate(   abfaccu%metot(abfaccu%tot_nbins), &
                    abfaccu%m2etot(abfaccu%tot_nbins), &
                    abfaccu%c11hh(abfaccu%tot_cvs,abfaccu%tot_nbins), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (entropy)!')
        endif
    end if

    if( fblock_size .gt. 0 ) then
        ! for block averages
        allocate(  &
            abfaccu%block_nsamples(abfaccu%tot_nbins), &
            abfaccu%block_micf(abfaccu%tot_cvs,abfaccu%tot_nbins), &
            abfaccu%block_mepot(abfaccu%tot_nbins), &
            stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (block_micf)!')
        endif

    end if

    if( fserver_enabled ) then
        ! ABF incremental force arrays
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
    ! --------------------------------------------------------------------------

    abfaccu%weights(:)      = 1.0d0
    abfaccu%nsamples(:)     = 0
    abfaccu%micf(:,:)       = 0.0d0
    abfaccu%m2icf(:,:)      = 0.0d0
    abfaccu%mepot(:)        = 0.0d0
    abfaccu%m2epot(:)       = 0.0d0

    if( fentropy ) then
        abfaccu%metot(:)        = 0.0d0
        abfaccu%m2etot(:)       = 0.0d0
        abfaccu%c11hh(:,:)      = 0.0d0
    end if

    if( fblock_size .gt. 0 ) then
        abfaccu%block_nsamples(:) = 0
        abfaccu%block_micf(:,:) = 0.0d0
        abfaccu%block_mepot(:) = 0.0d0
    end if

    if( fserver_enabled ) then
        abfaccu%inc_nsamples(:) = 0
        abfaccu%inc_micf(:,:) = 0.0d0
        abfaccu%inc_m2icf(:,:) = 0.0d0
        abfaccu%inc_mepot(:) = 0
        abfaccu%inc_m2epot(:) = 0
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
                    call pmf_accu_read_ibuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%nsamples)
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
                ! FIXME
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

subroutine abf_accu_write(iounit)

    use abf_dat

    implicit none
    integer                     :: iounit
    !---------------------------------------------------------------------------

    abfaccu%method = 'ABF'
    call pmf_accu_write_header(abfaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(abfaccu%PMFAccuType,iounit,'NSAMPLES','AD',abfaccu%nsamples)
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MICF','WA',abfaccu%micf)
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2ICF','AD',abfaccu%m2icf)

    if( fenthalpy ) then
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MEPOT','WA',abfaccu%mepot)
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2EPOT','AD',abfaccu%m2epot)
    end if

    if( fentropy ) then
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'METOT','WA',abfaccu%metot)
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2ETOT','AD',abfaccu%m2etot)
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11HH','AD',abfaccu%c11hh)
    end if

end subroutine abf_accu_write

!===============================================================================
! Subroutine:  abf_accu_add_data_online
!===============================================================================

subroutine abf_accu_add_data_online(cvs,gfx,epot,ekin,erst)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: ekin
    real(PMFDP)    :: erst
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, icf
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: detot1, detot2, etot
    real(PMFDP)    :: dicf1, dicf2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    abfaccu%nsamples(gi0) = abfaccu%nsamples(gi0) + 1
    invn = 1.0d0 / real(abfaccu%nsamples(gi0),PMFDP)

    ! potential energy
    depot1 = epot - abfaccu%mepot(gi0)
    abfaccu%mepot(gi0)  = abfaccu%mepot(gi0)  + depot1 * invn
    depot2 = epot - abfaccu%mepot(gi0)
    abfaccu%m2epot(gi0) = abfaccu%m2epot(gi0) + depot1 * depot2

    if( fentropy ) then
        etot = epot + ekin + erst
        ! potential energy
        detot1 = etot - abfaccu%metot(gi0)
        abfaccu%metot(gi0)  = abfaccu%metot(gi0)  + detot1 * invn
        detot2 = etot - abfaccu%metot(gi0)
        abfaccu%m2etot(gi0) = abfaccu%m2etot(gi0) + detot1 * detot2
    end if

    do i=1,abfaccu%tot_cvs
        icf = gfx(i)
        dicf1 = - icf - abfaccu%micf(i,gi0)
        abfaccu%micf(i,gi0)  = abfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = - icf -  abfaccu%micf(i,gi0)
        abfaccu%m2icf(i,gi0) = abfaccu%m2icf(i,gi0) + dicf1 * dicf2

        if( fentropy ) then
            abfaccu%c11hh(i,gi0)  = abfaccu%c11hh(i,gi0) + dicf1     * detot2
        end if

    end do


!    if( fserver_enabled ) then
!        abfaccu%inc_nsamples(gi0) = abfaccu%inc_nsamples(gi0) + 1
!        do i=1,NumOfABFCVs
!            a = gfx(i)
!            abfaccu%inc_micf(i,gi0)  = abfaccu%inc_micf(i,gi0) - a
!            abfaccu%inc_m2icf(i,gi0) = abfaccu%inc_m2icf(i,gi0) + a**2
!
!            a = a * lpotene
!            abfaccu%inc_icfmetot(i,gi0)  = abfaccu%inc_icfmetot(i,gi0) - a
!            abfaccu%inc_icfm2etot(i,gi0) = abfaccu%inc_icfm2etot(i,gi0) + a**2
!        end do
!        ! potential energy
!        abfaccu%inc_metot(gi0)  = abfaccu%inc_metot(gi0)  + lpotene
!        abfaccu%inc_m2etot(gi0) = abfaccu%inc_m2etot(gi0) + lpotene**2
!    end if

end subroutine abf_accu_add_data_online

!===============================================================================
! Subroutine:  abf_accu_add_data_blocked
!===============================================================================

subroutine abf_accu_add_data_blocked(cvs,gfx,epot,ekin,erst)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: ekin
    real(PMFDP)    :: erst
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, icf
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dicf1, dicf2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    end if
    insidesamples = insidesamples + 1

    abfaccu%block_nsamples(gi0)    = abfaccu%block_nsamples(gi0)    + 1
    invn = 1.0d0 / real(abfaccu%block_nsamples(gi0),PMFDP)

    ! potential energy
    depot1 = epot - abfaccu%block_mepot(gi0)
    abfaccu%block_mepot(gi0)  = abfaccu%block_mepot(gi0)  + depot1 * invn

    do i=1,abfaccu%tot_cvs
        icf = gfx(i)
        dicf1 = - icf - abfaccu%block_micf(i,gi0)
        abfaccu%block_micf(i,gi0)  = abfaccu%block_micf(i,gi0)  + dicf1 * invn
    end do

    ! is block filled?
    if( abfaccu%block_nsamples(gi0) .lt. fblock_size ) then
        ! NO -> exit
        return
    end if

   ! populate abf accumulator

    ! increase number of samples
    abfaccu%nsamples(gi0) = abfaccu%nsamples(gi0) + 1
    invn = 1.0d0 / real(abfaccu%nsamples(gi0),PMFDP)

    ! potential energy
    epot = abfaccu%block_mepot(gi0)
    depot1 = epot - abfaccu%mepot(gi0)
    abfaccu%mepot(gi0)  = abfaccu%mepot(gi0)  + depot1 * invn
    depot2 = epot - abfaccu%mepot(gi0)
    abfaccu%m2epot(gi0) = abfaccu%m2epot(gi0) + depot1 * depot2

    ! entropy cannot be calculated in blocked version

    do i=1,abfaccu%tot_cvs
        icf = abfaccu%block_micf(i,gi0)
        dicf1 = icf - abfaccu%micf(i,gi0)
        abfaccu%micf(i,gi0)  = abfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = icf -  abfaccu%micf(i,gi0)
        abfaccu%m2icf(i,gi0) = abfaccu%m2icf(i,gi0) + dicf1 * dicf2
    end do

!    ! increase number of samples
!    abfaccu%nsamples(gi0) = abfaccu%nsamples(gi0) + 1
!    do i=1,NumOfABFCVs
!        a = abfaccu%block_micf(i,gi0) / real(abfaccu%block_nsamples(gi0),PMFDP)
!        abfaccu%micf(i,gi0)  = abfaccu%micf(i,gi0)  - a
!        abfaccu%m2icf(i,gi0) = abfaccu%m2icf(i,gi0) + a**2
!
!        a = abfaccu%block_icfmetot(i,gi0) / real(abfaccu%block_nsamples(gi0),PMFDP)
!        abfaccu%icfmetot(i,gi0)  = abfaccu%icfmetot(i,gi0)  - a
!        abfaccu%icfm2etot(i,gi0) = abfaccu%icfm2etot(i,gi0) + a**2
!    end do
!    a = abfaccu%block_metot(gi0) / real(abfaccu%block_nsamples(gi0),PMFDP)
!    abfaccu%metot(gi0)  = abfaccu%metot(gi0) + a
!    abfaccu%m2etot(gi0) = abfaccu%m2etot(gi0) + a**2
!
!!    if( fserver_enabled ) then
!!        abfaccu%inc_nsamples(gi0) = abfaccu%inc_nsamples(gi0) + 1
!!        do i=1,NumOfABFCVs
!!            a = abfaccu%block_micf(i,gi0) / real(abfaccu%block_nsamples(gi0),PMFDP)
!!            abfaccu%inc_micf(i,gi0)  = abfaccu%inc_micf(i,gi0) - a
!!            abfaccu%inc_m2icf(i,gi0) = abfaccu%inc_m2icf(i,gi0) + a**2
!!
!!            a = abfaccu%block_icfmetot(i,gi0) / real(abfaccu%block_nsamples(gi0),PMFDP)
!!            abfaccu%inc_icfmetot(i,gi0)  = abfaccu%inc_icfmetot(i,gi0)  - a
!!            abfaccu%inc_icfm2etot(i,gi0) = abfaccu%inc_icfm2etot(i,gi0) + a**2
!!        end do
!!        a = abfaccu%block_metot(gi0) / real(abfaccu%block_nsamples(gi0),PMFDP)
!!        abfaccu%inc_metot(gi0)  = abfaccu%inc_metot(gi0) + a
!!        abfaccu%inc_m2etot(gi0) = abfaccu%inc_m2etot(gi0) + a**2
!!    end if
!
    ! reset the block
    abfaccu%block_nsamples(gi0)         = 0
    abfaccu%block_micf(:,gi0)           = 0.0d0
    abfaccu%block_mepot(gi0)            = 0.0d0

end subroutine abf_accu_add_data_blocked

!===============================================================================
! Subroutine:  abf_accu_get_data
!===============================================================================

subroutine abf_accu_get_data(values,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0
    real(PMFDP)    :: w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,values)
    if( gi0 .le. 0 ) return ! out of valid area

    w      = abfaccu%weights(gi0)
    gfx(:) = w * abfaccu%micf(:,gi0)

end subroutine abf_accu_get_data

!===============================================================================
! Subroutine:  abf_accu_get_data_lramp
!===============================================================================

subroutine abf_accu_get_data_lramp(values,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0, n
    real(PMFDP)    :: sc_ramp, w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,values)
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
        gfx(:) = w * sc_ramp * abfaccu%micf(:,gi0)
    end if

end subroutine abf_accu_get_data_lramp

!===============================================================================
! Subroutine:  abf_accu_get_data_gks
!===============================================================================

subroutine abf_accu_get_data_gks(values,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0
    real(PMFDP)    :: w
    ! --------------------------------------------------------------------------

    ! FIXME

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,values)
    if( gi0 .le. 0 ) return ! out of valid area

    w      = abfaccu%weights(gi0)
    gfx(:) = w * abfaccu%micf(:,gi0)

end subroutine abf_accu_get_data_gks

!===============================================================================

end module abf_accu

