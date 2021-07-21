!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_accu_init
!===============================================================================

subroutine abp_accu_init()

    use abp_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    allocate(abpaccu%sizes(NumOfABPCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABP] Unable to allocate memory for abp abpaccu!')
    endif

    abpaccu%tot_cvs     = NumOfABPCVs
    abpaccu%tot_nbins   = 1
    do i=1,NumOfABPCVs
        abpaccu%sizes(i)%min_value  = ABPCVList(i)%min_value
        abpaccu%sizes(i)%max_value  = ABPCVList(i)%max_value
        abpaccu%sizes(i)%nbins      = ABPCVList(i)%nbins
        abpaccu%sizes(i)%width      = abs(abpaccu%sizes(i)%max_value - abpaccu%sizes(i)%min_value)
        abpaccu%sizes(i)%bin_width  = abpaccu%sizes(i)%width / abpaccu%sizes(i)%nbins
        abpaccu%sizes(i)%cv         => ABPCVList(i)%cv
        abpaccu%tot_nbins           = abpaccu%tot_nbins * abpaccu%sizes(i)%nbins
    end do

    ! ABP force arrays
    allocate(   abpaccu%nsamples(abpaccu%tot_nbins),                &
                abpaccu%dpop(abpaccu%tot_cvs,abpaccu%tot_nbins),    &
                abpaccu%pop(abpaccu%tot_nbins),                     &
                abpaccu%binpos(abpaccu%tot_cvs,abpaccu%tot_nbins),  &
                abpaccu%widths(abpaccu%tot_cvs),                    &
                abpaccu%iwidths2(abpaccu%tot_cvs),                  &
                stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABP] Unable to allocate memory fo abpaccu!')
    endif

    if( fserver_enabled ) then
        ! ABP incremental arrays
        allocate(   abpaccu%inc_nsamples(abpaccu%tot_nbins), &
                    abpaccu%inc_dpop(abpaccu%tot_cvs,abpaccu%tot_nbins), &
                    abpaccu%inc_pop(abpaccu%tot_nbins), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABP] Unable to allocate memory for abpaccu (mwa)!')
        endif
    end if

    call abp_accu_clear()

end subroutine abp_accu_init

!===============================================================================
! Subroutine:  abp_accu_clear
!===============================================================================

subroutine abp_accu_clear()

    use abp_dat
    use pmf_dat

    implicit none
    integer         :: i
    ! --------------------------------------------------------------------------

    abpaccu%nsamples(:) = 0
    abpaccu%dpop(:,:)   = 0.0d0
    abpaccu%pop(:)      = 0.0d0 ! 1.0 is added explicitly in abp_accu_get_la_hramp
    abpaccu%m           = 1.0d0

    abpaccu%dnorm = 1.0d0  ! normalization factor for gaussian
    do i=1,abpaccu%tot_cvs
        abpaccu%widths(i)   = ABPCVList(i)%width
        abpaccu%iwidths2(i) = 1.0d0 / abpaccu%widths(i)**2
        abpaccu%dnorm = abpaccu%dnorm / (abpaccu%widths(i)*sqrt(2.0d0*PMF_PI))
    end do

    ! init binpos
    do i=1,abpaccu%tot_nbins
        ! get CV values on a grid point
        call pmf_accu_get_point(abpaccu%PMFAccuType,i,abpaccu%binpos(:,i))
    end do

    if( fserver_enabled ) then
        abpaccu%inc_nsamples(:) = 0
        abpaccu%inc_dpop(:,:)   = 0.0d0
        abpaccu%inc_pop(:)      = 0.0d0
    end if

end subroutine abp_accu_clear

!===============================================================================
! Subroutine:  abp_accu_read
!===============================================================================

subroutine abp_accu_read(iounit)

    use abp_dat
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
            call pmf_accu_read_header(abpaccu%PMFAccuType,iounit,'ABP',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_ibuf_B(abpaccu%PMFAccuType,iounit,keyline,abpaccu%nsamples)
            ! ------------------------------------
                case('DPOP')
                    call pmf_accu_read_rbuf_M(abpaccu%PMFAccuType,iounit,keyline,abpaccu%dpop)
            ! ------------------------------------
                case('POP')
                    call pmf_accu_read_rbuf_B(abpaccu%PMFAccuType,iounit,keyline,abpaccu%pop)
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,ABP_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[ABP] Unable to read from the accumulator - keyline!')

end subroutine abp_accu_read

!===============================================================================
! Subroutine:  abp_accu_write
!===============================================================================

subroutine abp_accu_write(iounit)

    use abp_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    integer                     :: iounit
    !---------------------------------------------------------------------------

    abpaccu%method = 'ABP'
    call pmf_accu_write_header(abpaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(abpaccu%PMFAccuType,iounit,'NSAMPLES','AD',abpaccu%nsamples)
    call pmf_accu_write_rbuf_M(abpaccu%PMFAccuType,iounit,'DPOP',    'WA',abpaccu%dpop)
    call pmf_accu_write_rbuf_B(abpaccu%PMFAccuType,iounit,'POP',     'WA',abpaccu%pop)
    call pmf_accu_write_rbuf_C(abpaccu%PMFAccuType,iounit,'WIDTHS',  'SA',abpaccu%widths)

end subroutine abp_accu_write

!===============================================================================
! Subroutine:  abp_accu_get_la_hramp
!===============================================================================

subroutine abp_accu_get_la_hramp

    use abp_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    integer        :: gi0
    real(PMFDP)    :: sc_ramp
    ! --------------------------------------------------------------------------

    la(:) = 0.0d0

    ! get global index to abpaccu for current CV values
    gi0 = pmf_accu_globalindex(abpaccu%PMFAccuType,cvvalues)

    if( gi0 .gt. 0 ) then
        sc_ramp = min(1.0d0, real(abpaccu%pop(gi0),PMFDP) / real(fhramp,PMFDP) )
        la(:) = sc_ramp * abpaccu%dpop(:,gi0) / (1.0d0 + abpaccu%pop(gi0))
    end if

end subroutine abp_accu_get_la_hramp

!===============================================================================
! Subroutine:  abp_accu_update_direct
!===============================================================================

subroutine abp_accu_update_direct

    use abp_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    integer     :: gi0,gi,i
    real(PMFDP) :: d,w,s,ew
    ! --------------------------------------------------------------------------

    ! get global index to abpaccu for current CV values
    gi0 = pmf_accu_globalindex(abpaccu%PMFAccuType,cvvalues)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return   ! out of range
    end if

    insidesamples = insidesamples + 1

    abpaccu%nsamples(gi0) = abpaccu%nsamples(gi0) + 1
    if( fserver_enabled ) then
        abpaccu%inc_nsamples(gi0) = abpaccu%inc_nsamples(gi0) + 1
    end if

    ! calculate W factor
    w = exp( cfac )*(abpaccu%pop(gi0)+1.0d0)/abpaccu%M

    ! direct mode - for each bin
    do gi=1,abpaccu%tot_nbins
        ! update pop
        s = 0.0d0
        do i=1,NumOfABPCVs
            ! calculate difference considering periodicity
            diffvalues(i) = ABPCVList(i)%cv%get_deviation(abpaccu%binpos(i,gi),cvvalues(i))
            ! argument for exp
            s = s + diffvalues(i)**2 * abpaccu%iwidths2(i)
        end do

        ew = w * abpaccu%dnorm * exp(- 0.5d0 * s)

        abpaccu%pop(gi) = abpaccu%pop(gi) + ew
        if( (abpaccu%pop(gi) + 1.0d0) .gt. abpaccu%M ) then
            abpaccu%M = abpaccu%pop(gi) + 1.0d0
        end if

        if( fserver_enabled ) then
            abpaccu%inc_pop(gi) = abpaccu%inc_pop(gi) + ew
        end if

        ! update dpop
        do i=1,NumOfABPCVs
            d = diffvalues(i) * abpaccu%iwidths2(i)
            abpaccu%dpop(i,gi) = abpaccu%dpop(i,gi) + d*ew
            if( fserver_enabled ) then
                abpaccu%inc_dpop(i,gi) = abpaccu%inc_dpop(i,gi) + d*ew
            end if
        end do

    end do

end subroutine abp_accu_update_direct

!===============================================================================

end module abp_accu

