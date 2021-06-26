!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module rst_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  rst_accu_init
!===============================================================================

subroutine rst_accu_init()

    use rst_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! only if we would like to update restart file and of if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    ! init dimensions ------------------------------
    rstaccu%tot_cvs = NumOfRSTItems
    allocate(rstaccu%sizes(rstaccu%tot_cvs ), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[RST] Unable to allocate memory for rst accumulator!')
    endif

    tot_nbins = 1
    do i=1,rstaccu%tot_cvs
        rstaccu%sizes(i)%min_value  = RSTCVList(i)%min_value
        rstaccu%sizes(i)%max_value  = RSTCVList(i)%max_value
        rstaccu%sizes(i)%nbins      = RSTCVList(i)%nbins
        rstaccu%sizes(i)%width      = abs(rstaccu%sizes(i)%max_value - rstaccu%sizes(i)%min_value)
        rstaccu%sizes(i)%bin_width  = rstaccu%sizes(i)%width / rstaccu%sizes(i)%nbins
        rstaccu%sizes(i)%cv         => RSTCVList(i)%cv
        tot_nbins = tot_nbins * rstaccu%sizes(i)%nbins
    end do

    rstaccu%tot_nbins = tot_nbins

    ! RST force arrays
    allocate( rstaccu%nsamples(rstaccu%tot_nbins), &
              stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[RST] Unable to allocate memory for rst accumulator (rstforce)!')
    endif

    call rst_accu_clear()

    return

end subroutine rst_accu_init

!===============================================================================
! Subroutine:  rst_accu_clear
!===============================================================================

subroutine rst_accu_clear()

    use rst_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    rstaccu%nsamples(:) = 0

end subroutine rst_accu_clear

!===============================================================================
! Subroutine:  rst_accu_read
!===============================================================================

subroutine rst_accu_read(iounit)

    use rst_dat
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
            call pmf_accu_read_header(rstaccu%PMFAccuType,iounit,'RST',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_ibuf_B(rstaccu%PMFAccuType,iounit,keyline,rstaccu%nsamples)
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,RST_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to read from the accumulator - keyline!')

end subroutine rst_accu_read

!===============================================================================
! Subroutine:  rst_accu_write
!===============================================================================

subroutine rst_accu_write(iounit)

    use rst_dat
    use pmf_dat
    use pmf_utils
    use pmf_unit

    implicit none
    integer                :: iounit
    !---------------------------------------------------------------------------

    rstaccu%method = 'RST'
    call pmf_accu_write_header(rstaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(rstaccu%PMFAccuType,iounit,'NSAMPLES','AD',rstaccu%nsamples)

end subroutine rst_accu_write

!===============================================================================
! Subroutine:  rst_accu_add_sample
!===============================================================================

subroutine rst_accu_add_sample(values)

    use rst_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: values(:)
    ! -----------------------------------------------
    integer        :: gi0
    ! --------------------------------------------------------------------------

    ! only if we would like to update restart file and of if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(rstaccu%PMFAccuType,values)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    rstaccu%nsamples(gi0) = rstaccu%nsamples(gi0) + 1

end subroutine rst_accu_add_sample

!===============================================================================

end module rst_accu

