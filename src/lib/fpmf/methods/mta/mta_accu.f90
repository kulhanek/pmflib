!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtc_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtc_accu_init
!===============================================================================

subroutine mtc_accu_init()

    use mtc_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    mtcaccu%tot_cvs = NumOfMTCCVs

    ! init dimensions ------------------------------
    allocate(mtcaccu%sizes(mtcaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTC] Unable to allocate memory for mtc accumulator!')
    endif

    mtcaccu%tot_nbins   = 1
    do i=1, mtcaccu%tot_cvs
        mtcaccu%sizes(i)%min_value  = MTCCVList(i)%min_value
        mtcaccu%sizes(i)%max_value  = MTCCVList(i)%max_value
        mtcaccu%sizes(i)%nbins      = MTCCVList(i)%nbins
        mtcaccu%sizes(i)%width      = abs(mtcaccu%sizes(i)%max_value - mtcaccu%sizes(i)%min_value)
        mtcaccu%sizes(i)%bin_width  = mtcaccu%sizes(i)%width / mtcaccu%sizes(i)%nbins
        mtcaccu%sizes(i)%cv         => MTCCVList(i)%cv
        mtcaccu%tot_nbins           = mtcaccu%tot_nbins * mtcaccu%sizes(i)%nbins
    end do

    ! MTC arrays
    allocate(   mtcaccu%nsamples(mtcaccu%tot_nbins),                    &
                mtcaccu%mmtc(mtcaccu%tot_nbins),        &
                mtcaccu%m2mtc(mtcaccu%tot_nbins),       &
                mtcaccu%mimtc(mtcaccu%tot_nbins),       &
                mtcaccu%m2imtc(mtcaccu%tot_nbins),      &
                stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTC] Unable to allocate memory for mtc accumulator (mtcforce)!')
    endif

    call mtc_accu_clear()

    return

end subroutine mtc_accu_init

!===============================================================================
! Subroutine:  mtc_accu_clear
!===============================================================================

subroutine mtc_accu_clear()

    use mtc_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    mtcaccu%nsamples(:)     = 0.0d0
    mtcaccu%mmtc(:)         = 0.0d0
    mtcaccu%m2mtc(:)        = 0.0d0
    mtcaccu%mimtc(:)        = 0.0d0
    mtcaccu%m2imtc(:)       = 0.0d0

end subroutine mtc_accu_clear

!===============================================================================
! Subroutine:  mtc_accu_read
!===============================================================================

subroutine mtc_accu_read(iounit)

    use mtc_dat
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
            call pmf_accu_read_header(mtcaccu%PMFAccuType,iounit,'MTC',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_rbuf_B(mtcaccu%PMFAccuType,iounit,keyline,mtcaccu%nsamples)
            ! ------------------------------------
                case('MMTC')
                    call pmf_accu_read_rbuf_B(mtcaccu%PMFAccuType,iounit,keyline,mtcaccu%mmtc)
            ! ------------------------------------
                case('M2MTC')
                    call pmf_accu_read_rbuf_B(mtcaccu%PMFAccuType,iounit,keyline,mtcaccu%m2mtc)
            ! ------------------------------------
                case('MIMTC')
                    call pmf_accu_read_rbuf_B(mtcaccu%PMFAccuType,iounit,keyline,mtcaccu%mimtc)
            ! ------------------------------------
                case('M2IMTC')
                    call pmf_accu_read_rbuf_B(mtcaccu%PMFAccuType,iounit,keyline,mtcaccu%m2imtc)
            
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,MTC_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[MTC] Unable to read from the accumulator - keyline!')

end subroutine mtc_accu_read

!===============================================================================
! Subroutine:  mtc_accu_write
!===============================================================================

subroutine mtc_accu_write(iounit)

    use mtc_dat

    implicit none
    integer  :: iounit
    !---------------------------------------------------------------------------

    mtcaccu%method = 'MTC'
    call pmf_accu_write_header(mtcaccu%PMFAccuType,iounit)
    call pmf_accu_write_rbuf_B(mtcaccu%PMFAccuType,iounit,'NSAMPLES',   'AD',mtcaccu%nsamples)
    call pmf_accu_write_rbuf_B(mtcaccu%PMFAccuType,iounit,'MMTC',       'WA',mtcaccu%mmtc,  'NSAMPLES')
    call pmf_accu_write_rbuf_B(mtcaccu%PMFAccuType,iounit,'M2MTC',      'M2',mtcaccu%m2mtc, 'NSAMPLES','MMTC')
    call pmf_accu_write_rbuf_B(mtcaccu%PMFAccuType,iounit,'MIMTC',      'WA',mtcaccu%mimtc, 'NSAMPLES')
    call pmf_accu_write_rbuf_B(mtcaccu%PMFAccuType,iounit,'M2IMTC',     'M2',mtcaccu%m2imtc,'NSAMPLES','MIMTC')

end subroutine mtc_accu_write

!===============================================================================
! Subroutine:  mtc_accu_add_data_online
!===============================================================================

subroutine mtc_accu_add_data_online

    use mtc_dat
    use pmf_dat

    implicit none
    integer        :: gi0
    real(PMFDP)    :: invn
    real(PMFDP)    :: dmtc1, dmtc2
    real(PMFDP)    :: dimtc1, dimtc2
    real(PMFDP)    :: mtc, imtc
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(mtcaccu%PMFAccuType,CVContext%CVsValues(:))
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    mtcaccu%nsamples(gi0) = mtcaccu%nsamples(gi0) + 1.0d0
    invn = 1.0d0 / mtcaccu%nsamples(gi0)

    mtc = sqrt(fzdet)
    imtc = 1.0d0 / mtc

    dmtc1 = mtc - mtcaccu%mmtc(gi0)
    mtcaccu%mmtc(gi0)  = mtcaccu%mmtc(gi0)  + dmtc1 * invn
    dmtc2 = mtc - mtcaccu%mmtc(gi0)
    mtcaccu%m2mtc(gi0) = mtcaccu%m2mtc(gi0) + dmtc1 * dmtc2

    dimtc1 = imtc - mtcaccu%mimtc(gi0)
    mtcaccu%mimtc(gi0)  = mtcaccu%mimtc(gi0)  + dimtc1 * invn
    dimtc2 = imtc - mtcaccu%mimtc(gi0)
    mtcaccu%m2imtc(gi0) = mtcaccu%m2imtc(gi0) + dimtc1 * dimtc2

end subroutine mtc_accu_add_data_online

!===============================================================================

end module mtc_accu

