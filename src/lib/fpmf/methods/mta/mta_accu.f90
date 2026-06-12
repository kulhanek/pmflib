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

module mta_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mta_accu_init
!===============================================================================

subroutine mta_accu_init()

    use mta_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    mtaaccu%tot_cvs = NumOfMTACVs

    ! init dimensions ------------------------------
    allocate(mtaaccu%sizes(mtaaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTA] Unable to allocate memory for mta accumulator!')
    endif

    mtaaccu%tot_nbins   = 1
    do i=1, mtaaccu%tot_cvs
        mtaaccu%sizes(i)%min_value  = MTACVList(i)%min_value
        mtaaccu%sizes(i)%max_value  = MTACVList(i)%max_value
        mtaaccu%sizes(i)%nbins      = MTACVList(i)%nbins
        mtaaccu%sizes(i)%width      = abs(mtaaccu%sizes(i)%max_value - mtaaccu%sizes(i)%min_value)
        mtaaccu%sizes(i)%bin_width  = mtaaccu%sizes(i)%width / mtaaccu%sizes(i)%nbins
        mtaaccu%sizes(i)%cv         => MTACVList(i)%cv
        mtaaccu%tot_nbins           = mtaaccu%tot_nbins * mtaaccu%sizes(i)%nbins
    end do

    ! MTA arrays
    allocate(   mtaaccu%nsamples(mtaaccu%tot_nbins),                    &
                mtaaccu%mmta(mtaaccu%tot_nbins),        &
                mtaaccu%m2mta(mtaaccu%tot_nbins),       &
                mtaaccu%mimta(mtaaccu%tot_nbins),       &
                mtaaccu%m2imta(mtaaccu%tot_nbins),      &
                stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTA] Unable to allocate memory for mta accumulator (mtaforce)!')
    endif

    call mta_accu_clear()

    return

end subroutine mta_accu_init

!===============================================================================
! Subroutine:  mta_accu_clear
!===============================================================================

subroutine mta_accu_clear()

    use mta_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    mtaaccu%nsamples(:)     = 0.0d0
    mtaaccu%mmta(:)         = 0.0d0
    mtaaccu%m2mta(:)        = 0.0d0
    mtaaccu%mimta(:)        = 0.0d0
    mtaaccu%m2imta(:)       = 0.0d0

end subroutine mta_accu_clear

!===============================================================================
! Subroutine:  mta_accu_read
!===============================================================================

subroutine mta_accu_read(iounit)

    use mta_dat
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
            call pmf_accu_read_header(mtaaccu%PMFAccuType,iounit,'MTA',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_rbuf_B(mtaaccu%PMFAccuType,iounit,keyline,mtaaccu%nsamples)
            ! ------------------------------------
                case('MMTA')
                    call pmf_accu_read_rbuf_B(mtaaccu%PMFAccuType,iounit,keyline,mtaaccu%mmta)
            ! ------------------------------------
                case('M2MTA')
                    call pmf_accu_read_rbuf_B(mtaaccu%PMFAccuType,iounit,keyline,mtaaccu%m2mta)
            ! ------------------------------------
                case('MIMTA')
                    call pmf_accu_read_rbuf_B(mtaaccu%PMFAccuType,iounit,keyline,mtaaccu%mimta)
            ! ------------------------------------
                case('M2IMTA')
                    call pmf_accu_read_rbuf_B(mtaaccu%PMFAccuType,iounit,keyline,mtaaccu%m2imta)
            
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,MTA_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[MTA] Unable to read from the accumulator - keyline!')

end subroutine mta_accu_read

!===============================================================================
! Subroutine:  mta_accu_write
!===============================================================================

subroutine mta_accu_write(iounit)

    use mta_dat

    implicit none
    integer  :: iounit
    !---------------------------------------------------------------------------

    mtaaccu%method = 'MTA'
    call pmf_accu_write_header(mtaaccu%PMFAccuType,iounit)
    call pmf_accu_write_rbuf_B(mtaaccu%PMFAccuType,iounit,'NSAMPLES',   'AD',mtaaccu%nsamples)
    call pmf_accu_write_rbuf_B(mtaaccu%PMFAccuType,iounit,'MMTA',       'WA',mtaaccu%mmta,  'NSAMPLES')
    call pmf_accu_write_rbuf_B(mtaaccu%PMFAccuType,iounit,'M2MTA',      'M2',mtaaccu%m2mta, 'NSAMPLES','MMTA')
    call pmf_accu_write_rbuf_B(mtaaccu%PMFAccuType,iounit,'MIMTA',      'WA',mtaaccu%mimta, 'NSAMPLES')
    call pmf_accu_write_rbuf_B(mtaaccu%PMFAccuType,iounit,'M2IMTA',     'M2',mtaaccu%m2imta,'NSAMPLES','MIMTA')

end subroutine mta_accu_write

!===============================================================================
! Subroutine:  mta_accu_add_data_online
!===============================================================================

subroutine mta_accu_add_data_online

    use mta_dat
    use pmf_dat

    implicit none
    integer        :: gi0
    real(PMFDP)    :: invn
    real(PMFDP)    :: dmta1, dmta2
    real(PMFDP)    :: dimta1, dimta2
    real(PMFDP)    :: mta, imta
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(mtaaccu%PMFAccuType,CVContext%CVsValues(:))
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    mtaaccu%nsamples(gi0) = mtaaccu%nsamples(gi0) + 1.0d0
    invn = 1.0d0 / mtaaccu%nsamples(gi0)

    mta = sqrt(fzdet)
    imta = 1.0d0 / mta

    dmta1 = mta - mtaaccu%mmta(gi0)
    mtaaccu%mmta(gi0)  = mtaaccu%mmta(gi0)  + dmta1 * invn
    dmta2 = mta - mtaaccu%mmta(gi0)
    mtaaccu%m2mta(gi0) = mtaaccu%m2mta(gi0) + dmta1 * dmta2

    dimta1 = imta - mtaaccu%mimta(gi0)
    mtaaccu%mimta(gi0)  = mtaaccu%mimta(gi0)  + dimta1 * invn
    dimta2 = imta - mtaaccu%mimta(gi0)
    mtaaccu%m2imta(gi0) = mtaaccu%m2imta(gi0) + dimta1 * dimta2

end subroutine mta_accu_add_data_online

!===============================================================================

end module mta_accu

