!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtd_accu

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_accu_init
!===============================================================================

subroutine mtd_accu_init()

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    integer              :: j,k,multibins
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    mtdaccu%tot_cvs = NumOfMTDCVs
    allocate(mtdaccu%sizes(mtdaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd accumulator!')
    endif

    tot_nbins = 1
    do i=1,mtdaccu%tot_cvs
        mtdaccu%sizes(i)%min_value  = MTDCVList(i)%min_value
        mtdaccu%sizes(i)%max_value  = MTDCVList(i)%max_value
        mtdaccu%sizes(i)%nbins      = MTDCVList(i)%nbins
        mtdaccu%sizes(i)%width      = abs(mtdaccu%sizes(i)%max_value - mtdaccu%sizes(i)%min_value)
        mtdaccu%sizes(i)%bin_width  = mtdaccu%sizes(i)%width / mtdaccu%sizes(i)%nbins
        mtdaccu%sizes(i)%cv         => MTDCVList(i)%cv
        tot_nbins = tot_nbins * MTDCVList(i)%nbins
    end do

    mtdaccu%tot_nbins = tot_nbins

    ! MTD potential and force array
    allocate( mtdaccu%nsamples(mtdaccu%tot_nbins), &
              mtdaccu%binpositions(mtdaccu%tot_cvs,mtdaccu%tot_nbins), &
              mtdaccu%mtdpotential(mtdaccu%tot_nbins), &
              mtdaccu%mtdforce(mtdaccu%tot_cvs,mtdaccu%tot_nbins), &
              stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd accumulator!')
    endif

    do i=1,mtdaccu%tot_cvs
        multibins = 1

        do j=i+1,mtdaccu%tot_cvs
            multibins = multibins * mtdaccu%sizes(j)%nbins
        end do

        do j=1,mtdaccu%tot_nbins
            k = mod(((j-1)/multibins),mtdaccu%sizes(i)%nbins)
            mtdaccu%binpositions(i,j) = mtdaccu%sizes(i)%min_value + &
                                            real(k)*mtdaccu%sizes(i)%bin_width + mtdaccu%sizes(i)%bin_width / 2.0d0
        end do
    end do

    if( fserver_enabled ) then
        allocate( mtdaccu%inc_nsamples(mtdaccu%tot_nbins), &
                  mtdaccu%inc_mtdpotential(mtdaccu%tot_nbins), &
                  mtdaccu%inc_mtdforce(mtdaccu%tot_cvs,mtdaccu%tot_nbins), &
                  stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd accumulator (MWA)!')
        endif
    end if

    call mtd_accu_clear()

    return

end subroutine mtd_accu_init

!===============================================================================
! Subroutine:  mtd_accu_clear
!===============================================================================

subroutine mtd_accu_clear()

    use mtd_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    mtdaccu%nsamples(:)     = 0.0d0
    mtdaccu%mtdpotential(:) = 0.0d0
    mtdaccu%mtdforce(:,:)   = 0.0d0

    if( fserver_enabled ) then
        mtdaccu%inc_nsamples(:) = 0
        mtdaccu%inc_mtdpotential(:) = 0.0d0
        mtdaccu%inc_mtdforce(:,:) = 0.0d0
    end if

end subroutine mtd_accu_clear

!===============================================================================
! Subroutine:  mtd_accu_read
!===============================================================================

subroutine mtd_accu_read(iounit)

    use mtd_dat
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
            call pmf_accu_read_header(mtdaccu%PMFAccuType,iounit,'MTD',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_ibuf_B(mtdaccu%PMFAccuType,iounit,keyline,mtdaccu%nsamples)
            ! ------------------------------------
                case('MTDPOT')
                    call pmf_accu_read_rbuf_B(mtdaccu%PMFAccuType,iounit,keyline,mtdaccu%mtdpotential)
            ! ------------------------------------
                case('MTDFORCE')
                    call pmf_accu_read_rbuf_M(mtdaccu%PMFAccuType,iounit,keyline,mtdaccu%mtdforce)
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,MTD_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from the accumulator - keyline!')

end subroutine mtd_accu_read

!===============================================================================
! Subroutine:  mtd_accu_write
!===============================================================================

subroutine mtd_accu_write(iounit)

    use mtd_dat
    use pmf_dat
    use pmf_utils
    use pmf_unit

    implicit none
    integer                     :: iounit
    !---------------------------------------------------------------------------

    mtdaccu%method = 'MTD'
    call pmf_accu_write_header(mtdaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(mtdaccu%PMFAccuType,iounit,'NSAMPLES','AD',mtdaccu%nsamples)
    call pmf_accu_write_rbuf_B(mtdaccu%PMFAccuType,iounit,'MTDPOT','WA',mtdaccu%mtdpotential)
    call pmf_accu_write_rbuf_M(mtdaccu%PMFAccuType,iounit,'MTDFORCE','WA',mtdaccu%mtdforce)

    return

end subroutine mtd_accu_write

!===============================================================================
! Subroutine:  mtd_accu_add_data
!===============================================================================

subroutine mtd_accu_add_data(cvs, height, widths)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: height
    real(PMFDP)    :: widths(:)
    ! -----------------------------------------------
    integer        :: n, i, gi0
    real(PMFDP)    :: diff, fexparg, fh
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(mtdaccu%PMFAccuType, cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
    else
        insidesamples = insidesamples + 1
        ! increase number of samples
        mtdaccu%nsamples(gi0) = mtdaccu%nsamples(gi0) + 1
    end if

    ! update grid data
    do n=1,mtdaccu%tot_nbins
        fexparg = 0.0d0
        do i=1,NumOfMTDCVs
            diff = MTDCVList(i)%cv%get_deviation(mtdaccu%binpositions(i,n), cvs(i))
            fexparg = fexparg + diff**2 / (2.0d0 * widths(i)**2)
        end do
        fh = height * exp(-fexparg)
        mtdaccu%mtdpotential(n) = mtdaccu%mtdpotential(n) + fh
        do i=1,NumOfMTDCVs
            diff = MTDCVList(i)%cv%get_deviation(mtdaccu%binpositions(i,n), cvs(i))
            mtdaccu%mtdforce(i,n) = mtdaccu%mtdforce(i,n) + fh * diff / widths(i)**2
        end do
    end do

    return

end subroutine mtd_accu_add_data

!===============================================================================
! Subroutine:  mtd_accu_get_data
!===============================================================================

subroutine mtd_accu_get_data(cvs,potential,forces)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: potential
    real(PMFDP)    :: forces(:)
    ! -----------------------------------------------
    integer        :: gi0
    ! --------------------------------------------------------------------------

    potential = 0.0d0
    forces(:) = 0.0d0

    ! get global index to grid for cvs
    gi0 = pmf_accu_globalindex(mtdaccu%PMFAccuType, cvs)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get potential
    potential = mtdaccu%mtdpotential(gi0)

    ! get forces
    forces(:) = mtdaccu%mtdforce(:,gi0)

end subroutine mtd_accu_get_data

!===============================================================================

end module mtd_accu
