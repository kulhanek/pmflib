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

module usabf_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  usabf_accu_init
!===============================================================================

subroutine usabf_accu_init()

    use usabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    usabfaccu%tot_cvs = NumOfUSABFCVs
    allocate(usabfaccu%sizes(usabfaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator!')
    endif

    tot_nbins = 1
    do i=1,usabfaccu%tot_cvs
        usabfaccu%sizes(i)%min_value  = USABFCVList(i)%min_value
        usabfaccu%sizes(i)%max_value  = USABFCVList(i)%max_value
        usabfaccu%sizes(i)%nbins      = USABFCVList(i)%nbins
        usabfaccu%sizes(i)%width      = abs(usabfaccu%sizes(i)%max_value - usabfaccu%sizes(i)%min_value)
        usabfaccu%sizes(i)%bin_width  = usabfaccu%sizes(i)%width / usabfaccu%sizes(i)%nbins
        usabfaccu%sizes(i)%cv         => USABFCVList(i)%cv
        tot_nbins = tot_nbins * usabfaccu%sizes(i)%nbins
    end do

    usabfaccu%tot_nbins = tot_nbins

    ! ABF force arrays
    allocate(  &
            usabfaccu%nsamples(usabfaccu%tot_nbins), &
            usabfaccu%micf(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
            usabfaccu%m2icf(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
            usabfaccu%tvalues(usabfaccu%tot_cvs), &
            usabfaccu%fcs(usabfaccu%tot_cvs), &
            stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (icf)!')
    endif

! enthalpy ---------------------------------------------------------------------
    if( fenthalpy ) then
        allocate(  &
                usabfaccu%mepot(usabfaccu%tot_nbins), &
                usabfaccu%m2epot(usabfaccu%tot_nbins), &
                stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (emthalpy)!')
        endif
    end if

! entropy ----------------------------------------------------------------------
    if( fentropy ) then
        allocate(  &
                usabfaccu%metot(usabfaccu%tot_nbins), &
                usabfaccu%m2etot(usabfaccu%tot_nbins), &
                usabfaccu%c11hh(usabfaccu%tot_cvs,usabfaccu%tot_nbins), &
                stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[US-ABF] Unable to allocate memory for abf accumulator (entropy)!')
        end if
    end if

    call usabf_accu_clear()

    return

end subroutine usabf_accu_init

!===============================================================================
! Subroutine:  usabf_accu_clear
!===============================================================================

subroutine usabf_accu_clear()

    use usabf_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    usabfaccu%nsamples(:)       = 0

    usabfaccu%micf(:,:)         = 0.0d0
    usabfaccu%m2icf(:,:)        = 0.0d0

    if( fenthalpy ) then
        usabfaccu%mepot(:)      = 0.0d0
        usabfaccu%m2epot(:)     = 0.0d0
    end if

    if( fentropy ) then
        usabfaccu%metot(:)      = 0.0d0
        usabfaccu%m2etot(:)     = 0.0d0
        usabfaccu%c11hh(:,:)    = 0.0d0
    end if

end subroutine usabf_accu_clear

!===============================================================================
! Subroutine:  usabf_accu_read
!===============================================================================

subroutine usabf_accu_read(iounit)

    use usabf_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    character(len=PMF_KEYLINE)      :: keyline
    integer                         :: i
    ! --------------------------------------------------------------------------

    do while(.true.)

        ! read keyline
        read(iounit,5,end=500,err=300) keyline

        ! process keyline
        if( pmf_accu_is_header_key(keyline) ) then
            call pmf_accu_read_header(usabfaccu%PMFAccuType,iounit,'US-ABF',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_ibuf_B(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%nsamples)
            ! ------------------------------------
                case('MICF')
                    call pmf_accu_read_rbuf_M(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%micf)
            ! ------------------------------------
                case('M2ICF')
                    call pmf_accu_read_rbuf_M(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%m2icf)
            ! ------------------------------------
                case('TVALUES')
                    call pmf_accu_read_rbuf_C(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%tvalues)
                    do i=1,usabfaccu%tot_cvs
                        USABFCVList(i)%target_value = usabfaccu%tvalues(i)
                    end do
            ! ------------------------------------
                case('FCS')
                    call pmf_accu_read_rbuf_C(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%fcs)
                    do i=1,usabfaccu%tot_cvs
                        USABFCVList(i)%force_constant = usabfaccu%fcs(i)
                    end do
            ! ------------------------------------
                case('MEPOT')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%mepot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,USABF_OUT)
                    end if
            ! ------------------------------------
                case('M2EPOT')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%m2epot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,USABF_OUT)
                    end if
            ! ------------------------------------
                case('METOT')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_B(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%metot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,USABF_OUT)
                    end if
            ! ------------------------------------
                case('M2ETOT')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_B(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%m2etot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,USABF_OUT)
                    end if
           ! ------------------------------------
                case('C11HH')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_M(usabfaccu%PMFAccuType,iounit,keyline,usabfaccu%c11hh)
                    else
                        call pmf_accu_skip_section(iounit,keyline,USABF_OUT)
                    end if
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,USABF_OUT)
            end select
        end if
    end do

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[US-ABF] Unable to read from the accumulator - keyline!')

end subroutine usabf_accu_read

!===============================================================================
! Subroutine:  usabf_accu_write
!===============================================================================

subroutine usabf_accu_write(iounit)

    use usabf_dat

    implicit none
    integer     :: iounit
    ! --------------------------------------------
    integer     :: i
    !---------------------------------------------------------------------------

    do i=1,usabfaccu%tot_cvs
        usabfaccu%tvalues(i)  = USABFCVList(i)%target_value
        usabfaccu%fcs(i)      = USABFCVList(i)%force_constant
    end do

    usabfaccu%method = 'US-ABF'
    call pmf_accu_write_header(usabfaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(usabfaccu%PMFAccuType,iounit,'NSAMPLES',     'AD',usabfaccu%nsamples)
    call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'MICF',         'WA',usabfaccu%micf)
    call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'M2ICF',        'M2',usabfaccu%m2icf,'MICF')
    call pmf_accu_write_rbuf_C(usabfaccu%PMFAccuType,iounit,'TVALUES',      'IG',usabfaccu%tvalues)
    call pmf_accu_write_rbuf_C(usabfaccu%PMFAccuType,iounit,'FCS',          'IG',usabfaccu%fcs)

    if( fenthalpy ) then
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'MEPOT',    'WA',usabfaccu%mepot)
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'M2EPOT',   'M2',usabfaccu%m2epot,'MEPOT')
    end if

    if( fentropy ) then
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'METOT',    'WA',usabfaccu%metot)
        call pmf_accu_write_rbuf_B(usabfaccu%PMFAccuType,iounit,'M2ETOT',   'M2',usabfaccu%m2etot,'METOT')
        call pmf_accu_write_rbuf_M(usabfaccu%PMFAccuType,iounit,'C11HH',    'CO',usabfaccu%c11hh,'MICF','METOT')
    end if

end subroutine usabf_accu_write

!===============================================================================
! Subroutine:  usabf_accu_add_data_online
!===============================================================================

subroutine usabf_accu_add_data_online(cvs,gfx,epot,etot)

    use usabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: etot
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, icf
    real(PMFDP)    :: detot1, detot2
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dicf1, dicf2
    ! --------------------------------------------------------------------------

    ! reset the accumulated data if requested
    if( (faccurst .ge. 0) .and. (faccurst .eq. fstep) ) then
        write(USABF_OUT,10) fstep
        call usabf_accu_clear
    end if

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(usabfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    usabfaccu%nsamples(gi0) = usabfaccu%nsamples(gi0) + 1
    invn = 1.0d0 / real(usabfaccu%nsamples(gi0),PMFDP)

    if( fenthalpy ) then
        ! potential energy
        depot1 = epot - usabfaccu%mepot(gi0)
        usabfaccu%mepot(gi0)  = usabfaccu%mepot(gi0)  + depot1 * invn
        depot2 = epot - usabfaccu%mepot(gi0)
        usabfaccu%m2epot(gi0) = usabfaccu%m2epot(gi0) + depot1 * depot2
    end if

    if( fentropy ) then
        ! total energy
        detot1 = etot - usabfaccu%metot(gi0)
        usabfaccu%metot(gi0)  = usabfaccu%metot(gi0)  + detot1 * invn
        detot2 = etot - usabfaccu%metot(gi0)
        usabfaccu%m2etot(gi0) = usabfaccu%m2etot(gi0) + detot1 * detot2
    end if

    do i=1,NumOfUSABFCVs
        icf = gfx(i)

        ! write(12458,*) icf

        dicf1 = - icf - usabfaccu%micf(i,gi0)
        usabfaccu%micf(i,gi0)  = usabfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = - icf -  usabfaccu%micf(i,gi0)
        usabfaccu%m2icf(i,gi0) = usabfaccu%m2icf(i,gi0) + dicf1 * dicf2

        if( fentropy ) then
            usabfaccu%c11hh(i,gi0)  = usabfaccu%c11hh(i,gi0) + dicf1 * detot2
        end if
    end do

 10 format('# Resetting the accumulator at: ', I9)

end subroutine usabf_accu_add_data_online

!===============================================================================

end module usabf_accu

