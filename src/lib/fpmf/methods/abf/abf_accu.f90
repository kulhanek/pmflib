!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2022-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    allocate(   abfaccu%weights(abfaccu%tot_nbins),                 &
                abfaccu%binpos(abfaccu%tot_cvs,abfaccu%tot_nbins),  &
                abfaccu%nsamples(abfaccu%tot_nbins),                &
                abfaccu%micf(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                abfaccu%m2icf(abfaccu%tot_cvs,abfaccu%tot_nbins),   &
                abfaccu%mgfx(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                abfaccu%m2gfx(abfaccu%tot_cvs,abfaccu%tot_nbins),   &
                abfaccu%bnsamples(abfaccu%tot_nbins),               &
                abfaccu%bmicf(abfaccu%tot_cvs,abfaccu%tot_nbins),   &
                stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (abfforce)!')
    endif

    if( fenthalpy ) then
        allocate(   abfaccu%meint(abfaccu%tot_nbins),   &
                    abfaccu%m2eint(abfaccu%tot_nbins),  &
                    abfaccu%mepot(abfaccu%tot_nbins),   &
                    abfaccu%m2epot(abfaccu%tot_nbins),  &
                    abfaccu%merst(abfaccu%tot_nbins),   &
                    abfaccu%m2erst(abfaccu%tot_nbins),  &
                    abfaccu%mekin(abfaccu%tot_nbins),   &
                    abfaccu%m2ekin(abfaccu%tot_nbins),  &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (enthalpy)!')
        endif
    end if

    if( fentropy ) then
        allocate(   abfaccu%metot(abfaccu%tot_nbins),                   &
                    abfaccu%m2etot(abfaccu%tot_nbins),                  &
                    abfaccu%mpp(abfaccu%tot_cvs,abfaccu%tot_nbins),     &
                    abfaccu%m2pp(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%mpn(abfaccu%tot_cvs,abfaccu%tot_nbins),     &
                    abfaccu%m2pn(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (entropy)!')
        endif
    end if

    if( fentropy .and. fentdecomp ) then
        allocate(   abfaccu%ntds(abfaccu%tot_nbins),                        &
                    abfaccu%mtdsepot(abfaccu%tot_nbins),                    &
                    abfaccu%m2tdsepot(abfaccu%tot_nbins),                   &
                    abfaccu%mtdserst(abfaccu%tot_nbins),                    &
                    abfaccu%m2tdserst(abfaccu%tot_nbins),                   &
                    abfaccu%mtdsekin(abfaccu%tot_nbins),                    &
                    abfaccu%m2tdsekin(abfaccu%tot_nbins),                   &
                    abfaccu%mtdsfx(abfaccu%tot_cvs,abfaccu%tot_nbins),      &
                    abfaccu%m2tdsfx(abfaccu%tot_cvs,abfaccu%tot_nbins),     &
                    abfaccu%mtdssx(abfaccu%tot_cvs,abfaccu%tot_nbins),      &
                    abfaccu%m2tdssx(abfaccu%tot_cvs,abfaccu%tot_nbins),     &
                    abfaccu%mtdsvx(abfaccu%tot_cvs,abfaccu%tot_nbins),      &
                    abfaccu%m2tdsvx(abfaccu%tot_cvs,abfaccu%tot_nbins),     &
                    abfaccu%mtdsbx(abfaccu%tot_cvs,abfaccu%tot_nbins),      &
                    abfaccu%m2tdsbx(abfaccu%tot_cvs,abfaccu%tot_nbins),     &
                    abfaccu%c11tdsfp(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsfk(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsfr(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdssp(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdssk(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdssr(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsvp(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsvk(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsvr(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsbp(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsbk(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%c11tdsbr(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%ntds_h(abfaccu%tot_nbins),                      &
                    abfaccu%mtdsfx_h(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%m2tdsfx_h(abfaccu%tot_cvs,abfaccu%tot_nbins),   &
                    abfaccu%mtdsekin_h(abfaccu%tot_nbins),                  &
                    abfaccu%m2tdsekin_h(abfaccu%tot_nbins),                 &
                    abfaccu%c11tdsfk_h(abfaccu%tot_cvs,abfaccu%tot_nbins),  &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (entropy_decompose)!')
        endif
    end if

    if( frecord ) then
        allocate(   abfaccu%tcvs(abfaccu%tot_cvs,fnstlim),      &
                    abfaccu%tbicf(abfaccu%tot_cvs,fnstlim),     &
                    abfaccu%tepot(fnstlim),                     &
                    abfaccu%terst(fnstlim),                     &
                    abfaccu%tekin(fnstlim),                     &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (record)!')
        endif
    end if

    if( fserver_enabled ) then
        ! ABF incremental force arrays - allocate all to simplify F90/C++ interface
        allocate(   abfaccu%inc_nsamples(abfaccu%tot_nbins),                &
                    abfaccu%inc_micf(abfaccu%tot_cvs,abfaccu%tot_nbins),    &
                    abfaccu%inc_m2icf(abfaccu%tot_cvs,abfaccu%tot_nbins),   &
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
    abfaccu%nsamples(:)     = 0.0d0
    abfaccu%micf(:,:)       = 0.0d0
    abfaccu%m2icf(:,:)      = 0.0d0
    abfaccu%mgfx(:,:)       = 0.0d0
    abfaccu%m2gfx(:,:)      = 0.0d0

    abfaccu%bnsamples(:)    = 0
    abfaccu%bmicf(:,:)      = 0.0d0

    if( fenthalpy ) then
        abfaccu%meint(:)        = 0.0d0
        abfaccu%m2eint(:)       = 0.0d0
        abfaccu%mepot(:)        = 0.0d0
        abfaccu%m2epot(:)       = 0.0d0
        abfaccu%merst(:)        = 0.0d0
        abfaccu%m2erst(:)       = 0.0d0
        abfaccu%mekin(:)        = 0.0d0
        abfaccu%m2ekin(:)       = 0.0d0
    end if

    if( fentropy ) then
        abfaccu%metot(:)        = 0.0d0
        abfaccu%m2etot(:)       = 0.0d0
        abfaccu%mpp(:,:)        = 0.0d0
        abfaccu%m2pp(:,:)       = 0.0d0
        abfaccu%mpn(:,:)        = 0.0d0
        abfaccu%m2pn(:,:)       = 0.0d0
    end if

    if( fentropy .and. fentdecomp ) then
        abfaccu%ntds(:)         = 0.0d0
        abfaccu%mtdsepot(:)     = 0.0d0
        abfaccu%m2tdsepot(:)    = 0.0d0
        abfaccu%mtdserst(:)     = 0.0d0
        abfaccu%m2tdserst(:)    = 0.0d0
        abfaccu%mtdsekin(:)     = 0.0d0
        abfaccu%m2tdsekin(:)    = 0.0d0
        abfaccu%mtdsfx(:,:)     = 0.0d0
        abfaccu%m2tdsfx(:,:)    = 0.0d0
        abfaccu%mtdssx(:,:)     = 0.0d0
        abfaccu%m2tdssx(:,:)    = 0.0d0
        abfaccu%mtdsvx(:,:)     = 0.0d0
        abfaccu%m2tdsvx(:,:)    = 0.0d0
        abfaccu%mtdsbx(:,:)     = 0.0d0
        abfaccu%m2tdsbx(:,:)    = 0.0d0
        abfaccu%c11tdsfp(:,:)   = 0.0d0
        abfaccu%c11tdsfr(:,:)   = 0.0d0
        abfaccu%c11tdsfk(:,:)   = 0.0d0
        abfaccu%c11tdssp(:,:)   = 0.0d0
        abfaccu%c11tdssr(:,:)   = 0.0d0
        abfaccu%c11tdssk(:,:)   = 0.0d0
        abfaccu%c11tdsvp(:,:)   = 0.0d0
        abfaccu%c11tdsvr(:,:)   = 0.0d0
        abfaccu%c11tdsvk(:,:)   = 0.0d0
        abfaccu%c11tdsbp(:,:)   = 0.0d0
        abfaccu%c11tdsbr(:,:)   = 0.0d0
        abfaccu%c11tdsbk(:,:)   = 0.0d0

        abfaccu%ntds_h(:)       = 0.0d0
        abfaccu%mtdsfx_h(:,:)   = 0.0d0
        abfaccu%m2tdsfx_h(:,:)  = 0.0d0
        abfaccu%mtdsekin_h(:)   = 0.0d0
        abfaccu%m2tdsekin_h(:)  = 0.0d0
        abfaccu%c11tdsfk_h(:,:) = 0.0d0
    end if

    if( frecord ) then
        abfaccu%tcvs(:,:)       = 0.0d0
        abfaccu%tbicf(:,:)      = 0.0d0
        abfaccu%tepot(:)        = 0.0d0
        abfaccu%terst(:)        = 0.0d0
        abfaccu%tekin(:)        = 0.0d0
    end if

    if( fserver_enabled ) then
        abfaccu%inc_nsamples(:) = 0
        abfaccu%inc_micf(:,:)   = 0.0d0
        abfaccu%inc_m2icf(:,:)  = 0.0d0
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
                case('MEINT')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%meint)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2EINT')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2eint)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
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
                case('MERST')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%merst)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2ERST')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2erst)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('MEKIN')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%mekin)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2EKIN')
                    if( fenthalpy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2ekin)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('METOT')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%metot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2ETOT')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_B(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2etot)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('MPP')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_M(abfaccu%PMFAccuType,iounit,keyline,abfaccu%mpp)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2PP')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_M(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2pp)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
           ! ------------------------------------
                case('MPN')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_M(abfaccu%PMFAccuType,iounit,keyline,abfaccu%mpn)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case('M2PN')
                    if( fentropy ) then
                        call pmf_accu_read_rbuf_M(abfaccu%PMFAccuType,iounit,keyline,abfaccu%m2pn)
                    else
                        call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
                    end if
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,ABF_OUT)
            end select
        end if
    end do

    ! copy MICF and NSAMPLES into
    abfaccu%bnsamples(:)    = abfaccu%nsamples(:)
    abfaccu%bmicf(:,:)      = abfaccu%micf(:,:)

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
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MGFX',       'WA',abfaccu%mgfx,  'NSAMPLES')
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2GFX',      'M2',abfaccu%m2gfx, 'NSAMPLES','MGFX')

    if( fenthalpy ) then
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MEINT',  'WA',abfaccu%meint, 'NSAMPLES')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2EINT', 'M2',abfaccu%m2eint,'NSAMPLES','MEINT')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MEPOT',  'WA',abfaccu%mepot, 'NSAMPLES')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2EPOT', 'M2',abfaccu%m2epot,'NSAMPLES','MEPOT')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MERST',  'WA',abfaccu%merst, 'NSAMPLES')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2ERST', 'M2',abfaccu%m2erst,'NSAMPLES','MERST')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MEKIN',  'WA',abfaccu%mekin, 'NSAMPLES')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2EKIN', 'M2',abfaccu%m2ekin,'NSAMPLES','MEKIN')
    end if

    if( fentropy ) then
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'METOT',      'WA',abfaccu%metot,     'NSAMPLES')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2ETOT',     'M2',abfaccu%m2etot,    'NSAMPLES','METOT')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MPP',        'WA',abfaccu%mpp,       'NSAMPLES')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2PP',       'M2',abfaccu%m2pp,      'NSAMPLES','MPP')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MPN',        'WA',abfaccu%mpn,       'NSAMPLES')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2PN',       'M2',abfaccu%m2pn,      'NSAMPLES','MPN')
    end if

    if( fentropy .and. fentdecomp ) then
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'NTDS',       'AD',abfaccu%ntds)
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MTDSEPOT',   'WA',abfaccu%mtdsepot,  'NTDS')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2TDSEPOT',  'M2',abfaccu%m2tdsepot, 'NTDS','MTDSEPOT')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MTDSERST',   'WA',abfaccu%mtdserst,  'NTDS')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2TDSERST',  'M2',abfaccu%m2tdserst, 'NTDS','MTDSERST')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MTDSEKIN',   'WA',abfaccu%mtdsekin,  'NTDS')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2TDSEKIN',  'M2',abfaccu%m2tdsekin, 'NTDS','MTDSEKIN')

        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MTDSFX',     'WA',abfaccu%mtdsfx,    'NTDS')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2TDSFX',    'M2',abfaccu%m2tdsfx,   'NTDS','MTDSFX')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MTDSSX',     'WA',abfaccu%mtdssx,    'NTDS')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2TDSSX',    'M2',abfaccu%m2tdssx,   'NTDS','MTDSSX')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MTDSVX',     'WA',abfaccu%mtdsvx,    'NTDS')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2TDSVX',    'M2',abfaccu%m2tdsvx,   'NTDS','MTDSVX')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MTDSBX',     'WA',abfaccu%mtdsbx,    'NTDS')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2TDSBX',    'M2',abfaccu%m2tdsbx,   'NTDS','MTDSBX')

        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSFP',   'CO',abfaccu%c11tdsfp,  'NTDS','MTDSFX','MTDSEPOT')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSFR',   'CO',abfaccu%c11tdsfr,  'NTDS','MTDSFX','MTDSERST')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSFK',   'CO',abfaccu%c11tdsfk,  'NTDS','MTDSFX','MTDSEKIN')

        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSSP',   'CO',abfaccu%c11tdssp,  'NTDS','MTDSSX','MTDSEPOT')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSSR',   'CO',abfaccu%c11tdssr,  'NTDS','MTDSSX','MTDSERST')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSSK',   'CO',abfaccu%c11tdssk,  'NTDS','MTDSSX','MTDSEKIN')

        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSVP',   'CO',abfaccu%c11tdsvp,  'NTDS','MTDSVX','MTDSEPOT')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSVR',   'CO',abfaccu%c11tdsvr,  'NTDS','MTDSVX','MTDSERST')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSVK',   'CO',abfaccu%c11tdsvk,  'NTDS','MTDSVX','MTDSEKIN')

        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSBP',   'CO',abfaccu%c11tdsbp,  'NTDS','MTDSBX','MTDSEPOT')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSBR',   'CO',abfaccu%c11tdsbr,  'NTDS','MTDSBX','MTDSERST')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSBK',   'CO',abfaccu%c11tdsbk,  'NTDS','MTDSBX','MTDSEKIN')

        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'NTDS_H',      'AD',abfaccu%ntds_h)
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'MTDSEKIN_H',  'WA',abfaccu%mtdsekin_h,  'NTDS_H')
        call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'M2TDSEKIN_H', 'M2',abfaccu%m2tdsekin_h, 'NTDS_H','MTDSEKIN_H')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MTDSFX_H',    'WA',abfaccu%mtdsfx_h,    'NTDS_H')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2TDSFX_H',   'M2',abfaccu%m2tdsfx_h,   'NTDS_H','MTDSFX_H')
        call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'C11TDSFK_H',  'CO',abfaccu%c11tdsfk_h,  &
                                                                             'NTDS_H','MTDSFX_H','MTDSEKIN_H')
    end if

    call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'BNSAMPLES',  'IG',abfaccu%bnsamples)
    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'BMICF',      'IG',abfaccu%bmicf, 'BNSAMPLES')

    if( full ) then
        if( frecord ) then
            call pmf_accu_write_rbuf_TC(abfaccu%PMFAccuType,iounit,     'TCVS', abfaccu%tcvs)
            call pmf_accu_write_rbuf_TC(abfaccu%PMFAccuType,iounit,     'TBICF',abfaccu%tbicf)
            call pmf_accu_write_rbuf_T(iounit,                          'TEPOT',abfaccu%tepot)
            call pmf_accu_write_rbuf_T(iounit,                          'TERST',abfaccu%terst)
            call pmf_accu_write_rbuf_T(iounit,                          'TEKIN',abfaccu%tekin)
        end if
    end if

end subroutine abf_accu_write

!===============================================================================
! Subroutine:  abf_accu_add_data_online
!===============================================================================

subroutine abf_accu_add_data_online(cvs,gfx,bfx,epot,erst,ekin,etot)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP),intent(in)  :: cvs(:)
    real(PMFDP),intent(in)  :: gfx(:)
    real(PMFDP),intent(in)  :: bfx(:)
    real(PMFDP),intent(in)  :: epot
    real(PMFDP),intent(in)  :: erst
    real(PMFDP),intent(in)  :: ekin
    real(PMFDP),intent(in)  :: etot
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, icf, igf, eint
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: derst1, derst2
    real(PMFDP)    :: dekin1, dekin2
    real(PMFDP)    :: detot1, detot2
    real(PMFDP)    :: deint1, deint2
    real(PMFDP)    :: dicf1, dicf2
    real(PMFDP)    :: dgfx1, dgfx2
    real(PMFDP)    :: dpp, dpp1, dpp2
    real(PMFDP)    :: dpn, dpn1, dpn2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! get total biasing ICF
    picf(:) = - (gfx(:) - bfx(:))

    ! increase number of samples
    abfaccu%nsamples(gi0) = abfaccu%nsamples(gi0) + 1.0d0
    invn = 1.0d0 / abfaccu%nsamples(gi0)

    if( fenthalpy ) then
        ! internal energy
        eint = epot + erst
        deint1 = eint - abfaccu%meint(gi0)
        abfaccu%meint(gi0)  = abfaccu%meint(gi0)  + deint1 * invn
        deint2 = eint - abfaccu%meint(gi0)
        abfaccu%m2eint(gi0) = abfaccu%m2eint(gi0) + deint1 * deint2

        ! potential energy
        depot1 = epot - abfaccu%mepot(gi0)
        abfaccu%mepot(gi0)  = abfaccu%mepot(gi0)  + depot1 * invn
        depot2 = epot - abfaccu%mepot(gi0)
        abfaccu%m2epot(gi0) = abfaccu%m2epot(gi0) + depot1 * depot2

        ! restraint energy
        derst1 = erst - abfaccu%merst(gi0)
        abfaccu%merst(gi0)  = abfaccu%merst(gi0)  + derst1 * invn
        derst2 = erst - abfaccu%merst(gi0)
        abfaccu%m2erst(gi0) = abfaccu%m2erst(gi0) + derst1 * derst2

        ! kinetic energy
        dekin1 = ekin - abfaccu%mekin(gi0)
        abfaccu%mekin(gi0)  = abfaccu%mekin(gi0)  + dekin1 * invn
        dekin2 = ekin - abfaccu%mekin(gi0)
        abfaccu%m2ekin(gi0) = abfaccu%m2ekin(gi0) + dekin1 * dekin2
    end if

    if( fentropy ) then
        ! total energy
        detot1 = etot - abfaccu%metot(gi0)
        abfaccu%metot(gi0)  = abfaccu%metot(gi0)  + detot1 * invn
        detot2 = etot - abfaccu%metot(gi0)
        abfaccu%m2etot(gi0) = abfaccu%m2etot(gi0) + detot1 * detot2
    end if

    do i=1,NumOfABFCVs
        icf = picf(i)
        dicf1 = icf - abfaccu%micf(i,gi0)
        abfaccu%micf(i,gi0)  = abfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = icf - abfaccu%micf(i,gi0)
        abfaccu%m2icf(i,gi0) = abfaccu%m2icf(i,gi0) + dicf1 * dicf2

        igf = -gfx(i)
        dgfx1 = igf - abfaccu%mgfx(i,gi0)
        abfaccu%mgfx(i,gi0)  = abfaccu%mgfx(i,gi0)  + dgfx1 * invn
        dgfx2 = igf - abfaccu%mgfx(i,gi0)
        abfaccu%m2gfx(i,gi0) = abfaccu%m2gfx(i,gi0) + dgfx1 * dgfx2

        if( fentropy ) then
            if( ftds_add_bias ) then
                icf = picf(i)
            else
                icf = - gfx(i)
            end if

            dpp = icf + etot
            dpp1 = dpp - abfaccu%mpp(i,gi0)
            abfaccu%mpp(i,gi0)  = abfaccu%mpp(i,gi0)  + dpp1 * invn
            dpp2 = dpp - abfaccu%mpp(i,gi0)
            abfaccu%m2pp(i,gi0) = abfaccu%m2pp(i,gi0) + dpp1 * dpp2

            dpn = icf - etot
            dpn1 = dpn - abfaccu%mpn(i,gi0)
            abfaccu%mpn(i,gi0)  = abfaccu%mpn(i,gi0)  + dpn1 * invn
            dpn2 = dpn - abfaccu%mpn(i,gi0)
            abfaccu%m2pn(i,gi0) = abfaccu%m2pn(i,gi0) + dpn1 * dpn2
        end if
    end do

    if( fupdate_abf ) then
        if( fserver_enabled ) then
            abfaccu%inc_nsamples(gi0) = abfaccu%inc_nsamples(gi0) + 1.0d0
            invn = 1.0d0 / abfaccu%inc_nsamples(gi0)

            do i=1,NumOfABFCVs
                icf = picf(i)
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
            icf = picf(i)
            dicf1 = icf - abfaccu%bmicf(i,gi0)
            abfaccu%bmicf(i,gi0)  = abfaccu%bmicf(i,gi0)  + dicf1 * invn
        end do
    end if

end subroutine abf_accu_add_data_online

!===============================================================================
! Subroutine:  abf_accu_add_data_entropy_decompose
!===============================================================================

subroutine abf_accu_add_data_entropy_decompose(cvs,fx,sx,vx,bx,epot,erst,ekin)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP),intent(in)  :: cvs(:)
    real(PMFDP),intent(in)  :: fx(:)
    real(PMFDP),intent(in)  :: sx(:)
    real(PMFDP),intent(in)  :: vx(:)
    real(PMFDP),intent(in)  :: bx(:)
    real(PMFDP),intent(in)  :: epot
    real(PMFDP),intent(in)  :: erst
    real(PMFDP),intent(in)  :: ekin
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, ifx, ivx, ibx, isx
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: derst1, derst2
    real(PMFDP)    :: dekin1, dekin2
    real(PMFDP)    :: difx1, difx2
    real(PMFDP)    :: divx1, divx2
    real(PMFDP)    :: dibx1, dibx2
    real(PMFDP)    :: disx1, disx2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out-of-valid area

    ! increase number of samples
    abfaccu%ntds(gi0) = abfaccu%ntds(gi0) + 1.0d0
    invn = 1.0d0 / abfaccu%ntds(gi0)

    depot1 = epot - abfaccu%mtdsepot(gi0)
    abfaccu%mtdsepot(gi0)  = abfaccu%mtdsepot(gi0)  + depot1 * invn
    depot2 = epot - abfaccu%mtdsepot(gi0)
    abfaccu%m2tdsepot(gi0) = abfaccu%m2tdsepot(gi0) + depot1 * depot2

    ! restraint energy
    derst1 = erst - abfaccu%mtdserst(gi0)
    abfaccu%mtdserst(gi0)  = abfaccu%mtdserst(gi0)  + derst1 * invn
    derst2 = erst - abfaccu%mtdserst(gi0)
    abfaccu%m2tdserst(gi0) = abfaccu%m2tdserst(gi0) + derst1 * derst2

    ! kinetic energy
    dekin1 = ekin - abfaccu%mtdsekin(gi0)
    abfaccu%mtdsekin(gi0)  = abfaccu%mtdsekin(gi0)  + dekin1 * invn
    dekin2 = ekin - abfaccu%mtdsekin(gi0)
    abfaccu%m2tdsekin(gi0) = abfaccu%m2tdsekin(gi0) + dekin1 * dekin2

    do i=1,NumOfABFCVs
        ifx = - fx(i)
        difx1 = ifx - abfaccu%mtdsfx(i,gi0)
        abfaccu%mtdsfx(i,gi0)  = abfaccu%mtdsfx(i,gi0)  + difx1 * invn
        difx2 = ifx - abfaccu%mtdsfx(i,gi0)
        abfaccu%m2tdsfx(i,gi0) = abfaccu%m2tdsfx(i,gi0) + difx1 * difx2

        isx = - sx(i)
        disx1 = isx - abfaccu%mtdssx(i,gi0)
        abfaccu%mtdssx(i,gi0)  = abfaccu%mtdssx(i,gi0)  + disx1 * invn
        disx2 = isx - abfaccu%mtdssx(i,gi0)
        abfaccu%m2tdssx(i,gi0) = abfaccu%m2tdssx(i,gi0) + disx1 * disx2

        ivx = - vx(i)
        divx1 = ivx - abfaccu%mtdsvx(i,gi0)
        abfaccu%mtdsvx(i,gi0)  = abfaccu%mtdsvx(i,gi0)  + divx1 * invn
        divx2 = ivx - abfaccu%mtdsvx(i,gi0)
        abfaccu%m2tdsvx(i,gi0) = abfaccu%m2tdsvx(i,gi0) + divx1 * divx2

        ibx =   bx(i)
        dibx1 = ibx - abfaccu%mtdsbx(i,gi0)
        abfaccu%mtdsbx(i,gi0)  = abfaccu%mtdsbx(i,gi0)  + dibx1 * invn
        dibx2 = ibx - abfaccu%mtdsbx(i,gi0)
        abfaccu%m2tdsbx(i,gi0) = abfaccu%m2tdsbx(i,gi0) + dibx1 * dibx2

        abfaccu%c11tdsfp(i,gi0)  = abfaccu%c11tdsfp(i,gi0) + difx1 * depot2
        abfaccu%c11tdsfr(i,gi0)  = abfaccu%c11tdsfr(i,gi0) + difx1 * derst2
        abfaccu%c11tdsfk(i,gi0)  = abfaccu%c11tdsfk(i,gi0) + difx1 * dekin2

        abfaccu%c11tdssp(i,gi0)  = abfaccu%c11tdssp(i,gi0) + disx1 * depot2
        abfaccu%c11tdssr(i,gi0)  = abfaccu%c11tdssr(i,gi0) + disx1 * derst2
        abfaccu%c11tdssk(i,gi0)  = abfaccu%c11tdssk(i,gi0) + disx1 * dekin2

        abfaccu%c11tdsvp(i,gi0)  = abfaccu%c11tdsvp(i,gi0) + divx1 * depot2
        abfaccu%c11tdsvr(i,gi0)  = abfaccu%c11tdsvr(i,gi0) + divx1 * derst2
        abfaccu%c11tdsvk(i,gi0)  = abfaccu%c11tdsvk(i,gi0) + divx1 * dekin2

        abfaccu%c11tdsbp(i,gi0)  = abfaccu%c11tdsbp(i,gi0) + dibx1 * depot2
        abfaccu%c11tdsbr(i,gi0)  = abfaccu%c11tdsbr(i,gi0) + dibx1 * derst2
        abfaccu%c11tdsbk(i,gi0)  = abfaccu%c11tdsbk(i,gi0) + dibx1 * dekin2
    end do

end subroutine abf_accu_add_data_entropy_decompose

!===============================================================================
! Subroutine:  abf_accu_add_data_entropy_decompose_hs
!===============================================================================

subroutine abf_accu_add_data_entropy_decompose_hs(cvs,fx,ekin)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP),intent(in)  :: cvs(:)
    real(PMFDP),intent(in)  :: fx(:)
    real(PMFDP),intent(in)  :: ekin
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, ifx
    real(PMFDP)    :: dekin1, dekin2
    real(PMFDP)    :: difx1, difx2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for cvs values
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out-of-valid area

    ! increase number of samples
    abfaccu%ntds_h(gi0) = abfaccu%ntds_h(gi0) + 1.0d0
    invn = 1.0d0 / abfaccu%ntds_h(gi0)

    ! kinetic energy
    dekin1 = ekin - abfaccu%mtdsekin_h(gi0)
    abfaccu%mtdsekin_h(gi0)  = abfaccu%mtdsekin_h(gi0)  + dekin1 * invn
    dekin2 = ekin - abfaccu%mtdsekin_h(gi0)
    abfaccu%m2tdsekin_h(gi0) = abfaccu%m2tdsekin_h(gi0) + dekin1 * dekin2

    do i=1,NumOfABFCVs
        ifx = - fx(i)
        difx1 = ifx - abfaccu%mtdsfx_h(i,gi0)
        abfaccu%mtdsfx_h(i,gi0)  = abfaccu%mtdsfx_h(i,gi0)  + difx1 * invn
        difx2 = ifx - abfaccu%mtdsfx_h(i,gi0)
        abfaccu%m2tdsfx_h(i,gi0) = abfaccu%m2tdsfx_h(i,gi0) + difx1 * difx2

        abfaccu%c11tdsfk_h(i,gi0)  = abfaccu%c11tdsfk_h(i,gi0) + difx1 * dekin2
    end do

end subroutine abf_accu_add_data_entropy_decompose_hs

!===============================================================================
! Subroutine:  abf_accu_add_data_record_lf
!===============================================================================

subroutine abf_accu_add_data_record_lf(cvs,bicf,epot,erst,ekin)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: bicf(:)       ! applied bias
    real(PMFDP)    :: epot
    real(PMFDP)    :: erst
    real(PMFDP)    :: ekin
    ! --------------------------------------------------------------------------

    abfaccu%tcvs(:,fstep)       = cvs(:)
    abfaccu%tbicf(:,fstep)      = bicf(:)
    abfaccu%tepot(fstep)        = epot
    abfaccu%terst(fstep)        = erst
    if( fstep .gt. 1 ) then
        abfaccu%tekin(fstep-1)  = ekin
    end if

end subroutine abf_accu_add_data_record_lf

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
        dx = ABFCVList(i)%cv%get_deviation(cvs1(i),cvs2(i)) / (ABFCVList(i)%wfac*abfaccu%PMFAccuType%sizes(i)%bin_width)
        u2 = u2 + dx**2
    end do

! get value of kernel
! https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
    kval = 0.0d0
    select case(fsmooth_kernel)
    ! Epanechnikov (parabolic)
        case(0)
            if( u2 .lt.  1.0d0 ) then
                kval = 3.0d0/4.0d0*(1.0d0-u2)
            end if
     ! Triweight
        case(1)
            if( u2 .lt.  1.0d0 ) then
                kval = 35.0d0/32.0d0*(1.0d0-u2)**3
            end if
        case default
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Not implemented smoothing kernel in abf_get_skernel!')
    end select

end function abf_get_skernel

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

    if( fswitch2zero ) then
        call abf_get_switching_factors(cvs)
    end if

    w      = abfaccu%weights(gi0)
    gfx(:) = w * abfaccu%bmicf(:,gi0) * sfac(:)

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

    if( fswitch2zero ) then
        call abf_get_switching_factors(cvs)
    end if

    ! get number of samples
    n = abfaccu%bnsamples(gi0)
    if( n .gt. 0 ) then
        sc_ramp = 1.0d0
        if( n .le. fhramp_max ) then
            sc_ramp = 0.0d0
            if( n .gt. fhramp_min ) then
                sc_ramp = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        w      = abfaccu%weights(gi0)
        gfx(:) = w * sc_ramp * abfaccu%bmicf(:,gi0) * sfac(:)
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
    real(PMFDP)    :: stot_weight,w,kw,rw,n
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out of valid area

    if( fswitch2zero ) then
        call abf_get_switching_factors(cvs)
    end if

! first calculate kernel values
    sweights(:) = 0.0d0
    stot_weight = 0.0d0
    do ni=1,max_snb_size
        si0 =  snb_list(ni,gi0)
        if( si0 .le. 0 ) cycle
        w = abf_get_skernel(abfaccu%binpos(:,si0),cvs(:))
        stot_weight  = stot_weight + w
        sweights(ni) = w
    end do
    if( stot_weight .eq. 0.0d0 ) return ! out of valid area

   ! write(4781,*) sweights

! normalize weights
    do ni=1,max_snb_size
        if( snb_list(ni,gi0) .le. 0 ) cycle
        sweights(ni) = sweights(ni) / stot_weight
    end do

! get smoothed mean forces
    do ni=1,max_snb_size
        si0 = snb_list(ni,gi0)
        if( si0 .le. 0 ) cycle
        kw = sweights(ni)
        rw = 1.0d0
        n = abfaccu%bnsamples(si0)
        if( n .le. fhramp_max ) then
            rw = 0.0d0
            if( n .gt. fhramp_min ) then
                rw = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        w = rw * kw * abfaccu%weights(si0)
        gfx(:) = gfx(:) + w * abfaccu%bmicf(:,si0)
    end do

    ! apply switching factors
    gfx(:) = gfx(:) * sfac(:)

end subroutine abf_accu_get_data_ksmooth

!===============================================================================
! Subroutine:  abf_accu_get_data_lsmooth
!===============================================================================

subroutine abf_accu_get_data_lsmooth(cvs,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils
    use pmf_accu

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0,gi1
    real(PMFDP)    :: f0,f1,dx,rw,n,c0,c1
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(abfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) return ! out of valid area

    if( fswitch2zero ) then
        call abf_get_switching_factors(cvs)
    end if

    rw = 1.0d0
    n = abfaccu%bnsamples(gi0)
    if( n .le. fhramp_max ) then
        rw = 0.0d0
        if( n .gt. fhramp_min ) then
            rw = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
        end if
    end if
    f0 = rw * abfaccu%bmicf(1,gi0)

    dx = cvs(1) - abfaccu%binpos(1,gi0)
    if( dx .gt. 0 ) then
        gi1 = gi0 + 1
        c1 = dx / abfaccu%PMFAccuType%sizes(1)%bin_width
        c0 = 1.0d0 - c1
    else
        gi1 = gi0 - 1
        c0 = 1.0d0 + dx / abfaccu%PMFAccuType%sizes(1)%bin_width
        c1 = 1.0d0 - c0
    end if

    if( (gi1 .gt. 1) .and. (gi1 .le. abfaccu%PMFAccuType%tot_nbins) ) then
        rw = 1.0d0
        n = abfaccu%bnsamples(gi1)
        if( n .le. fhramp_max ) then
            rw = 0.0d0
            if( n .gt. fhramp_min ) then
                rw = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        f1 = rw * abfaccu%bmicf(1,gi1)
    else
        f1 = 0.0d0
    end if

    gfx(1) = c0 * f0 + c1 * f1

    ! apply switching factors
    gfx(:) = gfx(:) * sfac(:)

end subroutine abf_accu_get_data_lsmooth

!===============================================================================
! Function:  abf_get_switching_factors
!===============================================================================

subroutine abf_get_switching_factors(vals)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP)     :: vals(:)
    ! -----------------------------------------------
    integer         :: i
    real(PMFDP)     :: r
    ! --------------------------------------------------------------------------

    sfac(:) = 0.0d0
    do i=1,NumOfABFCVs
        if( vals(i) .le. ABFCVList(i)%min_value ) cycle
        if( vals(i) .ge. ABFCVList(i)%max_value ) cycle

        if( vals(i) .lt. (ABFCVList(i)%min_value + ABFCVList(i)%buffer) ) then
            r = ABFCVList(i)%buffer - (vals(i) - ABFCVList(i)%min_value)
        else if ( vals(i) .ge. (ABFCVList(i)%max_value - ABFCVList(i)%buffer) ) then
            r = ABFCVList(i)%buffer + (vals(i) - ABFCVList(i)%max_value)
        else
            sfac(i) = 1.0d0
            cycle
        end if

    ! normalize
        r = r / ABFCVList(i)%buffer

    ! get switching function
        sfac(i) = 1.0d0 - 10.0d0*r**3 + 15.0d0*r**4 - 6.0d0*r**5

    end do

end subroutine abf_get_switching_factors

!===============================================================================

end module abf_accu

