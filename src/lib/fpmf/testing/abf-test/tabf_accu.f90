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

module tabf_accu

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  tabf_accu_init
!===============================================================================

subroutine tabf_accu_init()

    use tabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    tabfaccu%tot_cvs = NumOfTABFCVs
    allocate(tabfaccu%sizes(tabfaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[TABF] Unable to allocate memory for abf accumulator!')
    endif

    tot_nbins = 1
    do i=1,tabfaccu%tot_cvs
        tabfaccu%sizes(i)%min_value  = TABFCVList(i)%min_value
        tabfaccu%sizes(i)%max_value  = TABFCVList(i)%max_value
        tabfaccu%sizes(i)%nbins      = TABFCVList(i)%nbins
        tabfaccu%sizes(i)%width      = abs(tabfaccu%sizes(i)%max_value - tabfaccu%sizes(i)%min_value)
        tabfaccu%sizes(i)%bin_width  = tabfaccu%sizes(i)%width / tabfaccu%sizes(i)%nbins
        tabfaccu%sizes(i)%cv         => TABFCVList(i)%cv
        tot_nbins = tot_nbins * tabfaccu%sizes(i)%nbins
    end do

    tabfaccu%tot_nbins = tot_nbins

    ! ABF force arrays
    allocate(  &
            tabfaccu%nsamples(tabfaccu%tot_nbins), &
            tabfaccu%micf(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
            tabfaccu%m2icf(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
            stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[TABF] Unable to allocate memory for abf accumulator (icf)!')
    endif

! enthalpy ---------------------------------------------------------------------
    if( fenthalpy .or. fentropy ) then
        allocate(  &
                tabfaccu%metot(tabfaccu%tot_nbins), &
                tabfaccu%m2etot(tabfaccu%tot_nbins), &
                tabfaccu%mepot(tabfaccu%tot_nbins), &
                tabfaccu%m2epot(tabfaccu%tot_nbins), &
                tabfaccu%mekin(tabfaccu%tot_nbins), &
                tabfaccu%m2ekin(tabfaccu%tot_nbins), &
                tabfaccu%merst(tabfaccu%tot_nbins), &
                tabfaccu%m2erst(tabfaccu%tot_nbins), &
                stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[TABF] Unable to allocate memory for abf accumulator (emthalpy)!')
        endif
    end if

! entropy ----------------------------------------------------------------------
    if( fentropy ) then
        allocate(  &
                tabfaccu%micf_pot(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%m2icf_pot(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%micf_kin(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%m2icf_kin(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &

                tabfaccu%c11hh(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%c11pp(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%c11pk(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%c11pr(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%c11kp(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%c11kk(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                tabfaccu%c11kr(tabfaccu%tot_cvs,tabfaccu%tot_nbins), &
                stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[TABF] Unable to allocate memory for abf accumulator (abfforce)!')
        endif
    end if

    call tabf_accu_clear()

    return

end subroutine tabf_accu_init

!===============================================================================
! Subroutine:  tabf_accu_clear
!===============================================================================

subroutine tabf_accu_clear()

    use tabf_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    tabfaccu%nsamples(:)     = 0

    tabfaccu%micf(:,:)       = 0.0d0
    tabfaccu%m2icf(:,:)      = 0.0d0

    if( fenthalpy .or. fentropy ) then
        tabfaccu%metot(:)        = 0.0d0
        tabfaccu%m2etot(:)       = 0.0d0
        tabfaccu%mepot(:)        = 0.0d0
        tabfaccu%m2epot(:)       = 0.0d0
        tabfaccu%mekin(:)        = 0.0d0
        tabfaccu%m2ekin(:)       = 0.0d0
        tabfaccu%merst(:)        = 0.0d0
        tabfaccu%m2erst(:)       = 0.0d0
    end if

    if( fentropy ) then
        tabfaccu%micf_pot(:,:)   = 0.0d0
        tabfaccu%m2icf_pot(:,:)  = 0.0d0
        tabfaccu%micf_kin(:,:)   = 0.0d0
        tabfaccu%m2icf_kin(:,:)  = 0.0d0

        tabfaccu%c11hh(:,:)      = 0.0d0
        tabfaccu%c11pp(:,:)      = 0.0d0
        tabfaccu%c11pk(:,:)      = 0.0d0
        tabfaccu%c11pr(:,:)      = 0.0d0
        tabfaccu%c11kp(:,:)      = 0.0d0
        tabfaccu%c11kk(:,:)      = 0.0d0
        tabfaccu%c11kr(:,:)      = 0.0d0
    end if

end subroutine tabf_accu_clear

!===============================================================================
! Subroutine:  tabf_accu_write
!===============================================================================

subroutine tabf_accu_write(iounit)

    use tabf_dat

    implicit none
    integer                     :: iounit
    !---------------------------------------------------------------------------

    tabfaccu%method = 'TABF'
    call pmf_accu_write_header(tabfaccu%PMFAccuType,iounit)
    call pmf_accu_write_ibuf_B(tabfaccu%PMFAccuType,iounit,'NSAMPLES',  'AD',tabfaccu%nsamples)
    call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'MICF',      'WA',tabfaccu%micf)
    call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'M2ICF',     'AD',tabfaccu%m2icf)

    if( fentropy ) then
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'MICF_POT',  'WA',tabfaccu%micf_pot)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'M2ICF_POT', 'AD',tabfaccu%m2icf_pot)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'MICF_KIN',  'WA',tabfaccu%micf_kin)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'M2ICF_KIN', 'AD',tabfaccu%m2icf_kin)
    end if

    if( fenthalpy .or. fentropy ) then
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'METOT','WA',tabfaccu%metot)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'M2ETOT','AD',tabfaccu%m2etot)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'MEPOT','WA',tabfaccu%mepot)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'M2EPOT','AD',tabfaccu%m2epot)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'MEKIN','WA',tabfaccu%mekin)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'M2EKIN','AD',tabfaccu%m2ekin)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'MERST','WA',tabfaccu%merst)
        call pmf_accu_write_rbuf_B(tabfaccu%PMFAccuType,iounit,'M2ERST','AD',tabfaccu%m2erst)
    end if

    if( fentropy ) then
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11HH','AD',tabfaccu%c11hh)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11PP','AD',tabfaccu%c11pp)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11PK','AD',tabfaccu%c11pk)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11PR','AD',tabfaccu%c11pr)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11KP','AD',tabfaccu%c11kp)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11KK','AD',tabfaccu%c11kk)
        call pmf_accu_write_rbuf_M(tabfaccu%PMFAccuType,iounit,'C11KR','AD',tabfaccu%c11kr)
    end if

end subroutine tabf_accu_write

!===============================================================================
! Subroutine:  tabf_accu_add_data_online
!===============================================================================

subroutine tabf_accu_add_data_online(cvs,gfx_pot,gfx_kin,epot,ekin,erst)

    use tabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: gfx_pot(:)
    real(PMFDP)    :: gfx_kin(:)
    real(PMFDP)    :: epot
    real(PMFDP)    :: ekin
    real(PMFDP)    :: erst
    ! -----------------------------------------------
    integer        :: gi0, i
    real(PMFDP)    :: invn, etot, icf
    real(PMFDP)    :: detot1, detot2
    real(PMFDP)    :: depot1, depot2
    real(PMFDP)    :: dekin1, dekin2
    real(PMFDP)    :: derst1, derst2
    real(PMFDP)    :: dicf1, dicf2
    real(PMFDP)    :: dicf_pot1, dicf_pot2
    real(PMFDP)    :: dicf_kin1, dicf_kin2
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(tabfaccu%PMFAccuType,cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    tabfaccu%nsamples(gi0) = tabfaccu%nsamples(gi0) + 1
    invn = 1.0d0 / real(tabfaccu%nsamples(gi0),PMFDP)

    if( fenthalpy .or. fentropy ) then
        etot = epot + ekin + erst

        ! total energy
        detot1 = etot - tabfaccu%metot(gi0)
        tabfaccu%metot(gi0)  = tabfaccu%metot(gi0)  + detot1 * invn
        detot2 = etot - tabfaccu%metot(gi0)
        tabfaccu%m2etot(gi0) = tabfaccu%m2etot(gi0) + detot1 * detot2

        ! potential energy
        depot1 = epot - tabfaccu%mepot(gi0)
        tabfaccu%mepot(gi0)  = tabfaccu%mepot(gi0)  + depot1 * invn
        depot2 = epot - tabfaccu%mepot(gi0)
        tabfaccu%m2epot(gi0) = tabfaccu%m2epot(gi0) + depot1 * depot2

        ! potential energy
        dekin1 = ekin - tabfaccu%mekin(gi0)
        tabfaccu%mekin(gi0)  = tabfaccu%mekin(gi0)  + dekin1 * invn
        dekin2 = ekin - tabfaccu%mekin(gi0)
        tabfaccu%m2ekin(gi0) = tabfaccu%m2ekin(gi0) + dekin1 * dekin2

        ! restraint energy
        derst1 = erst - tabfaccu%merst(gi0)
        tabfaccu%merst(gi0)  = tabfaccu%merst(gi0)  + derst1 * invn
        derst2 = erst - tabfaccu%merst(gi0)
        tabfaccu%m2erst(gi0) = tabfaccu%m2erst(gi0) + derst1 * derst2
    end if

    do i=1,NumOfTABFCVs
        icf = gfx_pot(i) + gfx_kin(i)

        dicf1 = - icf - tabfaccu%micf(i,gi0)
        tabfaccu%micf(i,gi0)  = tabfaccu%micf(i,gi0)  + dicf1 * invn
        dicf2 = - icf -  tabfaccu%micf(i,gi0)
        tabfaccu%m2icf(i,gi0) = tabfaccu%m2icf(i,gi0) + dicf1 * dicf2

        if( fentropy ) then
            dicf_pot1 = - gfx_pot(i) - tabfaccu%micf_pot(i,gi0)
            tabfaccu%micf_pot(i,gi0)  = tabfaccu%micf_pot(i,gi0)  + dicf_pot1 * invn
            dicf_pot2 = - gfx_pot(i) -  tabfaccu%micf_pot(i,gi0)
            tabfaccu%m2icf_pot(i,gi0) = tabfaccu%m2icf_pot(i,gi0) + dicf_pot1 * dicf_pot2

            dicf_kin1 = - gfx_kin(i) - tabfaccu%micf_kin(i,gi0)
            tabfaccu%micf_kin(i,gi0)  = tabfaccu%micf_kin(i,gi0)  + dicf_kin1 * invn
            dicf_kin2 = - gfx_kin(i) -  tabfaccu%micf_kin(i,gi0)
            tabfaccu%m2icf_kin(i,gi0) = tabfaccu%m2icf_kin(i,gi0) + dicf_kin1 * dicf_kin2

            tabfaccu%c11hh(i,gi0)  = tabfaccu%c11hh(i,gi0) + dicf1     * detot2
            tabfaccu%c11pp(i,gi0)  = tabfaccu%c11pp(i,gi0) + dicf_pot1 * depot2
            tabfaccu%c11pk(i,gi0)  = tabfaccu%c11pk(i,gi0) + dicf_pot1 * dekin2
            tabfaccu%c11pr(i,gi0)  = tabfaccu%c11pr(i,gi0) + dicf_pot1 * derst2
            tabfaccu%c11kp(i,gi0)  = tabfaccu%c11kp(i,gi0) + dicf_kin1 * depot2
            tabfaccu%c11kk(i,gi0)  = tabfaccu%c11kk(i,gi0) + dicf_kin1 * dekin2
            tabfaccu%c11kr(i,gi0)  = tabfaccu%c11kr(i,gi0) + dicf_kin1 * derst2
        end if
    end do

end subroutine tabf_accu_add_data_online

!===============================================================================
! Subroutine:  tabf_accu_get_data_lramp
!===============================================================================

subroutine tabf_accu_get_data_lramp(values,gfx)

    use tabf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0, n
    real(PMFDP)    :: sc_ramp
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = pmf_accu_globalindex(tabfaccu%PMFAccuType,values)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get number of samples
    n = tabfaccu%nsamples(gi0)
    if( n .gt. 0 ) then
        sc_ramp = 1.0d0
        if( n .le. fhramp_max ) then
            sc_ramp = 0.0d0
            if( n .gt. fhramp_min ) then
                sc_ramp = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        gfx(:) = sc_ramp * tabfaccu%micf(:,gi0)
    end if

end subroutine tabf_accu_get_data_lramp

!===============================================================================

end module tabf_accu

