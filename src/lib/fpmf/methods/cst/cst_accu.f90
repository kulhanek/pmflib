!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module cst_accu

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_accu_alloc
!===============================================================================

subroutine cst_accu_alloc

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer     :: i,tot_nbins,alloc_failed
    ! --------------------------------------------------------------------------

! accumulator setup for free energy calculation
    allocate( mlambda(NumOfAllCONs),   &
              m2lambda(NumOfAllCONs),  &
              stat= alloc_failed )

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,&
                 '[CST] Unable to allocate memory for arrays used in lambda calculation in cst_accu_alloc!')
    end if


! fdhtds  ----------------------------------------------------------------------
! accumulator setup for entropy and enthalpy

    if( fentropy ) then
        allocate( mlamtds(NumOfAllCONs),        &
                  m2lamtds(NumOfAllCONs),       &
                  mlamtdsfw(NumOfAllCONs),      &
                  m2lamtdsfw(NumOfAllCONs),     &
                  c11lt(NumOfAllCONs),          &
                  c11ltfw(NumOfAllCONs),        &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
        end if

        if( fentropy_decomp ) then
            allocate(   c11li(NumOfAllCONs),        &
                        c11lp(NumOfAllCONs),        &
                        c11lr(NumOfAllCONs),        &
                        c11lk(NumOfAllCONs),        &
                        stat= alloc_failed )

            if( alloc_failed .ne. 0 ) then
                call pmf_utils_exit(PMF_OUT,1,&
                         '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
            end if
        end if
    end if

    if( fintene .and. fintene_der ) then
        allocate( micf(NumOfAllCONs),           &
                  m2icf(NumOfAllCONs),          &
                  micffw(NumOfAllCONs),         &
                  m2icffw(NumOfAllCONs),        &
                  micfpfw(NumOfAllCONs),        &
                  m2icfpfw(NumOfAllCONs),       &
                  micfkfw(NumOfAllCONs),        &
                  m2icfkfw(NumOfAllCONs),       &
                  c11ii(NumOfAllCONs),          &
                  c11iifw(NumOfAllCONs),        &
                  stat= alloc_failed )

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,&
                     '[CST] Unable to allocate memory for arrays used for enthalpy/entropy calculations!')
        end if
    end if


! init PMF accu
    cstaccu%tot_cvs = NumOfCONs

    allocate(cstaccu%sizes(cstaccu%tot_cvs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CST accumulator!')
    endif

    tot_nbins       = 1
    fallconstant    = .true.

    do i=1,cstaccu%tot_cvs
        if( CONList(i)%mode .ne. 'C' ) then
            fallconstant = .false.
        end if
        if( CONList(i)%ibin .eq. 0 ) then
            freadranges = .false.   ! without all cv bins this cannot be enabled
        end if
    end do

    do i=1,cstaccu%tot_cvs
        if( freadranges ) then
            cstaccu%sizes(i)%min_value  = CONList(i)%min_value
            cstaccu%sizes(i)%max_value  = CONList(i)%max_value
            cstaccu%sizes(i)%nbins      = CONList(i)%nbins
            cstaccu%sizes(i)%width      = abs(cstaccu%sizes(i)%max_value - cstaccu%sizes(i)%min_value)
            cstaccu%sizes(i)%bin_width  = cstaccu%sizes(i)%width / cstaccu%sizes(i)%nbins
        else
            cstaccu%sizes(i)%min_value  = CONList(i)%value
            cstaccu%sizes(i)%max_value  = CONList(i)%value
            cstaccu%sizes(i)%nbins      = 1
            cstaccu%sizes(i)%width      = 0.0d0
            cstaccu%sizes(i)%bin_width  = 0.0d0
        end if
        cstaccu%sizes(i)%cv => CONList(i)%cv
        tot_nbins = tot_nbins * cstaccu%sizes(i)%nbins
    end do

    cstaccu%tot_nbins = tot_nbins

    allocate(   rbuf_B(cstaccu%tot_nbins),                    &
                rbuf_M(cstaccu%tot_cvs,cstaccu%tot_nbins),    &
                stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[CST] Unable to allocate memory for CST accumulator!')
    endif

! result the accumulate data
    call cst_accu_clear

end subroutine cst_accu_alloc

!===============================================================================
! Subroutine:  cst_accu_clear
!===============================================================================

subroutine cst_accu_clear

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    faccustep   = 0

    rbuf_B(:)   = 0.0d0
    rbuf_M(:,:) = 0.0d0

! free energy calculation
    nsamples    = 0.0d0

    mfw         = 0.0d0
    m2fw        = 0.0d0
    mlambda(:)  = 0.0d0
    m2lambda(:) = 0.0d0

! fdhtds  = enthalpy/entropy calculations
    if( fintene .or. fentropy ) then
        ntds            = 0.0d0
        fwsum           = 0.0d0
        fwsum2          = 0.0d0
        mfwtds          = 0.0d0
        m2fwtds         = 0.0d0
    end if

    if( fentropy ) then
        mlamtds(:)      = 0.0d0
        m2lamtds(:)     = 0.0d0
        mlamtdsfw(:)    = 0.0d0
        m2lamtdsfw(:)   = 0.0d0

        c11lt(:)    = 0.0d0
        c11ltfw(:)  = 0.0d0

        c11zh       = 0.0d0

        if( fentropy_decomp ) then
            c11li(:)  = 0.0d0
            c11lp(:)  = 0.0d0
            c11lr(:)  = 0.0d0
            c11lk(:)  = 0.0d0
        end if
    end if

    if( fentropy .or. fintene ) then
        metot       = 0.0d0
        m2etot      = 0.0d0
        metotfw     = 0.0d0
        m2etotfw    = 0.0d0
    end if

    if( (fentropy .and. fentropy_decomp) .or. fintene ) then
        meint       = 0.0d0
        m2eint      = 0.0d0
        meintfw     = 0.0d0
        m2eintfw    = 0.0d0
    end if

    if( fentropy .and. fentropy_decomp ) then
        mepot       = 0.0d0
        m2epot      = 0.0d0
        merst       = 0.0d0
        m2erst      = 0.0d0
        mekin       = 0.0d0
        m2ekin      = 0.0d0

        mepotfw     = 0.0d0
        m2epotfw    = 0.0d0
        merstfw     = 0.0d0
        m2erstfw    = 0.0d0
        mekinfw     = 0.0d0
        m2ekinfw    = 0.0d0
    end if

    if( fintene .and. fintene_der ) then
        micf(:)     = 0.0d0
        m2icf(:)    = 0.0d0
        micffw(:)   = 0.0d0
        m2icffw(:)  = 0.0d0

        micfpfw(:)  = 0.0d0
        m2icfpfw(:) = 0.0d0
        micfkfw(:)  = 0.0d0
        m2icfkfw(:) = 0.0d0

        c11ii(:)    = 0.0d0
        c11iifw(:)  = 0.0d0
    end if

end subroutine cst_accu_clear

!===============================================================================
! Subroutine:  cst_accu_read
!===============================================================================

subroutine cst_accu_read(iounit)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                         :: iounit
    !---------------------------------------------
    character(len=PMF_KEYLINE)      :: keyline
    integer                         :: i,glbidx,idx_local
    ! --------------------------------------------------------------------------

    if( freadranges ) then
        glbidx = 0
        do i=1,cstaccu%tot_cvs
            idx_local = CONList(i)%ibin - 1
            glbidx = glbidx*cstaccu%sizes(i)%nbins + idx_local
        end do
        glbidx = glbidx + 1
    else
        glbidx = 1
    end if

    do while(.true.)

        ! read keyline
        read(iounit,5,end=500,err=300) keyline

        ! process keyline
        if( pmf_accu_is_header_key(keyline) ) then
            call pmf_accu_read_header(cstaccu,iounit,'CST',keyline)
        else
            select case( pmf_accu_get_key(keyline) )
            ! ------------------------------------
                case('NSAMPLES')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    nsamples = int(rbuf_B(glbidx))
!            ! ------------------------------------
!                case('MLAMBDA')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    do i=1,cstaccu%tot_cvs
!                        mlambda(i) = rbuf_M(i,glbidx)
!                    end do
!            ! ------------------------------------
!                case('M2LAMBDA')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    do i=1,cstaccu%tot_cvs
!                        m2lambda(i) = rbuf_M(i,glbidx)
!                    end do
!            ! ------------------------------------
!                case('MISRZ')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    misrz = rbuf_B(glbidx)
!            ! ------------------------------------
!                case('M2ISRZ')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    m2isrz = rbuf_B(glbidx)
!
!! ------------------------------------
!
!                case('NTDS')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. fentropy ) then
!                        nsamples = int(rbuf_B(glbidx))
!                    end if
!
!! ------------------------------------
!
!                case('MEINT')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        meint = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('M2EINT')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        m2eint = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('MEPOT')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        mepot = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('M2EPOT')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        m2epot = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('MEKIN')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        mekin = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('M2EKIN')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        m2ekin = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('MERST')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        merst = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('M2ERST')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
!                        m2erst = rbuf_B(glbidx)
!                    end if
!
!! ------------------------------------
!
!                case('MICFP')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fenthalpy .and. fenthalpy_der ) then
!                        do i=1,cstaccu%tot_cvs
!                            micfp(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!            ! ------------------------------------
!                case('M2ICFP')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fenthalpy .and. fenthalpy_der ) then
!                        do i=1,cstaccu%tot_cvs
!                            m2icfp(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('C11PP')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fenthalpy .and. fenthalpy_der ) then
!                        do i=1,cstaccu%tot_cvs
!                            c11pp(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!
!! ------------------------------------
!
!                case('METOT')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fentropy ) then
!                        metot = rbuf_B(glbidx)
!                    end if
!            ! ------------------------------------
!                case('M2ETOT')
!                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
!                    if( fentropy ) then
!                        m2etot = rbuf_B(glbidx)
!                    end if
!           ! ------------------------------------
!                case('MPP')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            mpp(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('M2PP')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            m2pp(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('MPN')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            mpn(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('M2PN')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            m2pn(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('MHICF')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            mhicf(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('M2HICF')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            m2hicf(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!
!! ------------------------------------
!
!                case('C11HP')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            c11hp(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('C11HR')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            c11hr(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
!           ! ------------------------------------
!                case('C11HK')
!                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
!                    if( fentropy .and. fentdecomp ) then
!                        do i=1,cstaccu%tot_cvs
!                            c11hk(i) = rbuf_M(i,glbidx)
!                        end do
!                    end if
            ! ------------------------------------
                case default
                    call pmf_accu_skip_section(iounit,keyline,MTD_OUT)
            end select
        end if
    end do

    close(CST_RST)

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from the accumulator - keyline!')

end subroutine cst_accu_read

!===============================================================================
! Subroutine:  cst_accu_write_mean_M
!===============================================================================

subroutine cst_accu_write_mean_M(iounit,glbidx,mkey,mval,m2key,m2val,skey)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer         :: iounit
    integer         :: glbidx
    character(*)    :: mkey
    real(PMFDP)     :: mval(:)
    character(*)    :: m2key
    real(PMFDP)     :: m2val(:)
    character(*)    :: skey
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = mval(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,mkey, 'WA',rbuf_M,skey)

    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = m2val(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,m2key,'M2',rbuf_M,skey,mkey)

end subroutine cst_accu_write_mean_M

!===============================================================================
! Subroutine:  cst_accu_write_cmom_M
!===============================================================================

subroutine cst_accu_write_cmom_M(iounit,glbidx,mkey,mval,nkey,akey,bkey)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer         :: iounit
    integer         :: glbidx
    character(*)    :: mkey
    real(PMFDP)     :: mval(:)
    character(*)    :: nkey
    character(*)    :: akey
    character(*)    :: bkey
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = mval(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,mkey, 'CO',rbuf_M,nkey,akey,bkey)

end subroutine cst_accu_write_cmom_M

!===============================================================================
! Subroutine:  cst_accu_write_mean_B
!===============================================================================

subroutine cst_accu_write_mean_B(iounit,glbidx,mkey,mval,m2key,m2val,skey)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer         :: iounit
    integer         :: glbidx
    character(*)    :: mkey
    real(PMFDP)     :: mval
    character(*)    :: m2key
    real(PMFDP)     :: m2val
    character(*)    :: skey
    ! --------------------------------------------------------------------------

    rbuf_B(:)       = 0.0d0
    rbuf_B(glbidx)  = mval
    call pmf_accu_write_rbuf_B(cstaccu,iounit,mkey, 'WA',rbuf_B,skey)

    rbuf_B(glbidx)  = m2val
    call pmf_accu_write_rbuf_B(cstaccu,iounit,m2key,'M2',rbuf_B,skey,mkey)

end subroutine cst_accu_write_mean_B

!===============================================================================
! Subroutine:  cst_accu_write_cmom_B
!===============================================================================

subroutine cst_accu_write_cmom_B(iounit,glbidx,mkey,mval,nkey,akey,bkey)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer         :: iounit
    integer         :: glbidx
    character(*)    :: mkey
    real(PMFDP)     :: mval
    character(*)    :: nkey
    character(*)    :: akey
    character(*)    :: bkey
    ! --------------------------------------------------------------------------

    rbuf_B(glbidx) = mval
    call pmf_accu_write_rbuf_B(cstaccu,iounit,mkey, 'CO',rbuf_B,nkey,akey,bkey)

end subroutine cst_accu_write_cmom_B

!===============================================================================
! Subroutine:  cst_accu_write_counter_B
!===============================================================================

subroutine cst_accu_write_counter_B(iounit,glbidx,ckey,cval)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer         :: iounit
    integer         :: glbidx
    character(*)    :: ckey
    real(PMFDP)     :: cval
    ! --------------------------------------------------------------------------

    rbuf_B(:)       = 0.0d0
    rbuf_B(glbidx)  = cval
    call pmf_accu_write_rbuf_B(cstaccu,iounit,ckey, 'AD',rbuf_B)

end subroutine cst_accu_write_counter_B

!===============================================================================
! Subroutine:  cst_accu_write
!===============================================================================

subroutine cst_accu_write(iounit)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer     :: iounit
    ! --------------------------------------------
    integer     :: i,glbidx,idx_local
    !---------------------------------------------------------------------------

    cstaccu%method = 'CST'
    call pmf_accu_write_header(cstaccu,iounit)

    if( freadranges ) then
        glbidx = 0
        do i=1,cstaccu%tot_cvs
            idx_local = CONList(i)%ibin - 1
            glbidx = glbidx*cstaccu%sizes(i)%nbins + idx_local
        end do
        glbidx = glbidx + 1
    else
        glbidx = 1
    end if

! ------------------------------------------------

    call cst_accu_write_counter_B(iounit,glbidx,'NSAMPLES', nsamples)
    call cst_accu_write_mean_M(iounit,glbidx,   'MLAMBDA',  mlambda,'M2LAMBDA', m2lambda,'NSAMPLES')
    call cst_accu_write_mean_B(iounit,glbidx,   'MFW',      mfw,    'M2FW',     m2fw,    'NSAMPLES')

! fdhtds  ----------------------------------------------------------------------

   if( fintene .or. fentropy ) then
        call cst_accu_write_counter_B(iounit,glbidx,'NTDS',     ntds)
        call cst_accu_write_counter_B(iounit,glbidx,'FWSUM',    fwsum)
        call cst_accu_write_counter_B(iounit,glbidx,'FWSUM2',   fwsum2)

        call cst_accu_write_mean_B(iounit,glbidx,'MFWTDS',mfwtds,'M2FWTDS',m2fwtds,   'NTDS')
    end if

    if( fentropy ) then
        call cst_accu_write_mean_M(iounit,glbidx,'MLAMTDS',   mlamtds,   'M2LAMTDS',   m2lamtds,   'NTDS')
        call cst_accu_write_mean_M(iounit,glbidx,'MLAMTDSFW', mlamtdsfw, 'M2LAMTDSFW', m2lamtdsfw, 'FWSUM')

        call cst_accu_write_cmom_M(iounit,glbidx,'C11LT',   c11lt,   'NTDS',  'MLAMTDS',   'METOT')
        call cst_accu_write_cmom_M(iounit,glbidx,'C11LTFW', c11ltfw, 'FWSUM', 'MLAMTDSFW', 'METOTFW')

        call cst_accu_write_cmom_B(iounit,glbidx,'C11ZH',   c11zh,   'NTDS',  'MFWTDS',    'METOT')

        if( fentropy_decomp ) then
            call cst_accu_write_cmom_M(iounit,glbidx,'C11LI', c11li, 'NTDS', 'MLAMTDS', 'MEINT')
            call cst_accu_write_cmom_M(iounit,glbidx,'C11LP', c11lp, 'NTDS', 'MLAMTDS', 'MEPOT')
            call cst_accu_write_cmom_M(iounit,glbidx,'C11LR', c11lr, 'NTDS', 'MLAMTDS', 'MERST')
            call cst_accu_write_cmom_M(iounit,glbidx,'C11LK', c11lk, 'NTDS', 'MLAMTDS', 'MEKIN')
        end if
    end if

    if( fentropy .or. fintene ) then
        call cst_accu_write_mean_B(iounit,glbidx,'METOT',metot,'M2ETOT',m2etot,'NTDS')
        call cst_accu_write_mean_B(iounit,glbidx,'METOTFW',metotfw,'M2ETOTFW',m2etotfw,'FWSUM')
    end if

    if( (fentropy .and. fentropy_decomp) .or. fintene ) then
        call cst_accu_write_mean_B(iounit,glbidx,'MEINT',meint,'M2EINT',m2eint,'NTDS')
        call cst_accu_write_mean_B(iounit,glbidx,'MEINTFW',meintfw,'M2EINTFW',m2eintfw,'FWSUM')
    end if

    if( fentropy .and. fentropy_decomp ) then
        call cst_accu_write_mean_B(iounit,glbidx,'MEPOT',mepot,'M2EPOT',m2epot,'NTDS')
        call cst_accu_write_mean_B(iounit,glbidx,'MERST',merst,'M2ERST',m2erst,'NTDS')
        call cst_accu_write_mean_B(iounit,glbidx,'MEKIN',mekin,'M2EKIN',m2ekin,'NTDS')

        call cst_accu_write_mean_B(iounit,glbidx,'MEPOTFW',mepotfw,'M2EPOTFW',m2epotfw,'FWSUM')
        call cst_accu_write_mean_B(iounit,glbidx,'MERSTFW',merstfw,'M2ERSTFW',m2erstfw,'FWSUM')
        call cst_accu_write_mean_B(iounit,glbidx,'MEKINFW',mekinfw,'M2EKINFW',m2ekinfw,'FWSUM')
    end if

    if( fintene .and. fintene_der ) then
        call cst_accu_write_mean_M(iounit,glbidx,'MICF',   micf,   'M2ICF',   m2icf,   'NTDS')
        call cst_accu_write_mean_M(iounit,glbidx,'MICFFW', micffw, 'M2ICFFW', m2icffw, 'FWSUM')
        call cst_accu_write_mean_M(iounit,glbidx,'MICFPFW',micfpfw,'M2ICFPFW',m2icfpfw,'FWSUM')
        call cst_accu_write_mean_M(iounit,glbidx,'MICFKFW',micfkfw,'M2ICFKFW',m2icfkfw,'FWSUM')

        call cst_accu_write_cmom_M(iounit,glbidx,'C11II',   c11ii,   'NTDS',  'MICF',    'MEINT')
        call cst_accu_write_cmom_M(iounit,glbidx,'C11IIFW', c11iifw, 'FWSUM', 'MICFFW',  'MEINTFW')
    end if

end subroutine cst_accu_write

!===============================================================================
! Subroutine:  cst_accu_add_data_OM
! online mean and M2
!===============================================================================

subroutine cst_accu_add_data_OM(ival,invn,mval,m2val)

    implicit none
    real(PMFDP) :: ival
    real(PMFDP) :: invn
    real(PMFDP) :: mval
    real(PMFDP) :: m2val
    ! --------------------------------------------
    real(PMFDP) :: dval1,dval2
    ! --------------------------------------------------------------------------

    dval1 = ival  - mval
    mval  = mval  + dval1 * invn
    dval2 = ival  - mval
    m2val = m2val + dval1 * dval2

end subroutine cst_accu_add_data_OM

!===============================================================================
! Subroutine:  cst_accu_add_data_WOM
! online weighted mean and M2
!===============================================================================

subroutine cst_accu_add_data_WOM(ival,invw,w,mval,m2val)

    implicit none
    real(PMFDP) :: ival
    real(PMFDP) :: invw   ! w / wsum
    real(PMFDP) :: w
    real(PMFDP) :: mval
    real(PMFDP) :: m2val
    ! --------------------------------------------
    real(PMFDP) :: dval1,dval2
    ! --------------------------------------------------------------------------

    dval1 = ival  - mval
    mval  = mval  + dval1 * invw
    dval2 = ival  - mval
    m2val = m2val + w * dval1 * dval2

end subroutine cst_accu_add_data_WOM

!===============================================================================
! Subroutine:  cst_accu_add_data_OM
! online mean and M2
!===============================================================================

subroutine cst_accu_add_data_OMI(ival,invn,mval,m2val,dval1,dval2)

    implicit none
    real(PMFDP) :: ival
    real(PMFDP) :: invn
    real(PMFDP) :: mval
    real(PMFDP) :: m2val
    real(PMFDP) :: dval1
    real(PMFDP) :: dval2
    ! --------------------------------------------------------------------------

    dval1 = ival  - mval
    mval  = mval  + dval1 * invn
    dval2 = ival  - mval
    m2val = m2val + dval1 * dval2

end subroutine cst_accu_add_data_OMI

!===============================================================================
! Subroutine:  cst_accu_add_data_WOMI
! online weighted mean and M2
!===============================================================================

subroutine cst_accu_add_data_WOMI(ival,invw,w,mval,m2val,dval1,dval2)

    implicit none
    real(PMFDP) :: ival
    real(PMFDP) :: invw   ! w / wsum
    real(PMFDP) :: w
    real(PMFDP) :: mval
    real(PMFDP) :: m2val
    real(PMFDP) :: dval1
    real(PMFDP) :: dval2
    ! --------------------------------------------------------------------------

    dval1 = ival  - mval
    mval  = mval  + dval1 * invw
    dval2 = ival  - mval
    m2val = m2val + w * dval1 * dval2

end subroutine cst_accu_add_data_WOMI

!===============================================================================
! Subroutine:  cst_accu_add_lam
! free energy
!===============================================================================

subroutine cst_accu_add_lam

    use pmf_dat
    use cst_dat

    implicit none
    integer         :: i
    real(PMFDP)     :: invn, llam, lfw
    ! --------------------------------------------------------------------------

    if( mod(fstep,flamsample) .ne. 0 ) return

    nsamples = nsamples + 1
    if( nsamples .le. 0 ) return
    invn = 1.0d0/nsamples

    do i=1,NumOfAllCONs
        llam = lambdahist(i,hist_len+hist_fidx)
        call cst_accu_add_data_OM(llam,invn,mlambda(i),m2lambda(i))
    end do

    lfw = fwhist(hist_len+hist_fidx)
    call cst_accu_add_data_OM(lfw,invn,mfw,m2fw)

end subroutine cst_accu_add_lam

!===============================================================================
! Subroutine:  cst_accu_add_dhTds
! enthalpy and entropy
!===============================================================================

subroutine cst_accu_add_dhTds

    use pmf_dat
    use cst_dat

    implicit none
    integer         :: i
    real(PMFDP)     :: invn,invw
    real(PMFDP)     :: lfw,llam,licf,licfp,licfk
    real(PMFDP)     :: letot,leint,lepot,lerst,lekin
    real(PMFDP)     :: detot1,detot2
    real(PMFDP)     :: deint1,deint2
    real(PMFDP)     :: detot1fw,detot2fw
    real(PMFDP)     :: deint1fw,deint2fw
    real(PMFDP)     :: depot1,depot2
    real(PMFDP)     :: derst1,derst2
    real(PMFDP)     :: dekin1,dekin2
    real(PMFDP)     :: dicf1,dicf2
    real(PMFDP)     :: dlam1,dlam2
    real(PMFDP)     :: dicf1fw,dicf2fw
    real(PMFDP)     :: dlam1fw,dlam2fw
    real(PMFDP)     :: dfw1,dfw2
    ! --------------------------------------------------------------------------

    if( .not. (fintene .or. fentropy) ) return
    if( enevalidhist(hist_len+hist_fidx) ) faccustep = faccustep + 1
    if( .not. ( (mod(faccustep,fenesample) .eq. 0) .and. enevalidhist(hist_len+hist_fidx) ) ) return
    if( fstep .le. 3*hist_len ) return

    ntds = ntds + 1.0d0
    invn = 1.0d0/ntds

    lfw = fwhist(hist_len+hist_fidx)

    fwsum   = fwsum + lfw
    invw    = lfw / fwsum
    fwsum2  = fwsum2 + lfw*lfw

    call cst_accu_add_data_OMI(lfw,invn,mfwtds,m2fwtds,dfw1,dfw2)

! other data
    lepot        = epothist(hist_len+hist_fidx)
    lerst        = ersthist(hist_len+hist_fidx)
    lekin        = ekinhist(hist_len+hist_fidx)
    letot        = lepot + lerst + lekin
    leint        = lepot + lerst

   ! write(12478,*) fstep, lepot, lekin

    call cst_accu_add_data_OMI(letot,invn,metot,m2etot,detot1,detot2)
    call cst_accu_add_data_OMI(leint,invn,meint,m2eint,deint1,deint2)
    call cst_accu_add_data_WOMI(letot,invw,lfw,metotfw,m2etotfw,detot1fw,detot2fw)
    call cst_accu_add_data_WOMI(leint,invw,lfw,meintfw,m2eintfw,deint1fw,deint2fw)

    if( fentropy ) then
        c11zh = c11zh + dfw1 * detot2

        if( fentropy_decomp ) then
            call cst_accu_add_data_OMI(lepot,invn,mepot,m2epot,depot1,depot2)
            call cst_accu_add_data_OMI(lerst,invn,merst,m2erst,derst1,derst2)
            call cst_accu_add_data_OMI(lekin,invn,mekin,m2ekin,dekin1,dekin2)

            call cst_accu_add_data_WOM(lepot,invw,lfw,mepotfw,m2epotfw)
            call cst_accu_add_data_WOM(lerst,invw,lfw,merstfw,m2erstfw)
            call cst_accu_add_data_WOM(lekin,invw,lfw,mekinfw,m2ekinfw)
        end if
    end if

    do i=1,NumOfAllCONs

        select case(flambda_src)
            case(0)
                llam  = lambdahist(i,hist_len+hist_fidx+flambda_lag)
            case(1)
                llam  = lambda1hist(i,hist_len+hist_fidx+flambda_lag)
        end select

        if( fentropy ) then
            call cst_accu_add_data_OMI(llam, invn, mlamtds(i), m2lamtds(i), dlam1, dlam2)
            call cst_accu_add_data_WOMI(llam, invw, lfw, mlamtdsfw(i), m2lamtdsfw(i), dlam1fw, dlam2fw)

            c11lt(i)    = c11lt(i)      + dlam1 * detot2
            c11ltfw(i)  = c11ltfw(i)    + lfw * dlam1fw * detot2fw

            if( fentropy_decomp ) then
                c11li(i)  = c11li(i)    +  dlam1 * deint2
                c11lp(i)  = c11lp(i)    +  dlam1 * depot2
                c11lr(i)  = c11lr(i)    +  dlam1 * derst2
                c11lk(i)  = c11lk(i)    +  dlam1 * dekin2
            end if
        end if

        if( fintene .and. fintene_der ) then
            licfp = icfphist(i,hist_len+hist_fidx)
            licfk = - PMF_Rgas*ftemp * icfkhist(i,hist_len+hist_fidx)
            licf  = licfp + licfk

            call cst_accu_add_data_OMI(licf, invn, micf(i),    m2icf(i),    dicf1, dicf2)
            call cst_accu_add_data_WOMI(licf, invw, lfw, micffw(i),    m2icffw(i),    dicf1fw, dicf2fw)
            call cst_accu_add_data_WOM(licfp, invw, lfw, micfpfw(i),   m2icfpfw(i))
            call cst_accu_add_data_WOM(licfk, invw, lfw, micfkfw(i),   m2icfkfw(i))

            c11ii(i)    = c11ii(i)      +  dicf1   * deint2
            c11iifw(i)  = c11iifw(i)    +  lfw * dicf1fw * deint2fw
        end if
    end do

end subroutine cst_accu_add_dhTds

!===============================================================================

end module cst_accu
