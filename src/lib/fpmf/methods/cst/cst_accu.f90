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
            ! ------------------------------------
                case('MLAMBDA')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    do i=1,cstaccu%tot_cvs
                        mlambda(i) = rbuf_M(i,glbidx)
                    end do
            ! ------------------------------------
                case('M2LAMBDA')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    do i=1,cstaccu%tot_cvs
                        m2lambda(i) = rbuf_M(i,glbidx)
                    end do
            ! ------------------------------------
                case('MISRZ')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    misrz = rbuf_B(glbidx)
            ! ------------------------------------
                case('M2ISRZ')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    m2isrz = rbuf_B(glbidx)

! ------------------------------------

                case('NTDS')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. fentropy ) then
                        nsamples = int(rbuf_B(glbidx))
                    end if

! ------------------------------------

                case('MEINT')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        meint = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('M2EINT')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        m2eint = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('MEPOT')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        mepot = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('M2EPOT')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        m2epot = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('MEKIN')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        mekin = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('M2EKIN')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        m2ekin = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('MERST')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        merst = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('M2ERST')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
                        m2erst = rbuf_B(glbidx)
                    end if

! ------------------------------------

                case('MICFP')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fenthalpy .and. fenthalpy_der ) then
                        do i=1,cstaccu%tot_cvs
                            micfp(i) = rbuf_M(i,glbidx)
                        end do
                    end if
            ! ------------------------------------
                case('M2ICFP')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fenthalpy .and. fenthalpy_der ) then
                        do i=1,cstaccu%tot_cvs
                            m2icfp(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('C11PP')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fenthalpy .and. fenthalpy_der ) then
                        do i=1,cstaccu%tot_cvs
                            c11pp(i) = rbuf_M(i,glbidx)
                        end do
                    end if

! ------------------------------------

                case('METOT')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fentropy ) then
                        metot = rbuf_B(glbidx)
                    end if
            ! ------------------------------------
                case('M2ETOT')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    if( fentropy ) then
                        m2etot = rbuf_B(glbidx)
                    end if
           ! ------------------------------------
                case('MPP')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            mpp(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('M2PP')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            m2pp(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('MPN')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            mpn(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('M2PN')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            m2pn(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('MHICF')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            mhicf(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('M2HICF')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            m2hicf(i) = rbuf_M(i,glbidx)
                        end do
                    end if

! ------------------------------------

                case('C11HP')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            c11hp(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('C11HR')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            c11hr(i) = rbuf_M(i,glbidx)
                        end do
                    end if
           ! ------------------------------------
                case('C11HK')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( fentropy .and. fentdecomp ) then
                        do i=1,cstaccu%tot_cvs
                            c11hk(i) = rbuf_M(i,glbidx)
                        end do
                    end if
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

    rbuf_B(:) = 0
    rbuf_B(glbidx) = nsamples
    call pmf_accu_write_rbuf_B(cstaccu,iounit,'NSAMPLES','AD',rbuf_B)

! ------------------------------------------------
    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = mlambda(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,'MLAMBDA','WA',rbuf_M,'NSAMPLES')

    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = m2lambda(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2LAMBDA','M2',rbuf_M,'NSAMPLES','MLAMBDA')

! ------------------------------------------------
    rbuf_B(:) = 0.0d0
    rbuf_B(glbidx) = misrz
    call pmf_accu_write_rbuf_B(cstaccu,iounit,'MISRZ','WA',rbuf_B,'NSAMPLES')

    rbuf_B(:) = 0.0d0
    rbuf_B(glbidx) = m2isrz
    call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2ISRZ','M2',rbuf_B,'NSAMPLES','MISRZ')

! ------------------------------------------------

    if( fenthalpy .or. fentropy ) then
        rbuf_B(:) = 0
        rbuf_B(glbidx) = ntds
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'NTDS',   'AD',rbuf_B)
    end if

    if( fenthalpy .or. (fentropy .and. fentdecomp) ) then
        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = meint
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MEINT',  'WA',rbuf_B, 'NTDS')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2eint
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2EINT', 'M2',rbuf_B, 'NTDS','MEINT')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = mepot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MEPOT',  'WA',rbuf_B, 'NTDS')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2epot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2EPOT', 'M2',rbuf_B, 'NTDS','MEPOT')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = merst
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MERST',  'WA',rbuf_B, 'NTDS')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2erst
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2ERST', 'M2',rbuf_B, 'NTDS','MERST')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = mekin
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MEKIN',  'WA',rbuf_B, 'NTDS')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2ekin
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2EKIN', 'M2',rbuf_B, 'NTDS','MEKIN')
    end if

    if( fenthalpy .and. fenthalpy_der ) then
        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = micfp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'MICFP',  'WA',rbuf_M, 'NTDS')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m2icfp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2ICFP', 'M2',rbuf_M, 'NTDS','MICFP')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = micfpz(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'MICFPZ', 'WA',rbuf_M, 'NTDS')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m2icfpz(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2ICFPZ','M2',rbuf_M, 'NTDS','MICFPZ')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = mfixmanw
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MFIXW',  'WA',rbuf_B, 'NTDS')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2fixmanw
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2FIXW', 'M2',rbuf_B, 'NTDS','MFIXW')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11pp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11PP',  'CO',rbuf_M, 'NTDS','MICFP','MEINT')
    end if

    if( fentropy ) then
        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = metot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'METOT',  'WA',rbuf_B, 'NTDS')

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2etot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2ETOT', 'M2',rbuf_B, 'NTDS','METOT')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = mpp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'MPP',    'WA',rbuf_M, 'NTDS')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m2pp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2PP',   'M2',rbuf_M, 'NTDS','MPP')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = mpn(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'MPN',    'WA',rbuf_M, 'NTDS')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m2pn(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2PN',   'M2',rbuf_M, 'NTDS','MPN')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = mhicf(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'MHICF',  'WA',rbuf_M, 'NTDS')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m2hicf(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2HICF', 'M2',rbuf_M, 'NTDS','MHICF')
    end if

    if( fentropy .and. fentdecomp ) then
        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HP','CO',rbuf_M, 'NTDS','MHICF','MEPOT')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hr(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HR','CO',rbuf_M, 'NTDS','MHICF','MERST')

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hk(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HK','CO',rbuf_M, 'NTDS','MHICF','MEKIN')
    end if

end subroutine cst_accu_write

!===============================================================================

end module cst_accu
