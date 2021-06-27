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
                    call pmf_accu_read_ibuf_B(cstaccu,iounit,keyline,ibuf_B)
                    faccumulation = ibuf_B(glbidx)
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
            ! ------------------------------------                             ! ------------------------------------
                case('MLAMBDAV')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( has_lambdav ) then
                        do i=1,cstaccu%tot_cvs
                            mlambdav(i) = rbuf_M(i,glbidx)
                        end do
                   end if
            ! ------------------------------------
                case('M2LAMBDAV')
                    call pmf_accu_read_rbuf_M(cstaccu,iounit,keyline,rbuf_M)
                    if( has_lambdav ) then
                        do i=1,cstaccu%tot_cvs
                            m2lambdav(i) = rbuf_M(i,glbidx)
                        end do
                    end if
            ! ------------------------------------
                case('MISRZ')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    misrz = rbuf_B(glbidx)
            ! ------------------------------------
                case('M2ISRZ')
                    call pmf_accu_read_rbuf_B(cstaccu,iounit,keyline,rbuf_B)
                    m2isrz = rbuf_B(glbidx)
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
! Subroutine:  cst_restart_write_data
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

    ibuf_B(:) = 0
    ibuf_B(glbidx) = faccumulation
    call pmf_accu_write_ibuf_B(cstaccu,iounit,'NSAMPLES','AD',ibuf_B)

! ------------------------------------------------
    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = mlambda(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,'MLAMBDA','WA',rbuf_M)

    rbuf_M(:,:) = 0.0d0
    do i=1,cstaccu%tot_cvs
        rbuf_M(i,glbidx) = m2lambda(i)
    end do
    call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2LAMBDA','AD',rbuf_M)

    if( has_lambdav ) then
        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = mlambdav(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'MLAMBDAV','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m2lambdav(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M2LAMBDAV','AD',rbuf_M)
    end if

! ------------------------------------------------
    rbuf_B(:) = 0.0d0
    rbuf_B(glbidx) = misrz
    call pmf_accu_write_rbuf_B(cstaccu,iounit,'MISRZ','WA',rbuf_B)

    rbuf_B(:) = 0.0d0
    rbuf_B(glbidx) = m2isrz
    call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2ISRZ','AD',rbuf_B)

! ------------------------------------------------
    if( fenthalpy .or. fentropy ) then
        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = metot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'METOT','WA',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2etot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2ETOT','AD',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = mepot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MEPOT','WA',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2epot
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2EPOT','AD',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = mekin
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MEKIN','WA',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2ekin
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2EKIN','AD',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = merst
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'MERST','WA',rbuf_B)

        rbuf_B(:) = 0.0d0
        rbuf_B(glbidx) = m2erst
        call pmf_accu_write_rbuf_B(cstaccu,iounit,'M2ERST','AD',rbuf_B)
    end if

! ------------------------------------------------
    if ( fentropy )  then
        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hh(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HH','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HP','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hk(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HK','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c11hr(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C11HR','WA',rbuf_M)

    ! errors
        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c22hh(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C22HH','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c22hp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C22HP','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c22hk(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C22HK','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = c22hr(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'C22HR','WA',rbuf_M)
    end if




end subroutine cst_accu_write

!===============================================================================

end module cst_accu
