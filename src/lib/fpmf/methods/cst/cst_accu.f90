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
    ! --------------------------------------------------------------------------

    do while(.true.)

        ! read keyline
        read(iounit,5,end=500,err=300) keyline

        ! process keyline
        if( pmf_accu_is_header_key(keyline) ) then
            call pmf_accu_read_header(cstaccu,iounit,'CST',keyline)
        else
!            select case( pmf_accu_get_key(keyline) )
!            ! ------------------------------------
!                case('NSAMPLES')
!                    call pmf_accu_read_ibuf_B(mtdaccu%PMFAccuType,iounit,keyline,mtdaccu%nsamples)
!            ! ------------------------------------
!                case('MTDPOT')
!                    call pmf_accu_read_rbuf_B(mtdaccu%PMFAccuType,iounit,keyline,mtdaccu%mtdpotential)
!            ! ------------------------------------
!                case('MTDFORCE')
!                    call pmf_accu_read_rbuf_M(mtdaccu%PMFAccuType,iounit,keyline,mtdaccu%mtdforce)
!            ! ------------------------------------
!                case default
!                    call pmf_accu_skip_section(iounit,keyline,MTD_OUT)
!            end select
        end if
    end do

    close(CST_RST)

500 return

  5 format(A80)

300 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from the accumulator - keyline!')


!    ! read restart - only master read file --------------------------------
!    if( .not. frestart ) then
!        write(CST_OUT,'(A)') '# RST: frestart = OFF'
!        return
!    end if
!
!    ! check if restart exist
!    open(CST_RST,FILE=fcstrst,STATUS='OLD',ERR=1000)
!    close(CST_RST)
!
!    ! open restart
!    call pmf_utils_open(CST_RST,fcstrst,'O')
!
!    write(CST_OUT,'(A)') '# RST: frestart = ON'
!
!    ! read header --------------------------
!    read(CST_RST,*,end=50,err=50) scst, sver, ibitem
!
!    if( trim(adjustl(scst)) .ne. 'CST' ) then
!        call pmf_utils_exit(PMF_OUT,1,'[CST] Missing CST key in CST restart!')
!    end if
!
!    nitems = NumOfCONs - NumOfSHAKECONs
!    if( ibitem .ne. nitems ) then
!        call pmf_utils_exit(PMF_OUT,1,'[CST] CST restart contains different number of CVs!')
!    end if
!
!    if( trim(adjustl(sver)) .ne. 'V6' ) then
!        call pmf_utils_exit(PMF_OUT,1,'[CST] Unsupported CST restart version: '''//trim(adjustl(sver))//'''!')
!    end if
!
!    do while(.true.)
!
!    ! read key
!        key     = 'ERROR'
!        flag    = 'N'
!        nsize   = 0
!
!        read(CST_RST,'(A80)',end=500,err=500) line
!
!        read(line,5,iostat=fstat) key, flag, nsize
!        if( fstat .ne. 0 ) then
!            read(line,5,err=51) key
!        end if
!
!        write(*,*) key, flag, nsize
!
!        select case( trim(key) )
!        ! ---------------------------------------
!            case('CVS')
!                do i=1, nitems
!                    ! read CV definition
!                    read(CST_RST,10,end=55,err=55) it, stype, sname
!                    ! check CV definition
!                    if( it .ne. i ) then
!                        call pmf_utils_exit(PMF_OUT,1,'[CST] Incorrect item in CST restart!')
!                    end if
!                    if( trim(adjustl(stype)) .ne. trim(CONList(i)%cv%ctype) ) then
!                        write(PMF_OUT,*) '[CST] CV type = [',trim(adjustl(stype)),'] should be [',trim(CONList(i)%cv%ctype),']'
!                        call pmf_utils_exit(PMF_OUT,1,'[CST] CV type was redefined in CST accumulator!')
!                    end if
!                    if( trim(adjustl(sname)) .ne. trim(CONList(i)%cv%name) ) then
!                        write(PMF_OUT,*) '[CST] CV name = [',trim(adjustl(sname)),'] should be [',trim(CONList(i)%cv%name),']'
!                        call pmf_utils_exit(PMF_OUT,1,'[CST] CV name was redefined in CST accumulator!')
!                    end if
!
!                    ! read unit but ignore them
!                    read(CST_RST,11,end=56,err=56) it, fconv, sunit
!                    if( it .ne. i ) then
!                        call pmf_utils_exit(PMF_OUT,1,'[CST] Incorrect item in CST accumulator!')
!                    end if
!                end do
!        ! ---------------------------------------
!            case('ENERGY-UNIT')
!                ! ignore
!                read(CST_RST,20,end=61,err=61) rdum, sunit
!        ! ---------------------------------------
!            case('TEMPERATURE')
!                ! ignore
!                read(CST_RST,30,end=62,err=62) rdum
!        ! ---------------------------------------
!            case('NSAMPLES')
!                read(CST_RST,40,end=100,err=100) faccumulation, faccurst
!        ! ---------------------------------------
!            case('MISRZ')
!                read(CST_RST,30,end=62,err=62) misrz
!        ! ---------------------------------------
!            case('M2ISRZ')
!                read(CST_RST,30,end=62,err=62) m2isrz
!        ! ---------------------------------------
!            case('MLAMBDA')
!                read(CST_RST,41,end=101,err=101) (mlambda(i),i=1,nitems)
!        ! ---------------------------------------
!            case('M2LAMBDA')
!                read(CST_RST,41,end=102,err=102) (m2lambda(i),i=1,nitems)
!                                        ! ---------------------------------------
!            case('MLAMBDAV')
!                read(CST_RST,41,end=103,err=103) (mlambdav(i),i=1,nitems)
!        ! ---------------------------------------
!            case('M2LAMBDAV')
!                read(CST_RST,41,end=104,err=104) (m2lambdav(i),i=1,nitems)
!        ! ---------------------------------------
!            case('NSAMPLES-ENTHALPY')
!                read(CST_RST,40,end=110,err=110) fentaccu
!        ! ---------------------------------------
!            case('METOT')
!                read(CST_RST,41,end=111,err=111) metot
!        ! ---------------------------------------
!            case('M2ETOT')
!                read(CST_RST,41,end=112,err=112) m2etot
!        ! ---------------------------------------
!            case('MEPOT')
!                read(CST_RST,41,end=113,err=113) mepot
!        ! ---------------------------------------
!            case('M2EPOT')
!                read(CST_RST,41,end=114,err=114) m2epot
!        ! ---------------------------------------
!            case('MEKIN')
!                read(CST_RST,41,end=115,err=115) mekin
!        ! ---------------------------------------
!            case('M2EKIN')
!                read(CST_RST,41,end=116,err=116) m2ekin
!        ! ---------------------------------------
!            case('NSAMPLES-ENTROPY')
!                read(CST_RST,40,end=120,err=120) fentaccu
!        ! ---------------------------------------
!!            case('M11HP')
!!                if( fentropy ) then
!!                    read(CST_RST,41,end=121,err=121) (cds_hp(i),i=1,nitems)
!!                else
!!                    read(CST_RST,41,end=121,err=121) (rdum,i=1,nitems)
!!                end if
!!        ! ---------------------------------------
!!            case('M11HK')
!!                if( fentropy ) then
!!                    read(CST_RST,41,end=122,err=122) (cds_hk(i),i=1,nitems)
!!                else
!!                    read(CST_RST,41,end=122,err=122) (rdum,i=1,nitems)
!!                end if
!
!            case default
!                select case(flag)
!                    case('I')
!                        write(CST_OUT,7) trim(key),flag,nsize
!                        read(CST_RST,40,end=71,err=71) (idum,i=1,nsize)
!                    case('R')
!                        write(CST_OUT,7) trim(key),flag,nsize
!                        read(CST_RST,41,end=72,err=72) (rdum,i=1,nsize)
!                    case default
!                        call pmf_utils_exit(PMF_OUT,1, &
!                             '[CST] Unable to read from CST restart - unrecognized keyword and flag: '//trim(key)//'/'//flag)
!                end select
!        end select
!    end do
!
! 500 close(CST_RST)
!     return
!
!1000 write(CST_OUT,'(A)') '# WARNING: Unable to load restart information but frestart = on'
!     frestart = .false.
!     return
!
! 5  format(A20,1X,A1,1X,I9)
! 7  format('# CST restart - section ignored: ',A30,1X,A1,1X,I9)
!10  format(I2,1X,A10,1X,A44)
!11  format(I2,1X,E18.11,1X,A36)
!20  format(E18.11,1X,A36)
!30  format(F10.4)
!40  format(8(I9,1X))
!41  format(4(E19.11,1X))
!
! 50 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - illegal header!')
! 51 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - illegal keyword!')
!
! 55 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - CV section - def/name!')
! 56 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - CV section - unit!')
!
! 61 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - energy unit!')
! 62 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - temperature!')
!
! 71 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - ignored section - I!')
! 72 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - ignored section - R!')
!
!100 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - NSAMPLES!')
!101 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDA!')
!102 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDA!')
!103 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDAV!')
!104 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDAV!')
!
!110 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - NSAMPLES-ENTHALPY!')
!111 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - METOT!')
!112 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - M2ETOT!')
!113 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MEPOT!')
!114 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - M2EPOT!')
!115 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MEKIN!')
!116 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - M2EKIN!')
!
!120 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - NSAMPLES-ENTROPY!')
!121 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - CDS-HP!')
!122 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - CDS-HK!')

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
    if( fenthalpy ) then
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
            rbuf_M(i,glbidx) = m11hp(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M11HP','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m11hk(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M11HK','WA',rbuf_M)

        rbuf_M(:,:) = 0.0d0
        do i=1,cstaccu%tot_cvs
            rbuf_M(i,glbidx) = m11hr(i)
        end do
        call pmf_accu_write_rbuf_M(cstaccu,iounit,'M11HR','WA',rbuf_M)
    end if

end subroutine cst_accu_write

!===============================================================================

end module cst_accu
