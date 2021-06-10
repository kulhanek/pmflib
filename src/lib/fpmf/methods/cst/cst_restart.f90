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

module cst_restart

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_restart_read
!===============================================================================

subroutine cst_restart_read

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                         :: i,it,fstat,nsize,nitems,idum
    character(len=PMF_MAX_TYPE)     :: stype,sunit
    character(len=PMF_MAX_KEY)      :: key
    character                       :: flag
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: fconv,rdum
    integer                         :: ibitem
    character(len=4)                :: scst
    character(len=2)                :: sver
    character(len=80)               :: line
    ! --------------------------------------------------------------------------

    ! read restart - only master read file --------------------------------
    if( .not. frestart ) then
        write(CST_OUT,'(A)') '# RST: frestart = OFF'
        return
    end if

    ! check if restart exist
    open(CST_RST,FILE=fcstrst,STATUS='OLD',ERR=1000)
    close(CST_RST)

    ! open restart
    call pmf_utils_open(CST_RST,fcstrst,'O')

    write(CST_OUT,'(A)') '# RST: frestart = ON'

    ! read header --------------------------
    read(CST_RST,*,end=50,err=50) scst, sver, ibitem

    if( trim(adjustl(scst)) .ne. 'CST' ) then
        call pmf_utils_exit(PMF_OUT,1,'[CST] Missing CST key in CST restart!')
    end if

    nitems = NumOfCONs - NumOfSHAKECONs
    if( ibitem .ne. nitems ) then
        call pmf_utils_exit(PMF_OUT,1,'[CST] CST restart contains different number of CVs!')
    end if

    if( trim(adjustl(sver)) .ne. 'V6' ) then
        call pmf_utils_exit(PMF_OUT,1,'[CST] Unsupported CST restart version: '''//trim(adjustl(sver))//'''!')
    end if

    do while(.true.)

    ! read key
        key     = 'ERROR'
        flag    = 'N'
        nsize   = 0

        read(CST_RST,'(A80)',end=500,err=500) line

        read(line,5,iostat=fstat) key, flag, nsize
        if( fstat .ne. 0 ) then
            read(line,5,err=51) key
        end if

        write(*,*) key, flag, nsize

        select case( trim(key) )
        ! ---------------------------------------
            case('CVS')
                do i=1, nitems
                    ! read CV definition
                    read(CST_RST,10,end=55,err=55) it, stype, sname
                    ! check CV definition
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[CST] Incorrect item in CST restart!')
                    end if
                    if( trim(adjustl(stype)) .ne. trim(CONList(i)%cv%ctype) ) then
                        write(PMF_OUT,*) '[CST] CV type = [',trim(adjustl(stype)),'] should be [',trim(CONList(i)%cv%ctype),']'
                        call pmf_utils_exit(PMF_OUT,1,'[CST] CV type was redefined in CST accumulator!')
                    end if
                    if( trim(adjustl(sname)) .ne. trim(CONList(i)%cv%name) ) then
                        write(PMF_OUT,*) '[CST] CV name = [',trim(adjustl(sname)),'] should be [',trim(CONList(i)%cv%name),']'
                        call pmf_utils_exit(PMF_OUT,1,'[CST] CV name was redefined in CST accumulator!')
                    end if

                    ! read unit but ignore them
                    read(CST_RST,11,end=56,err=56) it, fconv, sunit
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[CST] Incorrect item in CST accumulator!')
                    end if
                end do
        ! ---------------------------------------
            case('ENERGY-UNIT')
                ! ignore
                read(CST_RST,20,end=61,err=61) rdum, sunit
        ! ---------------------------------------
            case('TEMPERATURE')
                ! ignore
                read(CST_RST,30,end=62,err=62) rdum
        ! ---------------------------------------
            case('NSAMPLES')
                read(CST_RST,40,end=100,err=100) faccumulation, faccurst
        ! ---------------------------------------
            case('MISRZ')
                read(CST_RST,30,end=62,err=62) misrz
        ! ---------------------------------------
            case('M2ISRZ')
                read(CST_RST,30,end=62,err=62) m2isrz
        ! ---------------------------------------
            case('MLAMBDA')
                read(CST_RST,41,end=101,err=101) (mlambda(i),i=1,nitems)
        ! ---------------------------------------
            case('M2LAMBDA')
                read(CST_RST,41,end=102,err=102) (m2lambda(i),i=1,nitems)
                                        ! ---------------------------------------
            case('MLAMBDAV')
                read(CST_RST,41,end=103,err=103) (mlambdav(i),i=1,nitems)
        ! ---------------------------------------
            case('M2LAMBDAV')
                read(CST_RST,41,end=104,err=104) (m2lambdav(i),i=1,nitems)
        ! ---------------------------------------
            case('NSAMPLES-ENTHALPY')
                read(CST_RST,40,end=110,err=110) fentaccu
        ! ---------------------------------------
            case('METOT')
                read(CST_RST,41,end=111,err=111) metot
        ! ---------------------------------------
            case('M2ETOT')
                read(CST_RST,41,end=112,err=112) m2etot
        ! ---------------------------------------
            case('MEPOT')
                read(CST_RST,41,end=113,err=113) mepot
        ! ---------------------------------------
            case('M2EPOT')
                read(CST_RST,41,end=114,err=114) m2epot
        ! ---------------------------------------
            case('MEKIN')
                read(CST_RST,41,end=115,err=115) mekin
        ! ---------------------------------------
            case('M2EKIN')
                read(CST_RST,41,end=116,err=116) m2ekin
        ! ---------------------------------------
            case('NSAMPLES-ENTROPY')
                read(CST_RST,40,end=120,err=120) fentaccu
        ! ---------------------------------------
            case('CDS-HP')
                if( fentropy ) then
                    read(CST_RST,41,end=121,err=121) (cds_hp(i),i=1,nitems)
                else
                    read(CST_RST,41,end=121,err=121) (rdum,i=1,nitems)
                end if
        ! ---------------------------------------
            case('CDS-HK')
                if( fentropy ) then
                    read(CST_RST,41,end=122,err=122) (cds_hk(i),i=1,nitems)
                else
                    read(CST_RST,41,end=122,err=122) (rdum,i=1,nitems)
                end if

            case default
                select case(flag)
                    case('I')
                        write(CST_OUT,7) trim(key),flag,nsize
                        read(CST_RST,40,end=71,err=71) (idum,i=1,nsize)
                    case('R')
                        write(CST_OUT,7) trim(key),flag,nsize
                        read(CST_RST,41,end=72,err=72) (rdum,i=1,nsize)
                    case default
                        call pmf_utils_exit(PMF_OUT,1, &
                             '[CST] Unable to read from CST restart - unrecognized keyword and flag: '//trim(key)//'/'//flag)
                end select
        end select
    end do

 500 close(CST_RST)
     return

1000 write(CST_OUT,'(A)') '# WARNING: Unable to load restart information but frestart = on'
     frestart = .false.
     return

 5  format(A20,1X,A1,1X,I9)
 7  format('# CST restart - section ignored: ',A30,1X,A1,1X,I9)
10  format(I2,1X,A10,1X,A44)
11  format(I2,1X,E18.11,1X,A36)
20  format(E18.11,1X,A36)
30  format(F10.4)
40  format(8(I9,1X))
41  format(4(E19.11,1X))

 50 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - illegal header!')
 51 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - illegal keyword!')

 55 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - CV section - def/name!')
 56 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - CV section - unit!')

 61 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - energy unit!')
 62 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - temperature!')

 71 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - ignored section - I!')
 72 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - ignored section - R!')

100 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - NSAMPLES!')
101 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDA!')
102 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDA!')
103 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDAV!')
104 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MLAMBDAV!')

110 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - NSAMPLES-ENTHALPY!')
111 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - METOT!')
112 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - M2ETOT!')
113 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MEPOT!')
114 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - M2EPOT!')
115 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - MEKIN!')
116 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - M2EKIN!')

120 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - NSAMPLES-ENTROPY!')
121 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - CDS-HP!')
122 call pmf_utils_exit(PMF_OUT,1,'[CST] Unable to read from CST restart - data section - CDS-HK!')

end subroutine cst_restart_read

!===============================================================================
! Subroutine:  cst_restart_update
!===============================================================================

subroutine cst_restart_update

    use pmf_dat
    use pmf_utils
    use cst_dat

    implicit none
    !---------------------------------------------------------------------------

    if( frstupdate .le. 0 ) return ! restart is updated only of frstupdate > 0

    if( mod(fstep,frstupdate) .ne. 0 ) return

    call pmf_utils_open(CST_RST,fcstrst,'U')
    call cst_restart_write_data(CST_RST,.false.)
    close(CST_RST)

    write(CST_OUT,10) fstep

    return

 10 format('# [ACCU] Total steps     = ',I12)

end subroutine cst_restart_update

!===============================================================================
! Subroutine:  cst_restart_write
!===============================================================================

subroutine cst_restart_write

    use pmf_dat
    use pmf_utils

    implicit none
    !---------------------------------------------------------------------------

    call pmf_utils_open(CST_RST,fcstrst,'U')
    call cst_restart_write_data(CST_RST,.false.)
    close(CST_RST)

    call pmf_utils_open(CST_RST,fcstfrst,'U')
    call cst_restart_write_data(CST_RST,.true.)
    close(CST_RST)

    return

end subroutine cst_restart_write

!===============================================================================
! Subroutine:  cst_restart_write_data
!===============================================================================

subroutine cst_restart_write_data(iounit,full)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer                     :: iounit
    logical                     :: full
    ! -----------------------------------------------
    integer                     :: i,nitems
    character(len=PMF_MAX_KEY)  :: key
    ! --------------------------------------------------------------------------

! write header --------------------------
    nitems = NumOfCONs
    if( .not. full ) then
        nitems = NumOfCONs - NumOfSHAKECONs
    end if

    write(iounit,10) 'CST ', 'V6', nitems

    key = 'CVS'
    write(iounit,5) adjustl(key)
    do i=1, nitems
        write(iounit,25) i,trim(CONList(i)%cv%ctype), trim(CONList(i)%cv%name)
        write(iounit,26) i,pmf_unit_get_rvalue(CONList(i)%cv%unit,1.0d0),trim(pmf_unit_label(CONList(i)%cv%unit))
    end do

    key = 'TEMPERATURE'
    write(iounit,5) adjustl(key)
    write(iounit,6) ftemp

    key = 'ENERGY-UNIT'
    write(iounit,5) adjustl(key)
    write(iounit,27) pmf_unit_get_rvalue(EnergyUnit,1.0d0),trim(pmf_unit_label(EnergyUnit))

    key = 'NSAMPLES'
    write(iounit,5) adjustl(key), 'I', 2
    write(iounit,30) faccumulation,faccurst

    key = 'CVSVALUES'
    write(iounit,5) adjustl(key), 'R', nitems
    write(iounit,40) (CONList(i)%value, i=1,nitems)

    key = 'MISRZ'
    write(iounit,5) adjustl(key), 'R', 1
    write(iounit,40) misrz

    key = 'M2ISRZ'
    write(iounit,5) adjustl(key), 'R', 1
    write(iounit,40) m2isrz

    key = 'MLAMBDA'
    write(iounit,5) adjustl(key), 'R', nitems
    write(iounit,40) (mlambda(i), i=1,nitems)

    key = 'M2LAMBDA'
    write(iounit,5) adjustl(key), 'R', nitems
    write(iounit,40) (m2lambda(i), i=1,nitems)

    if( has_lambdav ) then
        key = 'MLAMBDAV'
        write(iounit,5) adjustl(key), 'R', nitems
        write(iounit,40) (mlambdav(i), i=1,nitems)

        key = 'M2LAMBDAV'
        write(iounit,5) adjustl(key), 'R', nitems
        write(iounit,40) (m2lambdav(i), i=1,nitems)
    end if

    if( fenthalpy ) then
        key = 'NSAMPLES-ENTHALPY'
        write(iounit,5) adjustl(key), 'I', 1
        write(iounit,30) fentaccu

        key = 'METOT'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) metot

        key = 'M2ETOT'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) m2etot

        key = 'MEPOT'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) mepot

        key = 'M2EPOT'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) m2epot

        key = 'MEKIN'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) mekin

        key = 'M2EKIN'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) m2ekin

        key = 'MERST'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) merst

        key = 'M2ERST'
        write(iounit,5) adjustl(key), 'R', 1
        write(iounit,40) m2erst
    end if

    if ( fentropy )  then
        key = 'NSAMPLES-ENTROPY'
        write(iounit,5) adjustl(key), 'I', 1
        write(iounit,30) fentaccu

        key = 'CDS-HP'
        write(iounit,5) adjustl(key), 'R', nitems
        write(iounit,40) (cds_hp(i), i=1,nitems)

        key = 'CDS-HK'
        write(iounit,5) adjustl(key), 'R', nitems
        write(iounit,40) (cds_hk(i), i=1,nitems)

        key = 'CDS-HR'
        write(iounit,5) adjustl(key), 'R', nitems
        write(iounit,40) (cds_hr(i), i=1,nitems)
    end if

    return

 5  format(A20,1X,A1,1X,I9)
 6  format(F10.4)
10  format(A4,1X,A2,1X,I2)
25  format(I2,1X,A10,1X,A44)
26  format(I2,1X,E18.11,1X,A36)
27  format(3X,E18.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

end subroutine cst_restart_write_data

!===============================================================================
! Subroutine:  cst_restart_trajectory_open
!===============================================================================

subroutine cst_restart_trajectory_open

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only if ftrjsample > 0

    call pmf_utils_open(CST_TRJ,fabftrj,'R')

    write(CST_TRJ,10)

    return

    10 format('# CSTTRAJ')

end subroutine cst_restart_trajectory_open

!===============================================================================
! Subroutine:  cst_restart_trajectory_write_snapshot
!===============================================================================

subroutine cst_restart_trajectory_write_snapshot

    use pmf_utils
    use pmf_dat
    use pmf_constants
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    if( mod(fstep,ftrjsample) .ne. 0 ) return

    ! write time
    write(CST_TRJ,10) fstep

    ! write accumulator
    call cst_restart_write_data(CST_TRJ,.false.)

    return

10 format('# CSTSNAP',I7)

end subroutine cst_restart_trajectory_write_snapshot

!===============================================================================
! Subroutine:  cst_restart_trajectory_close
!===============================================================================

subroutine cst_restart_trajectory_close

    use pmf_constants
    use pmf_dat
    use cst_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( ftrjsample .le. 0 ) return ! trajectory is written only of ftrjsample > 0

    close(CST_TRJ)

    return

end subroutine cst_restart_trajectory_close

!===============================================================================

end module cst_restart

