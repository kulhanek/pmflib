!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module pmf_accu

use pmf_constants
use pmf_cvs

implicit none

! accu CVs ---------------------------------------------------------------------
type CVInfoTypeACCU
    integer                 :: nbins            ! number of accumulator bins
    real(PMFDP)             :: min_value        ! left boundary of coordinate
    real(PMFDP)             :: max_value        ! left boundary of coordinate
    real(PMFDP)             :: bin_width        ! (right-left)/numbins
    real(PMFDP)             :: width            ! right - left
    class(CVType),pointer   :: cv               ! cv data
end type CVInfoTypeACCU

! common accu part -------------------------------------------------------------

type PMFAccuType
    character(PMF_MAX_TYPE)         :: method       ! PMF method
    integer                         :: tot_cvs      ! total number of independent CVs
    type(CVInfoTypeACCU), pointer   :: sizes(:)     ! CV information
    integer                         :: tot_nbins    ! number of total bins

end type PMFAccuType

contains

!===============================================================================
! Function:  pmf_accu_index
! compute accu index for one CV coordinate
! Arguments:
!               idxcoord    ... index of CV
!               accuvalue   ... CV value
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function pmf_accu_index(accu,idxcoord,accuvalue)

    use pmf_dat

    implicit none
    type(PMFAccuType)   :: accu
    integer             :: idxcoord
    real(PMFDP)         :: accuvalue
    ! --------------------------------------------------------------------------

    ! we need number from zero - therefore we use floor(x)
    pmf_accu_index = floor((accuvalue - accu%sizes(idxcoord)%min_value) / &
                                   accu%sizes(idxcoord)%bin_width)

    if( pmf_accu_index .lt. 0 .or. pmf_accu_index .ge.  accu%sizes(idxcoord)%nbins) then
        pmf_accu_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function pmf_accu_index

!===============================================================================
! Function:  pmf_accu_globalindex
! Description:  Compute globalindex for accumulator, based on accuvalues of all coordinates
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function pmf_accu_globalindex(accu,lvalues)

    use pmf_dat

    implicit none
    type(PMFAccuType)   :: accu
    real(PMFDP)         :: lvalues(:)
    ! -----------------------------------------------
    integer             :: idx_local,i
    ! --------------------------------------------------------------------------

    pmf_accu_globalindex = 0

    do i=1,accu%tot_cvs
        idx_local = pmf_accu_index(accu,i,lvalues(i))

        if (idx_local .eq. -1) then
            pmf_accu_globalindex = -1
            return
        end if

        pmf_accu_globalindex = pmf_accu_globalindex*accu%sizes(i)%nbins + idx_local
    end do

    pmf_accu_globalindex = pmf_accu_globalindex + 1

    return

end function pmf_accu_globalindex

!===============================================================================
! Subroutine:  pmf_accu_get_point
!===============================================================================

subroutine pmf_accu_get_point(accu,gidx,point)

    use pmf_dat

    implicit none
    type(PMFAccuType)   :: accu
    integer             :: gidx
    real(PMFDP)         :: point(:)
    ! --------------------------------------------
    integer             :: idx
    integer             :: k, cvbin
    ! --------------------------------------------------------------------------

    idx = gidx - 1

    do k=accu%tot_cvs,1,-1
        ! from 0 to nbins-1
        cvbin    = mod(idx,accu%sizes(k)%nbins)
        point(k) = accu%sizes(k)%min_value + (cvbin+0.5d0)*accu%sizes(k)%bin_width
        idx      = idx / accu%sizes(k)%nbins
    end do

end subroutine pmf_accu_get_point

!===============================================================================
! Subroutine:  pmf_accu_read_header
!===============================================================================

subroutine pmf_accu_read_header(accu,iounit,method,keyline)

    use pmf_dat
    use pmf_utils

    implicit none
    type(PMFAccuType)               :: accu
    integer                         :: iounit
    character(*)                    :: method
    character(*)                    :: keyline
    ! -----------------------------------------------
    integer                         :: i,it,nbins,ncvs,itmp
    character(len=PMF_MAX_TYPE)     :: stype,sunit
    character(len=PMF_MAX_CV_NAME)  :: sname
    character(len=PMF_MAX_PATH)     :: sbuff
    real(PMFDP)                     :: min_value,max_value,fconv
    ! --------------------------------------------------------------------------

    select case( pmf_accu_get_key(keyline) )
        case('PMFLIB-V6')
            ! read and check
            read(iounit,15,end=101,err=101) ncvs
            if( ncvs .ne. accu%tot_cvs ) then
                write(sbuff,110) ncvs,accu%tot_cvs
                call pmf_utils_exit(PMF_OUT,1,sbuff)
            end if
        case('VERSION')
            ! read but do not use
            read(iounit,10,end=102,err=102) sbuff
        case('SYSTYPE')
            ! read but do not use
            read(iounit,40,end=102,err=102) itmp
            if( itmp .ne. fsystype ) then
                call pmf_utils_exit(PMF_OUT,1, &
                                    '[PMFAccu] Inconsistent systype!')
            end if
        case('METHOD')
            ! read and check
            read(iounit,10,end=103,err=103) sbuff
            if( trim(sbuff) .ne. trim(method) ) then
                call pmf_utils_exit(PMF_OUT,1, &
                                    '[PMFAccu] Provided method "'//trim(sbuff)//'" is not requested one "'//trim(method)//'"!')
            end if
        case('DRIVER')
            ! read but do not use
            read(iounit,10,end=104,err=104) sbuff
        case('CVS')
            ! read header --------------------------
            do i=1, accu%tot_cvs
                ! read CV definition
                read(iounit,20,end=201,err=201) it, min_value, max_value, nbins, stype
                ! check CV definition
                if( it .ne. i ) then
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Incorrect item in the accumulator (CVS section)!')
                end if
                if( trim(adjustl(stype)) .ne. trim(accu%sizes(i)%cv%ctype) ) then
                    write(PMF_OUT,*) '[PMFAccu] CV type = [',trim(adjustl(stype)),'] should be [',trim(accu%sizes(i)%cv%ctype),']'
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] CV type was redefined in the accumulator (CVS section)!')
                end if
                if( abs(min_value-accu%sizes(i)%min_value) .gt. abs(accu%sizes(i)%min_value/100000.0d0) ) then
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Minimal value of CV was redefined in the accumulator (CVS section)!')
                end if
                if( abs(max_value-accu%sizes(i)%max_value) .gt. abs(accu%sizes(i)%max_value/100000.0d0) ) then
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Maximal value of CV was redefined in the accumulator (CVS section)!')
                end if
                if( nbins .ne. accu%sizes(i)%nbins ) then
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Number of CV bins was redefined in the accumulator (CVS section)!')
                end if

                ! read names
                read(iounit,25,end=202,err=202) it, sname
                ! check names
                if( it .ne. i ) then
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Incorrect item in the accumulator (CVS section)!')
                end if
                if( trim(adjustl(sname)) .ne. trim(accu%sizes(i)%cv%name) ) then
                    write(PMF_OUT,*) '[PMFAccu] CV name = [',trim(adjustl(sname)),'] should be [',trim(accu%sizes(i)%cv%name),']'
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] CV name was redefined in the accumulator (CVS section)!')
                end if

                ! read names
                read(iounit,26,end=203,err=203) it, fconv, sunit
                ! check names
                if( it .ne. i ) then
                    call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Incorrect item in the accumulator (CVS section)!')
                end if
                ! ignore values fconv and sunit
            end do
        case('ENERGY-UNIT')
            ! read but do not use
            read(iounit,27,end=301,err=301) fconv, sunit
        case('TEMPERATURE')
            ! read and test
            read(iounit,6,end=302,err=302) fconv
            if( (fconv - ftemp) .gt. 1.0d0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Inconsistent temperatures!')
            end if
        case('TEMPERATURE-UNIT')
            ! read but do not use
            read(iounit,27,end=303,err=303) fconv, sunit
        case('PRESSURE')
            ! read and test
            read(iounit,6,end=306,err=306) fconv
            if( (fconv - fpressure) .gt. 1.0d0 ) then
                call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Inconsistent pressures!')
            end if
        case('PRESSURE-UNIT')
            ! read but do not use
            read(iounit,27,end=307,err=307) fconv, sunit
        case('NSTLIMIT')
            ! read but do not use
            read(iounit,30,end=304,err=304) itmp
        case('TIMESTEP')
            ! read but do not use
            read(iounit,6,end=305,err=305) fconv
        case('CURRSTEP')
            ! read but do not use
            read(iounit,30,end=304,err=304) itmp
        case default
            call pmf_utils_exit(PMF_OUT,1, &
                 '[PMFAccu] Unable to read from the accumulator - unrecognized header keyword on keyline: "'//trim(keyline)//'"')
    end select

    return

 6  format(F10.4)

10  format(A)
15  format(I2)

20  format(I2,1X,E18.11,1X,E18.11,1X,I6,1X,A10)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
27  format(3X,E18.11,1X,A36)
30  format(I15)

40  format(I10)

101 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - number of CVS!')
102 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - version!')
103 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - method!')
104 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - driver!')

110 format('[PMFAccu] Incorrect number of CVs: ',I3,' but ',I3,' is expected!')

201 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - CVS section - def!')
202 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - CVS section - name!')
203 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - CVS section - unit!')

301 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - energy unit!')
302 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - temperature!')
303 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - temperature unit!')
304 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - number of MD steps!')
305 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - time step!')
306 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - pressure!')
307 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - pressure unit!')

end subroutine pmf_accu_read_header

!===============================================================================
! Subroutine:  pmf_accu_write_header
!===============================================================================

subroutine pmf_accu_write_header(accu,iounit)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    ! -----------------------------------------------
    integer                     :: i
    character(len=PMF_MAX_KEY)  :: key
    !---------------------------------------------------------------------------

    ! write header --------------------------
    key = '%PMFLIB-V6'
    write(iounit,5) adjustl(key)
    write(iounit,15) accu%tot_cvs

    key = '%VERSION'
    write(iounit,5) adjustl(key)
    write(iounit,10) trim(PMFLIBVER)

    key = '%METHOD'
    write(iounit,5) adjustl(key)
    write(iounit,10) trim(accu%method)

    key = '%DRIVER'
    write(iounit,5) adjustl(key)
    write(iounit,10) trim(DriverName)

    key = '%TEMPERATURE'
    write(iounit,5) adjustl(key)
    write(iounit,20) ftemp

    key = '%SYSTYPE'
    write(iounit,5) adjustl(key)
    write(iounit,26) fsystype

    if( fsystype .eq. 2 ) then
    key = '%PRESSURE'
    write(iounit,5) adjustl(key)
    write(iounit,25) fpressure
    end if

    key = '%CVS'
    write(iounit,5) adjustl(key)
    do i=1, accu%tot_cvs
        write(iounit,30) i, &
                          accu%sizes(i)%min_value,accu%sizes(i)%max_value, &
                          accu%sizes(i)%nbins, &
                          trim(accu%sizes(i)%cv%ctype)
        write(iounit,31) i,trim(accu%sizes(i)%cv%name)
        write(iounit,32) i,pmf_unit_get_rvalue(accu%sizes(i)%cv%unit,1.0d0),trim(pmf_unit_label(accu%sizes(i)%cv%unit))
    end do

    key = '%ENERGY-UNIT'
    write(iounit,5) adjustl(key)
    write(iounit,40) pmf_unit_get_rvalue(EnergyUnit,1.0d0),trim(pmf_unit_label(EnergyUnit))

    key = '%TEMPERATURE-UNIT'
    write(iounit,5) adjustl(key)
    write(iounit,40) pmf_unit_get_rvalue(TemperatureUnit,1.0d0),trim(pmf_unit_label(TemperatureUnit))

    if( fsystype .eq. 2 ) then
    key = '%PRESSURE-UNIT'
    write(iounit,5) adjustl(key)
    write(iounit,40) pmf_unit_get_rvalue(PressureUnit,1.0d0),trim(pmf_unit_label(PressureUnit))
    end if

    key = '%NSTLIMIT'
    write(iounit,5) adjustl(key)
    write(iounit,41) fnstlim

    key = '%TIMESTEP'
    write(iounit,5) adjustl(key)
    write(iounit,20) fdt

    key = '%CURRSTEP'
    write(iounit,5) adjustl(key)
    write(iounit,41) fstep

    return

 5  format(A20)
10  format(A)
15  format(I2)
20  format(F10.4)
25  format(F10.1)
26  format(I10)

30  format(I2,1X,E18.11,1X,E18.11,1X,I6,1X,A10)
31  format(I2,1X,A55)
32  format(I2,1X,E18.11,1X,A36)

40  format(3X,E18.11,1X,A36)
41  format(I15)

end subroutine pmf_accu_write_header

!===============================================================================
! Function:  pmf_accu_is_header_key
!===============================================================================

function pmf_accu_is_header_key(keyline) result(rst)

    use pmf_dat
    use pmf_utils

    implicit none
    character(*)                :: keyline
    logical                     :: rst
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: key
    !---------------------------------------------------------------------------

    rst = .false.
    read(keyline,*,end=200,err=200) key

    if( len(key) .ge. 1 ) then
        ! header keys start with '%'
        if( key(1:1) .eq. '%' ) then
            rst = .true.
        end if
        ! while data keys start with '@'
    end if

    return

200 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_is_header_key!')

end function pmf_accu_is_header_key

!===============================================================================
! Function:  pmf_accu_get_key
!===============================================================================

function pmf_accu_get_key(keyline) result(key)

    use pmf_dat
    use pmf_utils

    implicit none
    character(*)                :: keyline
    character(len=PMF_MAX_KEY)  :: key
    !---------------------------------------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    !---------------------------------------------------------------------------

    key = ''
    read(keyline,*,end=200,err=200) skey

    if( len(skey) .gt. 1 ) then
        key = skey(2:len(skey))
    end if

    return

200 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_get_key!')

end function pmf_accu_get_key

!===============================================================================
! Subroutine:  pmf_accu_skip_section
!===============================================================================

subroutine pmf_accu_skip_section(iounit,keyline,notify)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    integer                     :: iounit
    character(*)                :: keyline
    integer                     :: notify
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    character(len=PMF_MAX_KEY)  :: dop
    character(len=PMF_MAX_KEY)  :: dtype
    character(len=PMF_MAX_KEY)  :: dmode
    integer                     :: dsize
    integer                     :: i
    integer                     :: ibuf
    real(PMFDP)                 :: rbuf
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dop,dtype,dmode,dsize

    write(notify,50) trim(skey),dtype,dsize

    select case(dtype)
        case('I')
            read(iounit,10,end=110,err=110) (ibuf,i=1,dsize)
        case('R')
            read(iounit,20,end=110,err=110) (rbuf,i=1,dsize)
        case default
            call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require I or R but got "'//dtype//'" in pmf_accu_skip_section!')
    end select

    return

10  format(8(I9,1X))
20  format(4(ES23.15E3,1X))

50  format('# [PMFAccu] Ignoring section [',A,'] of type "',A1,'" with length: ',I10)

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_skip_section!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_skip_section!')

end subroutine pmf_accu_skip_section

!===============================================================================
! Subroutine:  pmf_accu_read_ibuf_B
! per bin
!===============================================================================

subroutine pmf_accu_read_ibuf_B(accu,iounit,keyline,ibuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: keyline
    integer                     :: ibuf(:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    character(len=PMF_MAX_KEY)  :: dop
    character(len=PMF_MAX_KEY)  :: dtype
    character(len=PMF_MAX_KEY)  :: dmode
    integer                     :: dsize
    integer                     :: i
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dop,dtype,dmode,dsize

    if( dtype .ne. 'I' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require I but got "'//trim(dtype)//'" in pmf_accu_read_ibuf_B!')
    end if
    if( dmode .ne. 'B' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Mode mismatch: require "B" but got "'//trim(dmode)//'" in pmf_accu_read_ibuf_B!')
    end if

    if( dsize .ne. accu%tot_nbins ) then
        write(serr,120) accu%tot_nbins, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    read(iounit,10,end=110,err=110) (ibuf(i),i=1,accu%tot_nbins)

    return

10  format(8(I9,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_ibuf_B!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_ibuf_B!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_ibuf_B!')

end subroutine pmf_accu_read_ibuf_B

!===============================================================================
! Subroutine:  pmf_accu_write_ibuf
!===============================================================================

subroutine pmf_accu_write_ibuf_B(accu,iounit,key,op,ibuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    character(*)                :: op
    integer                     :: ibuf(:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), trim(op), 'I', 'B', accu%tot_nbins
    write(iounit,10) (ibuf(i),i=1,accu%tot_nbins)

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(8(I9,1X))

end subroutine pmf_accu_write_ibuf_B

!===============================================================================
! Subroutine:  pmf_accu_read_ibuf_B
! per bin
!===============================================================================

subroutine pmf_accu_read_ibuf_D(iounit,keyline,ibuf,nsize)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    integer                     :: iounit
    character(*)                :: keyline
    integer                     :: ibuf(:)
    integer                     :: nsize
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    character(len=PMF_MAX_KEY)  :: dop
    character(len=PMF_MAX_KEY)  :: dtype
    character(len=PMF_MAX_KEY)  :: dmode
    integer                     :: dsize
    integer                     :: i
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dop,dtype,dmode,dsize

    if( dtype .ne. 'I' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require I but got "'//trim(dtype)//'" in pmf_accu_read_ibuf_D!')
    end if
    if( dmode .ne. 'D' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Mode mismatch: require "D" but got "'//trim(dmode)//'" in pmf_accu_read_ibuf_D!')
    end if

    if( dsize .ne. nsize ) then
        write(serr,120) nsize, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    read(iounit,10,end=110,err=110) (ibuf(i),i=1,nsize)

    return

10  format(8(I9,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_ibuf_D!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_ibuf_D!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_ibuf_D!')

end subroutine pmf_accu_read_ibuf_D

!===============================================================================
! Subroutine:  pmf_accu_write_ibuf_D
!===============================================================================

subroutine pmf_accu_write_ibuf_D(iounit,key,op,ibuf,nsize)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    integer                     :: iounit
    character(*)                :: key
    character(*)                :: op
    integer                     :: ibuf(:)
    integer                     :: nsize
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), trim(op), 'I', 'D', nsize
    write(iounit,10) (ibuf(i),i=1,nsize)

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(8(I9,1X))

end subroutine pmf_accu_write_ibuf_D

!===============================================================================
! Subroutine:  pmf_accu_read_rbuf_B
!===============================================================================

subroutine pmf_accu_read_rbuf_B(accu,iounit,keyline,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: keyline
    real(PMFDP)                 :: rbuf(:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    character(len=PMF_MAX_KEY)  :: dop
    character(len=PMF_MAX_KEY)  :: dtype
    character(len=PMF_MAX_KEY)  :: dmode
    integer                     :: dsize
    integer                     :: i
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dop,dtype,dmode,dsize

    if( dtype .ne. 'R' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require "R" but got "'//trim(dtype)//'" in pmf_accu_read_rbuf_B!')
    end if
    if( dmode .ne. 'B' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Mode mismatch: require "B" but got "'//trim(dmode)//'" in pmf_accu_read_rbuf_B!')
    end if

    if( dsize .ne. accu%tot_nbins ) then
        write(serr,120) accu%tot_nbins, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    read(iounit,10,end=110,err=110) (rbuf(i),i=1,accu%tot_nbins)

    return

10  format(4(ES23.15E3,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_rbuf_B!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_rbuf_B!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_rbuf_B')

end subroutine pmf_accu_read_rbuf_B

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_B
! data per bin
!===============================================================================

subroutine pmf_accu_write_rbuf_B(accu,iounit,key,op,rbuf,xmean,ymean)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    character(*)                :: op
    real(PMFDP)                 :: rbuf(:)
    character(*),optional       :: xmean
    character(*),optional       :: ymean
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5,ADVANCE='NO') adjustl(skey), trim(op), 'R', 'B', accu%tot_nbins

    if( present(xmean) ) then
        skey = xmean
        write(iounit,6,ADVANCE='NO') adjustl(skey)
    end if
    if( present(ymean) ) then
        skey = ymean
        write(iounit,6,ADVANCE='NO') adjustl(skey)
    end if
    write(iounit,*)

    write(iounit,10) (rbuf(i),i=1,accu%tot_nbins)

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
 6  format(1X,A20)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_B

!===============================================================================
! Subroutine:  pmf_accu_read_rbuf_M
! mixed - data per bin and CV
!===============================================================================

subroutine pmf_accu_read_rbuf_M(accu,iounit,keyline,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: keyline
    real(PMFDP)                 :: rbuf(:,:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    character(len=PMF_MAX_KEY)  :: dop
    character(len=PMF_MAX_KEY)  :: dtype
    character(len=PMF_MAX_KEY)  :: dmode
    integer                     :: dsize
    integer                     :: i,j
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dop,dtype,dmode,dsize

    if( dtype .ne. 'R' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require "R" but got "'//dtype//'" in pmf_accu_read_rbuf_M!')
    end if
    if( dmode .ne. 'M' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Mode mismatch: require "M" but got "'//trim(dmode)//'" in pmf_accu_read_rbuf_M!')
    end if

    if( dsize .ne. accu%tot_nbins*accu%tot_cvs ) then
        write(serr,120) accu%tot_nbins*accu%tot_cvs, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    do i=1,accu%tot_cvs
        read(iounit,10,end=110,err=110) (rbuf(i,j),j=1,accu%tot_nbins)
    end do

    return

10  format(4(ES23.15E3,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_rbuf_M!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_rbuf_M!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_rbuf_M!')

end subroutine pmf_accu_read_rbuf_M

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_M
!===============================================================================

subroutine pmf_accu_write_rbuf_M(accu,iounit,key,op,rbuf,nsamples,xmean,ymean)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    character(*)                :: op
    real(PMFDP)                 :: rbuf(:,:)
    character(*),optional       :: nsamples
    character(*),optional       :: xmean
    character(*),optional       :: ymean
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i,j
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5,ADVANCE='NO') adjustl(skey), trim(op), 'R', 'M', accu%tot_nbins*accu%tot_cvs

    if( present(nsamples) ) then
        skey = nsamples
        write(iounit,6,ADVANCE='NO') adjustl(skey)
    end if
    if( present(xmean) ) then
        skey = xmean
        write(iounit,6,ADVANCE='NO') adjustl(skey)
    end if
    if( present(ymean) ) then
        skey = ymean
        write(iounit,6,ADVANCE='NO') adjustl(skey)
    end if
    write(iounit,*)

    do i=1,accu%tot_cvs
        write(iounit,10) (rbuf(i,j),j=1,accu%tot_nbins)
    end do

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
 6  format(1X,A20)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_M

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_C
! data per CVs
!===============================================================================

subroutine pmf_accu_read_rbuf_C(accu,iounit,keyline,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: keyline
    real(PMFDP)                 :: rbuf(:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    character(len=PMF_MAX_KEY)  :: dop
    character(len=PMF_MAX_KEY)  :: dtype
    character(len=PMF_MAX_KEY)  :: dmode
    integer                     :: dsize
    integer                     :: i
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dop,dtype,dmode,dsize

    if( dtype .ne. 'R' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require "R" but got "'//dtype//'" in pmf_accu_read_rbuf_C!')
    end if
    if( dmode .ne. 'C' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Mode mismatch: require "C" but got "'//trim(dmode)//'" in pmf_accu_read_rbuf_C!')
    end if

    if( dsize .ne. accu%tot_cvs ) then
        write(serr,120) accu%tot_cvs, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    read(iounit,10,end=110,err=110) (rbuf(i),i=1,accu%tot_cvs)

    return

10  format(4(ES23.15E3,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_rbuf_C!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_rbuf_C!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_rbuf_C!')

end subroutine pmf_accu_read_rbuf_C

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_C
! data per CVs
!===============================================================================

subroutine pmf_accu_write_rbuf_C(accu,iounit,key,op,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    character(*)                :: op
    real(PMFDP)                 :: rbuf(:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), trim(op), 'R', 'C', accu%tot_cvs
    write(iounit,10) (rbuf(i),i=1,accu%tot_cvs)

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_C

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_D
! data per CVs
!===============================================================================

subroutine pmf_accu_write_rbuf_D(iounit,key,op,rbuf,nsize)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    integer                     :: iounit
    character(*)                :: key
    character(*)                :: op
    real(PMFDP)                 :: rbuf(:)
    integer                     :: nsize
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), trim(op), 'R', 'D', nsize
    write(iounit,10) (rbuf(i),i=1,nsize)

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_D

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_T
! data per time
!===============================================================================

subroutine pmf_accu_write_rbuf_T(iounit,key,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    integer                     :: iounit
    character(*)                :: key
    real(PMFDP)                 :: rbuf(:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), 'IG', 'R', 'T', fnstlim
    write(iounit,10) (rbuf(i),i=1,fnstlim)

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_T

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_TC
! data per time and cvs
!===============================================================================

subroutine pmf_accu_write_rbuf_TC(accu,iounit,key,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    real(PMFDP)                 :: rbuf(:,:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i,j
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), 'IG', 'R', 'S', fnstlim
    do i=1,accu%tot_cvs
        write(iounit,10) (rbuf(i,j),j=1,fnstlim)
    end do

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_TC

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf_TC
! data per time and cvs x cvs
!===============================================================================

subroutine pmf_accu_write_rbuf_TCC(accu,iounit,key,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    real(PMFDP)                 :: rbuf(:,:,:)
    ! -----------------------------------------------
    character(len=PMF_MAX_KEY)  :: skey
    integer                     :: i,j,t
    !---------------------------------------------------------------------------

    skey = key
    write(iounit,5) adjustl(skey), 'IG', 'R', 'Z', fnstlim
    do i=1,accu%tot_cvs
        do j=1,accu%tot_cvs
            write(iounit,10) (rbuf(i,j,t),t=1,fnstlim)
        end do
    end do

 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
10  format(4(ES23.15E3,1X))

end subroutine pmf_accu_write_rbuf_TCC

!===============================================================================

end module pmf_accu

