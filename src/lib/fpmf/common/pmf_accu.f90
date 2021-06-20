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
! Subroutine:  pmf_accu_read_header
!===============================================================================

subroutine pmf_accu_read_header(accu,iounit,keyline)

    use pmf_dat
    use pmf_utils

    implicit none
    type(PMFAccuType)               :: accu
    integer                         :: iounit
    character(*)                    :: keyline
    ! -----------------------------------------------
    integer                         :: i,it,nbins
    character(len=PMF_MAX_TYPE)     :: stype,sunit
    character(len=PMF_MAX_KEY)      :: key
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: min_value,max_value,fconv
    ! --------------------------------------------------------------------------

    select case( pmf_accu_get_key(keyline) )
        case('CVS')
            ! read header --------------------------
            do i=1, accu%tot_cvs
                ! read CV definition
                read(iounit,20,end=201,err=201) it, stype, min_value, max_value, nbins
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
        case('TEMP')
            ! read but do not use
            read(iounit,6,end=302,err=302) fconv
        case('ENERGY-UNIT')
            ! read but do not use
            read(iounit,27,end=301,err=301) fconv, sunit
        case default
            call pmf_utils_exit(PMF_OUT,1, &
                 '[PMFAccu] Unable to read from the accumulator - unrecognized header keyword: "'//trim(key)//'"')
    end select

    return

 6  format(F10.4)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
27  format(E18.11,1X,A36)

201 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - CVS section - def!')
202 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - CVS section - name!')
203 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - CVS section - unit!')

301 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - energy unit!')
302 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to read from the accumulator - temperature!')

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
    write(iounit,1) 'PMFLIB', 'V6', accu%tot_cvs

    key = '%VERSION'
    write(iounit,5) adjustl(key)
    write(iounit,10) trim(PMFLIBVER)

    key = '%METHOD'
    write(iounit,5) adjustl(key)
    write(iounit,10) trim(accu%method)

    key = '%TEMPERATURE'
    write(iounit,5) adjustl(key)
    write(iounit,20) ftemp

    key = '%CVS'
    write(iounit,5) adjustl(key)
    do i=1, accu%tot_cvs
        write(iounit,30) i,trim(accu%sizes(i)%cv%ctype), &
                          accu%sizes(i)%min_value,accu%sizes(i)%max_value, &
                          accu%sizes(i)%nbins
        write(iounit,31) i,trim(accu%sizes(i)%cv%name)
        write(iounit,32) i,pmf_unit_get_rvalue(accu%sizes(i)%cv%unit,1.0d0),trim(pmf_unit_label(accu%sizes(i)%cv%unit))
    end do

    key = '%ENERGY-UNIT'
    write(iounit,5) adjustl(key)
    write(iounit,40) pmf_unit_get_rvalue(EnergyUnit,1.0d0),trim(pmf_unit_label(EnergyUnit))

    return

 1  format(A6,1X,A2,1X,I2)
 5  format(A20)
10  format(A)
20  format(F10.4)

30  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
31  format(I2,1X,A55)
32  format(I2,1X,E18.11,1X,A36)

40  format(3X,E18.11,1X,A36)

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
    character(len=PMF_MAX_KEY)  :: dtype
    integer                     :: dsize
    integer                     :: i
    integer                     :: ibuf
    real(PMFDP)                 :: rbuf
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dtype,dsize
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
20  format(4(E19.11,1X))

50  format('[PMFAccu] Ignoring section [',A,'] of type "',A1,'" with length: ',I10)

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_skip_section!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_skip_section!')

end subroutine pmf_accu_skip_section

!===============================================================================
! Subroutine:  pmf_accu_read_ibuf
!===============================================================================

subroutine pmf_accu_read_ibuf(accu,iounit,keyline,ibuf)

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
    character(len=PMF_MAX_KEY)  :: dtype
    integer                     :: dsize
    integer                     :: i
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dtype,dsize

    if( dtype .ne. 'I' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require I but got "'//dtype//'" in pmf_accu_read_ibuf!')
    end if

    if( dsize .ne. accu%tot_nbins ) then
        write(serr,120) accu%tot_nbins, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    read(iounit,10,end=110,err=110) (ibuf(i),i=1,accu%tot_nbins)

    return

10  format(8(I9,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_ibuf!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_ibuf!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_ibuf!')

end subroutine pmf_accu_read_ibuf

!===============================================================================
! Subroutine:  pmf_accu_write_ibuf
!===============================================================================

subroutine pmf_accu_write_ibuf(accu,iounit,key,ibuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    integer                     :: ibuf(:)
    ! -----------------------------------------------
    integer                     :: i
    !---------------------------------------------------------------------------

    write(iounit,5) adjustl(key), 'I', accu%tot_nbins
    write(iounit,10) (ibuf(i),i=1,accu%tot_nbins)

 5  format(A20,1X,A1,1X,I10)
10  format(8(I9,1X))

end subroutine pmf_accu_write_ibuf

!===============================================================================
! Subroutine:  pmf_accu_read_rbuf
!===============================================================================

subroutine pmf_accu_read_rbuf(accu,iounit,keyline,rbuf)

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
    character(len=PMF_MAX_KEY)  :: dtype
    integer                     :: dsize
    integer                     :: i
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dtype,dsize

    if( dtype .ne. 'R' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require R but got "'//dtype//'" in pmf_accu_read_rbuf!')
    end if

    if( dsize .ne. accu%tot_nbins ) then
        write(serr,120) accu%tot_nbins, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    read(iounit,10,end=110,err=110) (rbuf(i),i=1,accu%tot_nbins)

    return

10  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_rbuf!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_rbuf!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_rbuf!')

end subroutine pmf_accu_read_rbuf

!===============================================================================
! Subroutine:  pmf_accu_write_rbuf
!===============================================================================

subroutine pmf_accu_write_rbuf(accu,iounit,key,rbuf)

    use pmf_dat
    use pmf_utils
    use pmf_unit
    use pmf_ver

    implicit none
    type(PMFAccuType)           :: accu
    integer                     :: iounit
    character(*)                :: key
    real(PMFDP)                 :: rbuf(:)
    ! -----------------------------------------------
    integer                     :: i
    !---------------------------------------------------------------------------

    write(iounit,5) adjustl(key), 'R', accu%tot_nbins
    write(iounit,10) (rbuf(i),i=1,accu%tot_nbins)

 5  format(A20,1X,A1,1X,I10)
10  format(4(E19.11,1X))

end subroutine pmf_accu_write_rbuf

!===============================================================================
! Subroutine:  pmf_accu_write_rbufCVs
!===============================================================================

subroutine pmf_accu_read_rbufCVs(accu,iounit,keyline,rbuf)

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
    character(len=PMF_MAX_KEY)  :: dtype
    integer                     :: dsize
    integer                     :: i,j
    character(len=PMF_MAX_PATH) :: serr
    !---------------------------------------------------------------------------

    read(keyline,*,end=100,err=100) skey,dtype,dsize

    if( dtype .ne. 'R' ) then
        call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Type mismatch: require R but got "'//dtype//'" in pmf_accu_read_rbufCVs!')
    end if

    if( dsize .ne. accu%tot_nbins*accu%tot_cvs ) then
        write(serr,120) accu%tot_nbins*accu%tot_cvs, dsize
        call pmf_utils_exit(PMF_OUT,1,serr)
    end if

    do i=1,accu%tot_cvs
        read(iounit,10,end=110,err=110) (rbuf(i,j),j=1,accu%tot_nbins)
    end do

    return

10  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Unable to parse keyline "'//keyline//'" in pmf_accu_read_rbufCVs!')
110 call pmf_utils_exit(PMF_OUT,1,'[PMFAccu] Premature end of data block for keyline "'//keyline//'" in pmf_accu_read_rbufCVs!')
120 format('[PMFAccu] Size mismatch: require ',I10,' but got ',I10,' in pmf_accu_read_rbufCVs!')

end subroutine pmf_accu_read_rbufCVs

!===============================================================================
! Subroutine:  pmf_accu_write_rbufCVs
!===============================================================================

subroutine pmf_accu_write_rbufCVs(accu,iounit,key,rbuf)

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
    integer                     :: i,j
    !---------------------------------------------------------------------------

    write(iounit,5) adjustl(key), 'R', accu%tot_nbins*accu%tot_cvs
    do i=1,accu%tot_cvs
        write(iounit,10) (rbuf(i,j),j=1,accu%tot_nbins)
    end do

 5  format(A20,1X,A1,1X,I10)
10  format(4(E19.11,1X))

end subroutine pmf_accu_write_rbufCVs

!===============================================================================

end module pmf_accu
