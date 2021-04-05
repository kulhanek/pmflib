!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
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

module abf_accumulator

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abf_accumulator_init
!===============================================================================

subroutine abf_accumulator_init()

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    allocate(accumulator%sizes(NumOfABFCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator!')
    endif

    tot_nbins = 1
    do i=1,NumOfABFCVs
        accumulator%sizes(i)%min_value  = ABFCVList(i)%min_value
        accumulator%sizes(i)%max_value  = ABFCVList(i)%max_value
        accumulator%sizes(i)%nbins      = ABFCVList(i)%nbins
        accumulator%sizes(i)%width      = abs(accumulator%sizes(i)%max_value - accumulator%sizes(i)%min_value)
        accumulator%sizes(i)%bin_width  = accumulator%sizes(i)%width / accumulator%sizes(i)%nbins
        tot_nbins = tot_nbins * accumulator%sizes(i)%nbins
    end do

    accumulator%tot_nbins = tot_nbins
    accumulator%tot_cvs = NumOfABFCVs

    ! ABF force arrays
    allocate(  &
            accumulator%weights(accumulator%tot_nbins), &
            accumulator%nsamples(accumulator%tot_nbins), &
            accumulator%icfsum(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%icfsum2(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%etotsum(accumulator%tot_nbins), &
            accumulator%etotsum2(accumulator%tot_nbins), &
            accumulator%icfetotsum(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%icfetotsum2(NumOfABFCVs,accumulator%tot_nbins), &
            stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (abfforce)!')
    endif

    if( feimode .eq. 3 ) then
        ! for block averages
        allocate(  &
            accumulator%block_nsamples(accumulator%tot_nbins), &
            accumulator%block_icfsum(NumOfABFCVs,accumulator%tot_nbins), &
            accumulator%block_etotsum(accumulator%tot_nbins), &
            accumulator%block_icfetotsum(NumOfABFCVs,accumulator%tot_nbins), &
            stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (block_icfsum)!')
        endif

    end if

    if( fserver_enabled ) then
        ! ABF incremental force arrays
        allocate(   accumulator%inc_nsamples(accumulator%tot_nbins), &
                    accumulator%inc_icfsum(NumOfABFCVs,accumulator%tot_nbins), &
                    accumulator%inc_icfsum2(NumOfABFCVs,accumulator%tot_nbins), &
                    accumulator%inc_etotsum(accumulator%tot_nbins), &
                    accumulator%inc_etotsum2(accumulator%tot_nbins), &
                    accumulator%inc_icfetotsum(NumOfABFCVs,accumulator%tot_nbins), &
                    accumulator%inc_icfetotsum2(NumOfABFCVs,accumulator%tot_nbins), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABF] Unable to allocate memory for abf accumulator (inc_icfsum)!')
        endif
    end if

    call abf_accumulator_clear()

    return

end subroutine abf_accumulator_init

!===============================================================================
! Subroutine:  abf_accumulator_clear
!===============================================================================

subroutine abf_accumulator_clear()

    use abf_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    accumulator%weights(:) = 1.0d0
    accumulator%nsamples(:) = 0
    accumulator%icfsum(:,:) = 0.0d0
    accumulator%icfsum2(:,:) = 0.0d0
    accumulator%etotsum(:) = 0.0d0
    accumulator%etotsum2(:) = 0.0d0
    accumulator%icfetotsum(:,:) = 0.0d0
    accumulator%icfetotsum2(:,:) = 0.0d0

    if( feimode .eq. 3 ) then
        accumulator%block_nsamples(:) = 0
        accumulator%block_icfsum(:,:) = 0.0d0
        accumulator%block_etotsum(:) = 0.0d0
        accumulator%block_icfetotsum(:,:) = 0.0d0
    end if

    if( fserver_enabled ) then
        accumulator%inc_nsamples(:) = 0
        accumulator%inc_icfsum(:,:) = 0.0d0
        accumulator%inc_icfsum2(:,:) = 0.0d0
        accumulator%inc_etotsum(:) = 0
        accumulator%inc_etotsum2(:) = 0
        accumulator%inc_icfetotsum(:,:) = 0
        accumulator%inc_icfetotsum2(:,:) = 0
    end if

end subroutine abf_accumulator_clear

!===============================================================================
! Function:  accumulator_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               accuvalue    ... value that is used to compute the bin index
! compute index for one accumulator coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function abf_accumulator_index(idxcoord,accuvalue)

    use abf_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: accuvalue
    ! --------------------------------------------------------------------------

    ! we need number from zero - therefore we use floor(x)
    abf_accumulator_index = floor((accuvalue - accumulator%sizes(idxcoord)%min_value) / &
                               accumulator%sizes(idxcoord)%bin_width)

    if( abf_accumulator_index .lt. 0 .or. abf_accumulator_index .ge.  accumulator%sizes(idxcoord)%nbins) then
        abf_accumulator_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function abf_accumulator_index

!===============================================================================
! Function:  accumulator_globalindex
! Description:  Compute globalindex for accumulator, based on accuvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function abf_accumulator_globalindex(lvalues)

    use abf_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    abf_accumulator_globalindex = 0

    do i=1,NumOfABFCVs
        idx_local = abf_accumulator_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            abf_accumulator_globalindex = -1
            return
        end if

        abf_accumulator_globalindex = abf_accumulator_globalindex*accumulator%sizes(i)%nbins + idx_local
    end do

    abf_accumulator_globalindex = abf_accumulator_globalindex + 1

    return

end function abf_accumulator_globalindex

!===============================================================================
! Subroutine:  abf_accumulator_read
!===============================================================================

subroutine abf_accumulator_read(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer             :: iounit
    ! -----------------------------------------------
    integer             :: ibitem
    character(len=4)    :: sabf
    character(len=2)    :: sver
    character(len=80)   :: buffer
    ! --------------------------------------------------------------------------

    read(iounit,'(A80)') buffer

    ! read header --------------------------
    read(buffer,*,end=200,err=200) sabf, sver, ibitem

    if( trim(adjustl(sabf)) .ne. 'ABF' ) then
        write(PMF_OUT,*) '[ABF] header v3/v4 key = [',sabf,']'
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Missing ABF key in ABF accumulator v3/v4 header!')
    end if

    if( ibitem .ne. NumOfABFCVs ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] ABF accumulator v3/v4 contains different number of CVs!')
    end if

    if( trim(adjustl(sver)) .eq. 'V3' ) then
        call abf_accumulator_read_v3(iounit)
        return
    end if

    if( trim(adjustl(sver)) .eq. 'V4' ) then
        call abf_accumulator_read_v4(iounit)
        return
    end if

    if( trim(adjustl(sver)) .eq. 'V5' ) then
        call abf_accumulator_read_v5(iounit)
        return
    end if

    write(PMF_OUT,*) '[ABF] header ver = [',sver,']'
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported version key in ABF accumulator header!')

    return

    ! unsupported header - process error
200 continue
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Illegal header in ABF accumulator (Only V3 header is supported)!')

end subroutine abf_accumulator_read

!===============================================================================
! Subroutine:  abf_accumulator_read_v4
!===============================================================================

subroutine abf_accumulator_read_v3(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    integer                         :: i,j,it,nbins
    character(len=PMF_MAX_TYPE)     :: stype
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: min_value,max_value
    ! --------------------------------------------------------------------------

    ! read header --------------------------
    do i=1, NumOfABFCVs
        ! read CV definition
        read(iounit,20,end=102,err=102) it, stype, min_value, max_value, nbins
        ! check CV definition
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Inccorect item in ABF accumulator v3!')
        end if
        if( trim(adjustl(stype)) .ne. trim(ABFCVList(i)%cv%ctype) ) then
            write(PMF_OUT,*) '[ABF] CV type = [',trim(adjustl(stype)),'] should be [',trim(ABFCVList(i)%cv%ctype),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV type was redefined in ABF accumulator!')
        end if
        if( abs(min_value-accumulator%sizes(i)%min_value) .gt. abs(accumulator%sizes(i)%min_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Minimal value of CV was redefined in ABF accumulator!')
        end if
        if( abs(max_value-accumulator%sizes(i)%max_value) .gt. abs(accumulator%sizes(i)%max_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Maximal value of CV was redefined in ABF accumulator!')
        end if
        if( nbins .ne. accumulator%sizes(i)%nbins ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Number of CV bins was redefined in ABF accumulator!')
        end if

        ! read names
        read(iounit,25,end=102,err=102) it, sname
        ! check names
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Inccorect item in ABF accumulator!')
        end if
        if( trim(adjustl(sname)) .ne. trim(ABFCVList(i)%cv%name) ) then
            write(PMF_OUT,*) '[ABF] CV name = [',trim(adjustl(sname)),'] should be [',trim(ABFCVList(i)%cv%name),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV name was redefined in ABF accumulator!')
        end if
    end do

    ! read accumulator data - ABF forces
    if( accumulator%tot_nbins .gt. 0 ) then
        read(iounit,30,end=100,err=100) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)

        do i=1,NumOfABFCVs
            read(iounit,40,end=100,err=100) (accumulator%icfsum(i,j),j=1,accumulator%tot_nbins)
        end do
        do i=1,NumOfABFCVs
            read(iounit,40,end=100,err=100) (accumulator%icfsum2(i,j),j=1,accumulator%tot_nbins)
        end do
    end if

    return

20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v3 - data section!')
102 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v3 - CV section!')

end subroutine abf_accumulator_read_v3

!===============================================================================
! Subroutine:  abf_accumulator_read_v4
!===============================================================================

subroutine abf_accumulator_read_v4(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    integer                         :: i,j,it,nbins
    character(len=PMF_MAX_TYPE)     :: stype,sunit
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: min_value,max_value,fconv
    ! --------------------------------------------------------------------------

    ! read header --------------------------
    do i=1, NumOfABFCVs
        ! read CV definition
        read(iounit,20,end=201,err=201) it, stype, min_value, max_value, nbins
        ! check CV definition
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF accumulator v4!')
        end if
        if( trim(adjustl(stype)) .ne. trim(ABFCVList(i)%cv%ctype) ) then
            write(PMF_OUT,*) '[ABF] CV type = [',trim(adjustl(stype)),'] should be [',trim(ABFCVList(i)%cv%ctype),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV type was redefined in ABF accumulator!')
        end if
        if( abs(min_value-accumulator%sizes(i)%min_value) .gt. abs(accumulator%sizes(i)%min_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Minimal value of CV was redefined in ABF accumulator!')
        end if
        if( abs(max_value-accumulator%sizes(i)%max_value) .gt. abs(accumulator%sizes(i)%max_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Maximal value of CV was redefined in ABF accumulator!')
        end if
        if( nbins .ne. accumulator%sizes(i)%nbins ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Number of CV bins was redefined in ABF accumulator!')
        end if

        ! read names
        read(iounit,25,end=202,err=202) it, sname
        ! check names
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF accumulator!')
        end if
        if( trim(adjustl(sname)) .ne. trim(ABFCVList(i)%cv%name) ) then
            write(PMF_OUT,*) '[ABF] CV name = [',trim(adjustl(sname)),'] should be [',trim(ABFCVList(i)%cv%name),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV name was redefined in ABF accumulator!')
        end if

        ! read names
        read(iounit,26,end=203,err=203) it, fconv, sunit
        ! check names
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF accumulator!')
        end if
        ! ignore values fconv and sunit
    end do

    ! read accumulator data - ABF forces
    if( accumulator%tot_nbins .gt. 0 ) then
        read(iounit,30,end=100,err=100) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)

        do i=1,NumOfABFCVs
            read(iounit,40,end=101,err=101) (accumulator%icfsum(i,j),j=1,accumulator%tot_nbins)
        end do
        do i=1,NumOfABFCVs
            read(iounit,40,end=102,err=102) (accumulator%icfsum2(i,j),j=1,accumulator%tot_nbins)
        end do

        read(iounit,40,end=103,err=103) (accumulator%etotsum(i),i=1,accumulator%tot_nbins)
        read(iounit,40,end=104,err=104) (accumulator%etotsum2(i),i=1,accumulator%tot_nbins)
    end if

    return

20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - data section - bins!')
101 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - data section - sabf!')
102 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - data section - sabf2!')
103 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - data section - spot!')
104 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - data section - spot2!')

201 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - CV section - def!')
202 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - CV section - name!')
203 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v4 - CV section - unit!')

end subroutine abf_accumulator_read_v4

!===============================================================================
! Subroutine:  abf_accumulator_read_v5
!===============================================================================

subroutine abf_accumulator_read_v5(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    integer                         :: i,j,it,nbins
    character(len=PMF_MAX_TYPE)     :: stype,sunit
    character(len=PMF_MAX_KEY)      :: key
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: min_value,max_value,fconv
    ! --------------------------------------------------------------------------

    do while(.true.)

        ! read key
        read(iounit,5,end=500,err=300) key

        select case( trim(key) )
            case('CVS')
                ! read header --------------------------
                do i=1, NumOfABFCVs
                    ! read CV definition
                    read(iounit,20,end=201,err=201) it, stype, min_value, max_value, nbins
                    ! check CV definition
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF accumulator v4!')
                    end if
                    if( trim(adjustl(stype)) .ne. trim(ABFCVList(i)%cv%ctype) ) then
                        write(PMF_OUT,*) '[ABF] CV type = [',trim(adjustl(stype)),'] should be [',trim(ABFCVList(i)%cv%ctype),']'
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] CV type was redefined in ABF accumulator!')
                    end if
                    if( abs(min_value-accumulator%sizes(i)%min_value) .gt. abs(accumulator%sizes(i)%min_value/100000.0d0) ) then
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Minimal value of CV was redefined in ABF accumulator!')
                    end if
                    if( abs(max_value-accumulator%sizes(i)%max_value) .gt. abs(accumulator%sizes(i)%max_value/100000.0d0) ) then
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Maximal value of CV was redefined in ABF accumulator!')
                    end if
                    if( nbins .ne. accumulator%sizes(i)%nbins ) then
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Number of CV bins was redefined in ABF accumulator!')
                    end if

                    ! read names
                    read(iounit,25,end=202,err=202) it, sname
                    ! check names
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF accumulator!')
                    end if
                    if( trim(adjustl(sname)) .ne. trim(ABFCVList(i)%cv%name) ) then
                        write(PMF_OUT,*) '[ABF] CV name = [',trim(adjustl(sname)),'] should be [',trim(ABFCVList(i)%cv%name),']'
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] CV name was redefined in ABF accumulator!')
                    end if

                    ! read names
                    read(iounit,26,end=203,err=203) it, fconv, sunit
                    ! check names
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF accumulator!')
                    end if
                    ! ignore values fconv and sunit
                end do
            case('TEMP')
                ! read but do not use
                read(iounit,6,end=302,err=302) fconv
            case('ENERGY')
                ! read but do not use
                read(iounit,27,end=301,err=301) fconv, sunit
            case('NSAMPLES')
                read(iounit,30,end=100,err=100) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)
            case('ABF')
                do i=1,NumOfABFCVs
                    read(iounit,40,end=101,err=101) (accumulator%icfsum(i,j),j=1,accumulator%tot_nbins)
                end do
                do i=1,NumOfABFCVs
                    read(iounit,40,end=102,err=102) (accumulator%icfsum2(i,j),j=1,accumulator%tot_nbins)
                end do
            case('ABF*EPOT')
                do i=1,NumOfABFCVs
                    read(iounit,40,end=101,err=105) (accumulator%icfetotsum(i,j),j=1,accumulator%tot_nbins)
                end do
                do i=1,NumOfABFCVs
                    read(iounit,40,end=102,err=106) (accumulator%icfetotsum2(i,j),j=1,accumulator%tot_nbins)
                end do
            case('EPOT')
                read(iounit,40,end=103,err=103) (accumulator%etotsum(i),i=1,accumulator%tot_nbins)
                read(iounit,40,end=104,err=104) (accumulator%etotsum2(i),i=1,accumulator%tot_nbins)
            case default
                call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - unrecognized keyword: '//trim(key))
        end select
    end do

500 return

 5  format(A20)
 6  format(F10.4)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
27  format(E18.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

300 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - keyword!')

100 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - bins!')
101 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - sabf!')
102 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - sabf2!')
103 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - spot!')
104 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - spot2!')
105 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - sabf*etot!')
106 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - data section - sabf*etot2!')

201 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - CV section - def!')
202 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - CV section - name!')
203 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - CV section - unit!')

301 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - energy unit!')
302 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABF accumulator v5 - temperature!')

end subroutine abf_accumulator_read_v5

!===============================================================================
! Subroutine:  abf_accumulator_read_mask
!===============================================================================

subroutine abf_accumulator_read_mask(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    integer                         :: ibitem,i,it,nbins
    character(len=5)                :: sabf,sver
    character(len=PMF_MAX_TYPE)     :: stype
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: min_value,max_value
    ! --------------------------------------------------------------------------

    ! read header - V3 --------------------------
    read(iounit,10,end=100,err=100) sabf, sver, ibitem

    if( sabf .ne. 'MABF' ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Missing MABF key in ABF mask header!')
    end if

    if( sver .ne. 'V3' ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] Unsupported version in ABF mask header!')
    end if

    if( ibitem .ne. NumOfABFCVs ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABF] ABF mask contains different number of CVs!')
    end if

    ! read CV header --------------------------
    do i=1, NumOfABFCVs
        read(iounit,20,end=100,err=100) it, stype, min_value, max_value, nbins
        ! check fingerprint
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF mask!')
        end if
        if( trim(adjustl(stype)) .ne. trim(ABFCVList(i)%cv%ctype) ) then
            write(PMF_OUT,*) '[ABF] CV type = [',trim(adjustl(stype)),'] should be [',trim(ABFCVList(i)%cv%ctype),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV type was redefined in ABF mask!')
        end if
        if( abs(min_value-accumulator%sizes(i)%min_value) .gt. abs(accumulator%sizes(i)%min_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Minimal value of CV was redefined in ABF mask!')
        end if
        if( abs(max_value-accumulator%sizes(i)%max_value) .gt. abs(accumulator%sizes(i)%max_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Maximal value of CV was redefined in ABF mask!')
        end if
        if( nbins .ne. accumulator%sizes(i)%nbins ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Number of CV bins was redefined in ABF mask!')
        end if
        ! read names
        read(iounit,25,end=100,err=100) it, sname
        ! check fingerprint
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Incorrect item in ABF mask!')
        end if
        if( trim(adjustl(sname)) .ne. trim(ABFCVList(i)%cv%name) ) then
            write(PMF_OUT,*) '[ABF] CV name = [',trim(adjustl(sname)),'] should be [',trim(ABFCVList(i)%cv%name),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABF] CV name was redefined in ABF mask!')
        end if
    end do

    ! read mask weights
    if( accumulator%tot_nbins .gt. 0 ) then
        read(iounit,40,end=100,err=100) (accumulator%weights(i),i=1,accumulator%tot_nbins)
    end if

    return

    ! unsupported header - process error
100 continue
    call pmf_utils_exit(PMF_OUT,1,'[ABF] Illegal header in ABF mask!')

10  format(A4,1X,A2,1X,I2)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
40  format(4(E19.11,1X))

end subroutine abf_accumulator_read_mask

!===============================================================================
! Subroutine:  abf_accumulator_write
!===============================================================================

subroutine abf_accumulator_write(iounit)

    use abf_dat
    use pmf_dat
    use pmf_utils
    use pmf_unit

    implicit none
    integer                     :: iounit
    ! -----------------------------------------------
    integer                     :: i,j
    character(len=PMF_MAX_KEY)  :: key
    !---------------------------------------------------------------------------

    ! write header --------------------------
    write(iounit,10) 'ABF ', 'V5', NumOfABFCVs

    key = 'CVS'
    write(iounit,5) adjustl(key)
    do i=1, NumOfABFCVs
        write(iounit,20) i,trim(ABFCVList(i)%cv%ctype), &
                          accumulator%sizes(i)%min_value,accumulator%sizes(i)%max_value, &
                          accumulator%sizes(i)%nbins
        write(iounit,25) i,trim(ABFCVList(i)%cv%name)
        write(iounit,26) i,pmf_unit_get_rvalue(ABFCVList(i)%cv%unit,1.0d0),trim(pmf_unit_label(ABFCVList(i)%cv%unit))
    end do

    key = 'TEMPERATURE'
    write(iounit,5) adjustl(key)
    write(iounit,6) ftemp

    key = 'ENERGY'
    write(iounit,5) adjustl(key)
    write(iounit,27) pmf_unit_get_rvalue(EnergyUnit,1.0d0),trim(pmf_unit_label(EnergyUnit))

    key = 'NSAMPLES'
    write(iounit,5) adjustl(key)
    write(iounit,30) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)

    key = 'ICF'
    write(iounit,5) adjustl(key)
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%icfsum(i,j),j=1,accumulator%tot_nbins)
    end do
    do i=1,NumOfABFCVs
        write(iounit,40) (accumulator%icfsum2(i,j),j=1,accumulator%tot_nbins)
    end do

    if( faccuepot .or. faccuekin ) then
        if( faccuekin ) then
            key = 'ETOT=EPOT+EKIN'
        else
            key = 'ETOT=EPOT'
        end if
        write(iounit,5) adjustl(key)
        key = 'ETOT'
        write(iounit,5) adjustl(key)
        write(iounit,40) (accumulator%etotsum(i),i=1,accumulator%tot_nbins)
        write(iounit,40) (accumulator%etotsum2(i),i=1,accumulator%tot_nbins)

        key = 'ICF*ETOT'
        write(iounit,5) adjustl(key)
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%icfetotsum(i,j),j=1,accumulator%tot_nbins)
        end do
        do i=1,NumOfABFCVs
            write(iounit,40) (accumulator%icfetotsum2(i,j),j=1,accumulator%tot_nbins)
        end do
    end if

    return

 5  format(A20)
 6  format(F10.4)
10  format(A4,1X,A2,1X,I2)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
27  format(3X,E18.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

end subroutine abf_accumulator_write

!===============================================================================
! Subroutine:  abf_accumulator_add_data12
!===============================================================================

subroutine abf_accumulator_add_data12(cvs,lpotene,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: lpotene
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0,i
    real(PMFDP)    :: a
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = abf_accumulator_globalindex(cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    accumulator%nsamples(gi0) = accumulator%nsamples(gi0) + 1
    do i=1,NumOfABFCVs
        a = gfx(i)
        accumulator%icfsum(i,gi0)  = accumulator%icfsum(i,gi0) - a
        accumulator%icfsum2(i,gi0) = accumulator%icfsum2(i,gi0) + a**2

        a = a * lpotene
        accumulator%icfetotsum(i,gi0)  = accumulator%icfetotsum(i,gi0) - a
        accumulator%icfetotsum2(i,gi0) = accumulator%icfetotsum2(i,gi0) + a**2
    end do

    ! potential energy
    accumulator%etotsum(gi0)  = accumulator%etotsum(gi0)  + lpotene
    accumulator%etotsum2(gi0) = accumulator%etotsum2(gi0) + lpotene**2

    if( fserver_enabled ) then
        accumulator%inc_nsamples(gi0) = accumulator%inc_nsamples(gi0) + 1
        do i=1,NumOfABFCVs
            a = gfx(i)
            accumulator%inc_icfsum(i,gi0)  = accumulator%inc_icfsum(i,gi0) - a
            accumulator%inc_icfsum2(i,gi0) = accumulator%inc_icfsum2(i,gi0) + a**2

            a = a * lpotene
            accumulator%inc_icfetotsum(i,gi0)  = accumulator%inc_icfetotsum(i,gi0) - a
            accumulator%inc_icfetotsum2(i,gi0) = accumulator%inc_icfetotsum2(i,gi0) + a**2
        end do
        ! potential energy
        accumulator%inc_etotsum(gi0)  = accumulator%inc_etotsum(gi0)  + lpotene
        accumulator%inc_etotsum2(gi0) = accumulator%inc_etotsum2(gi0) + lpotene**2
    end if

end subroutine abf_accumulator_add_data12

!===============================================================================
! Subroutine:  abf_accumulator_add_data3
!===============================================================================

subroutine abf_accumulator_add_data3(cvs,lpotene,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: lpotene
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0,i
    real(PMFDP)    :: a
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = abf_accumulator_globalindex(cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    end if

    insidesamples = insidesamples + 1

    ! update block
    accumulator%block_nsamples(gi0)      = accumulator%block_nsamples(gi0)      + 1
    accumulator%block_icfsum(:,gi0)      = accumulator%block_icfsum(:,gi0)      + gfx(:)
    accumulator%block_etotsum(gi0)       = accumulator%block_etotsum(gi0)       + lpotene
    accumulator%block_icfetotsum(:,gi0)  = accumulator%block_icfetotsum(:,gi0)  + lpotene * gfx(:)

    ! is block filled?
    if( accumulator%block_nsamples(gi0) .lt. fblock_size ) then
        ! NO -> exit
        return
    end if

    ! populate abf accumulator

    ! increase number of samples
    accumulator%nsamples(gi0) = accumulator%nsamples(gi0) + 1
    do i=1,NumOfABFCVs
        a = accumulator%block_icfsum(i,gi0) / real(accumulator%block_nsamples(gi0),PMFDP)
        accumulator%icfsum(i,gi0)  = accumulator%icfsum(i,gi0)  - a
        accumulator%icfsum2(i,gi0) = accumulator%icfsum2(i,gi0) + a**2

        a = accumulator%block_icfetotsum(i,gi0) / real(accumulator%block_nsamples(gi0),PMFDP)
        accumulator%icfetotsum(i,gi0)  = accumulator%icfetotsum(i,gi0)  - a
        accumulator%icfetotsum2(i,gi0) = accumulator%icfetotsum2(i,gi0) + a**2
    end do
    a = accumulator%block_etotsum(gi0) / real(accumulator%block_nsamples(gi0),PMFDP)
    accumulator%etotsum(gi0)  = accumulator%etotsum(gi0) + a
    accumulator%etotsum2(gi0) = accumulator%etotsum2(gi0) + a**2

    if( fserver_enabled ) then
        accumulator%inc_nsamples(gi0) = accumulator%inc_nsamples(gi0) + 1
        do i=1,NumOfABFCVs
            a = accumulator%block_icfsum(i,gi0) / real(accumulator%block_nsamples(gi0),PMFDP)
            accumulator%inc_icfsum(i,gi0)  = accumulator%inc_icfsum(i,gi0) - a
            accumulator%inc_icfsum2(i,gi0) = accumulator%inc_icfsum2(i,gi0) + a**2

            a = accumulator%block_icfetotsum(i,gi0) / real(accumulator%block_nsamples(gi0),PMFDP)
            accumulator%inc_icfetotsum(i,gi0)  = accumulator%inc_icfetotsum(i,gi0)  - a
            accumulator%inc_icfetotsum2(i,gi0) = accumulator%inc_icfetotsum2(i,gi0) + a**2
        end do
        a = accumulator%block_etotsum(gi0) / real(accumulator%block_nsamples(gi0),PMFDP)
        accumulator%inc_etotsum(gi0)  = accumulator%inc_etotsum(gi0) + a
        accumulator%inc_etotsum2(gi0) = accumulator%inc_etotsum2(gi0) + a**2
    end if

    ! reset the block
    accumulator%block_nsamples(gi0)         = 0
    accumulator%block_icfsum(:,gi0)         = 0.0d0
    accumulator%block_etotsum(gi0)          = 0.0d0
    accumulator%block_icfetotsum(:,gi0)     = 0.0d0

end subroutine abf_accumulator_add_data3

!===============================================================================
! Subroutine:  abf_accumulator_get_data1
!===============================================================================

subroutine abf_accumulator_get_data1(values,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0,n
    real(PMFDP)    :: sc_ramp, w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = abf_accumulator_globalindex(values)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get number of samples
    n = accumulator%nsamples(gi0)
    if( n .gt. 0 ) then
        sc_ramp = 1.0d0
        if( n .le. fhramp ) then
            sc_ramp = real(n)/real(fhramp)
        end if
        w      = accumulator%weights(gi0)
        gfx(:) = w * sc_ramp * accumulator%icfsum(:,gi0) / real(n)
    end if

end subroutine abf_accumulator_get_data1

!===============================================================================
! Subroutine:  abf_accumulator_get_data2
!===============================================================================

subroutine abf_accumulator_get_data2(values,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0,n
    real(PMFDP)    :: sc_ramp, w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = abf_accumulator_globalindex(values)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get number of samples
    n = accumulator%nsamples(gi0)
    if( n .gt. 0 ) then
        sc_ramp = 1.0d0
        if( n .le. fhramp_max ) then
            sc_ramp = 0.0d0
            if( n .gt. fhramp_min ) then
                sc_ramp = real(n-fhramp_min)/real(fhramp_max-fhramp_min)
            end if
        end if
        w      = accumulator%weights(gi0)
        gfx(:) = w * sc_ramp * accumulator%icfsum(:,gi0) / real(n)
    end if

end subroutine abf_accumulator_get_data2

!===============================================================================
! Subroutine:  abf_accumulator_get_data3
!===============================================================================

subroutine abf_accumulator_get_data3(values,gfx)

    use abf_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: gfx(:)
    ! -----------------------------------------------
    integer        :: gi0, n
    real(PMFDP)    :: w
    ! --------------------------------------------------------------------------

    gfx(:) = 0.0d0

    ! get global index to accumulator for average values within the set
    gi0 = abf_accumulator_globalindex(values)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get accumulated mean forces
    n = accumulator%nsamples(gi0)
    if( n .gt. 0 ) then
        ! get mean forces
        w       = accumulator%weights(gi0)
        gfx(:)  = w * accumulator%icfsum(:,gi0) / real(n,PMFDP)
    end if

end subroutine abf_accumulator_get_data3

!===============================================================================

end module abf_accumulator

