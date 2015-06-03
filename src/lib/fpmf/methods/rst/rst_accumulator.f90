!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module rst_accumulator

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  rst_accumulator_init
!===============================================================================

subroutine rst_accumulator_init()

    use rst_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! only if we would like to update restart file and of if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    ! init dimensions ------------------------------
    allocate(Accumulator%sizes(NumOfRSTItems), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[RST] Unable to allocate memory for rst accumulator!')
    endif

    tot_nbins = 1
    do i=1,NumOfRSTItems
        Accumulator%sizes(i)%min_value  = RSTCVList(i)%min_value
        Accumulator%sizes(i)%max_value  = RSTCVList(i)%max_value
        Accumulator%sizes(i)%nbins      = RSTCVList(i)%nbins
        Accumulator%sizes(i)%width      = abs(Accumulator%sizes(i)%max_value - Accumulator%sizes(i)%min_value)
        Accumulator%sizes(i)%bin_width  = Accumulator%sizes(i)%width / Accumulator%sizes(i)%nbins
        tot_nbins = tot_nbins * Accumulator%sizes(i)%nbins
    end do

    Accumulator%tot_nbins = tot_nbins

    ! RST force arrays
    allocate( Accumulator%nsamples(Accumulator%tot_nbins), &
              stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[RST] Unable to allocate memory for rst accumulator (rstforce)!')
    endif

    call rst_accumulator_clear()

    return

end subroutine rst_accumulator_init

!===============================================================================
! Subroutine:  rst_accumulator_clear
!===============================================================================

subroutine rst_accumulator_clear()

    use rst_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    Accumulator%nsamples(:) = 0

end subroutine rst_accumulator_clear

!===============================================================================
! Function:  accumulator_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               accuvalue    ... value that is used to compute the bin index
! compute index for one accumulator coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function rst_accumulator_index(idxcoord,accuvalue)

    use rst_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: accuvalue
    ! --------------------------------------------------------------------------

    ! we need number from zero - therefore we use floor(x)
    rst_accumulator_index = floor((accuvalue - Accumulator%sizes(idxcoord)%min_value) / &
                               Accumulator%sizes(idxcoord)%bin_width)

    if( rst_accumulator_index .lt. 0 .or. rst_accumulator_index .ge.  Accumulator%sizes(idxcoord)%nbins) then
        rst_accumulator_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function rst_accumulator_index

!===============================================================================
! Function:  accumulator_globalindex
! Description:  Compute globalindex for accumulator, based on accuvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function rst_accumulator_globalindex(lvalues)

    use rst_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    rst_accumulator_globalindex = 0

    do i=1,NumOfRSTItems
        idx_local = rst_accumulator_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            rst_accumulator_globalindex = -1
            return
        end if

        rst_accumulator_globalindex = rst_accumulator_globalindex*Accumulator%sizes(i)%nbins + idx_local
    end do

    rst_accumulator_globalindex = rst_accumulator_globalindex + 1

    return

end function rst_accumulator_globalindex

!===============================================================================
! Subroutine:  rst_accumulator_read
!===============================================================================

subroutine rst_accumulator_read(iounit)

    use rst_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer             :: iounit
    ! -----------------------------------------------
    integer             :: ibitem
    character(len=4)    :: srst
    character(len=2)    :: sver
    character(len=80)   :: buffer
    ! --------------------------------------------------------------------------

    read(iounit,'(A80)') buffer

    ! read header - V1 --------------------------
    read(buffer,20,end=200,err=200) srst, sver, ibitem

    if( trim(adjustl(srst)) .ne. 'RST' ) then
        write(PMF_OUT,*) '[RST] header v1 key = [',srst,']'
        call pmf_utils_exit(PMF_OUT,1,'[RST] Missing RST key in RST accumulator v1 header!')
    end if

    if( trim(adjustl(sver)) .ne. 'V1' ) then
        write(PMF_OUT,*) '[RST] header v1 ver = [',sver,']'
        call pmf_utils_exit(PMF_OUT,1,'[RST] Unsupported version key in RST accumulator v1 header!')
    end if

    if( ibitem .ne. NumOfRSTItems ) then
        call pmf_utils_exit(PMF_OUT,1,'[RST] RST accumulator v1 contains different number of CVs!')
    end if

    call rst_accumulator_read_v1(iounit)
    return

    ! unsupported header - process error
200 continue
    call pmf_utils_exit(PMF_OUT,1,'[RST] Illegal header in RST accumulator (Only V1 header is supported)!')

20  format(A3,1X,A2,1X,I2)

end subroutine rst_accumulator_read

!===============================================================================
! Subroutine:  rst_accumulator_read_v1
!===============================================================================

subroutine rst_accumulator_read_v1(iounit)

    use rst_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                         :: iounit
    ! -----------------------------------------------
    integer                         :: i,it,nbins
    character(len=PMF_MAX_TYPE)     :: stype
    character(len=PMF_MAX_CV_NAME)  :: sname
    real(PMFDP)                     :: min_value,max_value
    ! --------------------------------------------------------------------------

    ! read header --------------------------
    do i=1, NumOfRSTItems
        ! read CV definition
        read(iounit,20,end=102,err=102) it, stype, min_value, max_value, nbins
        ! check CV definition
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] Inccorect item in RST accumulator!')
        end if
        if( trim(adjustl(stype)) .ne. trim(RSTCVList(i)%cv%ctype) ) then
            write(PMF_OUT,*) '[RST] CV type = [',trim(adjustl(stype)),'] should be [',trim(RSTCVList(i)%cv%ctype),']'
            call pmf_utils_exit(PMF_OUT,1,'[RST] CV type was redefined in RST accumulator!')
        end if
        if( abs(min_value-Accumulator%sizes(i)%min_value) .gt. abs(Accumulator%sizes(i)%min_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] Minimal value of CV was redefined in RST accumulator!')
        end if
        if( abs(max_value-Accumulator%sizes(i)%max_value) .gt. abs(Accumulator%sizes(i)%max_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] Maximal value of CV was redefined in RST accumulator!')
        end if
        if( nbins .ne. Accumulator%sizes(i)%nbins ) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] Number of CV bins was redefined in RST accumulator!')
        end if

        ! read names
        read(iounit,25,end=102,err=102) it, sname
        ! check names
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[RST] Inccorect item in RST accumulator!')
        end if
        if( trim(adjustl(sname)) .ne. trim(RSTCVList(i)%cv%name) ) then
            write(PMF_OUT,*) '[RST] CV name = [',trim(adjustl(sname)),'] should be [',trim(RSTCVList(i)%cv%name),']'
            call pmf_utils_exit(PMF_OUT,1,'[RST] CV name was redefined in RST accumulator!')
        end if
    end do

    ! read accumulator data - RST forces
    if( Accumulator%tot_nbins .gt. 0 ) then
        read(iounit,30,end=100,err=100) (Accumulator%nsamples(i),i=1,Accumulator%tot_nbins)
    end if

    return

20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
30  format(8(I9,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to read from RST accumulator v1 - data section!')
102 call pmf_utils_exit(PMF_OUT,1,'[RST] Unable to read from RST accumulator v1 - cv section!')

end subroutine rst_accumulator_read_v1

!===============================================================================
! Subroutine:  rst_accumulator_write
!===============================================================================

subroutine rst_accumulator_write(iounit)

    use rst_dat
    use pmf_dat
    use pmf_utils
    use pmf_unit

    implicit none
    integer                :: iounit
    ! -----------------------------------------------
    integer                :: i
    !---------------------------------------------------------------------------

    ! write header --------------------------
    write(iounit,10) 'RST', 'V1', NumOfRSTItems
    do i=1, NumOfRSTItems
        write(iounit,20) i,trim(RSTCVList(i)%cv%ctype), &
                          Accumulator%sizes(i)%min_value,Accumulator%sizes(i)%max_value, &
                          Accumulator%sizes(i)%nbins
        write(iounit,25) i,trim(RSTCVList(i)%cv%name)
    end do

    ! write accumulator data - RST samples
    write(iounit,30) (Accumulator%nsamples(i),i=1,Accumulator%tot_nbins)

    return

10  format(A3,1X,A2,1X,I2)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
30  format(8(I9,1X))

end subroutine rst_accumulator_write

!===============================================================================
! Subroutine:  rst_accumulator_add_sample
!===============================================================================

subroutine rst_accumulator_add_sample(values)

    use rst_dat
    use pmf_dat

    implicit none
    real(PMFDP)    :: values(:)
    ! -----------------------------------------------
    integer        :: gi0
    ! --------------------------------------------------------------------------

    ! only if we would like to update restart file and of if restart is explicitly required
    if( (fhistupdate .eq. 0) .and. (frestart .eqv. .false.) ) return

    ! get global index to accumulator for average values within the set
    gi0 = rst_accumulator_globalindex(values)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
        return ! out of valid area
    else
        insidesamples = insidesamples + 1
    end if

    ! increase number of samples
    Accumulator%nsamples(gi0) = Accumulator%nsamples(gi0) + 1

end subroutine rst_accumulator_add_sample

!===============================================================================

end module rst_accumulator

