!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abp_accumulator

use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  abp_accumulator_init
!===============================================================================

subroutine abp_accumulator_init()

    use abp_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins,gi
    integer              :: alloc_failed
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    allocate(accumulator%sizes(NumOfABPCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABP] Unable to allocate memory for abp accumulator!')
    endif

    tot_nbins = 1
    do i=1,NumOfABPCVs
        accumulator%sizes(i)%min_value  = ABPCVList(i)%min_value
        accumulator%sizes(i)%max_value  = ABPCVList(i)%max_value
        accumulator%sizes(i)%nbins      = ABPCVList(i)%nbins
        accumulator%sizes(i)%width      = abs(accumulator%sizes(i)%max_value - accumulator%sizes(i)%min_value)
        accumulator%sizes(i)%bin_width  = accumulator%sizes(i)%width / accumulator%sizes(i)%nbins
        tot_nbins = tot_nbins * accumulator%sizes(i)%nbins
    end do

    accumulator%tot_nbins = tot_nbins

    ! ABP force arrays
    allocate(  accumulator%nsamples(accumulator%tot_nbins), &
            accumulator%dpop(NumOfABPCVs,accumulator%tot_nbins), &
            accumulator%pop(accumulator%tot_nbins), &
            accumulator%binpos(NumOfABPCVs,accumulator%tot_nbins), &
            stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[ABP] Unable to allocate memory for abp accumulator!')
    endif

    if( fserver_enabled ) then
        ! ABP incremental arrays
        allocate(  accumulator%nisamples(accumulator%tot_nbins), &
                    accumulator%idpop(NumOfABPCVs,accumulator%tot_nbins), &
                    accumulator%ipop(accumulator%tot_nbins), &
                    stat = alloc_failed)

        if( alloc_failed .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT, 1,'[ABP] Unable to allocate memory for abf accumulator (iabfforce)!')
        endif
    end if

    ! init binpos
    do gi=1,accumulator%tot_nbins
        ! get CV values on a grid point
        call abp_accumulator_get_values(gi,accumulator%binpos(:,gi))
    end do

    call abp_accumulator_clear()

end subroutine abp_accumulator_init

!===============================================================================
! Subroutine:  abp_accumulator_clear
!===============================================================================

subroutine abp_accumulator_clear()

    use abp_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    accumulator%nsamples(:) = 0
    accumulator%dpop(:,:) = 0.0d0
    accumulator%pop(:) = 1.0d0
    accumulator%m = 1.0d0

    if( fserver_enabled ) then
        accumulator%nisamples(:) = 0
        accumulator%idpop(:,:) = 0.0d0
        accumulator%ipop(:) = 0.0d0
    end if

end subroutine abp_accumulator_clear

!===============================================================================
! Function:  abp_accumulator_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               accuvalue    ... value that is used to compute the bin index
! compute index for one accumulator coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function abp_accumulator_index(idxcoord,accuvalue)

    use abp_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: accuvalue
    ! --------------------------------------------------------------------------

    ! we need number from zero - therefore we use floor(x)
    abp_accumulator_index = floor((accuvalue - accumulator%sizes(idxcoord)%min_value) / &
                               accumulator%sizes(idxcoord)%bin_width)

    if( abp_accumulator_index .lt. 0 .or. abp_accumulator_index .ge.  accumulator%sizes(idxcoord)%nbins) then
        abp_accumulator_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

end function abp_accumulator_index

!===============================================================================
! Function:  abp_accumulator_globalindex
! Description:  Compute globalindex for accumulator, based on accuvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function abp_accumulator_globalindex(lvalues)

    use abp_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    abp_accumulator_globalindex = 0

    do i=1,NumOfABPCVs
        idx_local = abp_accumulator_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            abp_accumulator_globalindex = -1
            return
        end if

        abp_accumulator_globalindex = abp_accumulator_globalindex*accumulator%sizes(i)%nbins + idx_local
    end do

    abp_accumulator_globalindex = abp_accumulator_globalindex + 1

    return

end function abp_accumulator_globalindex

!===============================================================================
! Function:  abp_accumulator_bin_pos
! Arguments:
!               idxcoord ... number of ksi coordinate
!               bindex    ... bin index (0,1,2,...)
!===============================================================================

real(PMFDP) function abp_accumulator_bin_pos(idxcoord,bindex)

    use abp_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    integer        :: bindex
    ! --------------------------------------------------------------------------

    ! 0.5 factor is for the middle of the bin
    abp_accumulator_bin_pos = accumulator%sizes(idxcoord)%min_value + &
                               accumulator%sizes(idxcoord)%bin_width*(real(bindex,PMFDP)+0.5d0)

end function abp_accumulator_bin_pos

!===============================================================================
! Subroutine:  abp_accumulator_get_values
!                       gindx: 1,2, ..., totalbins
!                       values: bin position
!===============================================================================

subroutine abp_accumulator_get_values(gindx,lvalues)

    use abp_dat

    implicit none
    integer        :: gindx
    real(PMFDP)    :: lvalues(:)
    ! -----------------------------------------------
    integer        :: i,indx,bindx
    ! --------------------------------------------------------------------------

    indx = gindx - 1

    do i=NumOfABPCVs,1
        bindx = MOD(indx,accumulator%sizes(i)%nbins)
        indx  = indx / accumulator%sizes(i)%nbins
        lvalues(i) = abp_accumulator_bin_pos(i,bindx)
    end do

end subroutine abp_accumulator_get_values

!===============================================================================
! Subroutine:  abp_accumulator_read
!===============================================================================

subroutine abp_accumulator_read(iounit)

    use abp_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer            :: iounit
    ! -----------------------------------------------
    character(len=4)                :: sabp
    character(len=2)                :: sver
    character(len=PMF_MAX_TYPE)     :: stype,sunit
    character(len=PMF_MAX_CV_NAME)  :: sname
    integer                         :: i,j,ibitem,it,nbins
    real(PMFDP)                     :: min_value,max_value,fconv,falpha
    ! --------------------------------------------------------------------------

    ! read header --------------------------
    read(iounit,10,end=100,err=100) sabp, sver, ibitem

    if( sabp .ne. ' ABP' ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABP] Attempt to read non-ABP accumulator!')
    end if

    if( sver .ne. 'V1' ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABP] Attempt to read ABP accumulator with unsupported version!')
    end if

    if( ibitem .ne. NumOfABPCVs ) then
        call pmf_utils_exit(PMF_OUT,1,'[ABP] ABP accumulator contains different number of coordinates!')
    end if

! read header --------------------------
    do i=1, NumOfABPCVs
        ! read CV definition
        read(iounit,20,end=102,err=102) it, stype, min_value, max_value, nbins
        ! check CV definition
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Incorrect item in ABP accumulator!')
        end if
        if( trim(adjustl(stype)) .ne. trim(ABPCVList(i)%cv%ctype) ) then
            write(PMF_OUT,*) '[ABP] CV type = [',trim(adjustl(stype)),'] should be [',trim(ABPCVList(i)%cv%ctype),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABP] CV type was redefined in ABP accumulator!')
        end if
        if( abs(min_value-accumulator%sizes(i)%min_value) .gt. abs(accumulator%sizes(i)%min_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Minimal value of CV was redefined in ABP accumulator!')
        end if
        if( abs(max_value-accumulator%sizes(i)%max_value) .gt. abs(accumulator%sizes(i)%max_value/100000.0d0) ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Maximal value of CV was redefined in ABP accumulator!')
        end if
        if( nbins .ne. accumulator%sizes(i)%nbins ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Number of CV bins was redefined in ABP accumulator!')
        end if

        ! read names
        read(iounit,25,end=102,err=102) it, sname
        ! check names
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Incorrect item in ABP accumulator!')
        end if
        if( trim(adjustl(sname)) .ne. trim(ABPCVList(i)%cv%name) ) then
            write(PMF_OUT,*) '[ABP] CV name = [',trim(adjustl(sname)),'] should be [',trim(ABPCVList(i)%cv%name),']'
            call pmf_utils_exit(PMF_OUT,1,'[ABP] CV name was redefined in ABP accumulator!')
        end if

        ! read names
        read(iounit,26,end=102,err=102) it, fconv, falpha, sunit
        ! check names
        if( it .ne. i ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABP] Incorrect item in ABP accumulator!')
        end if
        ! ignore values fconv, falpha and sunit
    end do

    ! read accumulator data - ABP forces
    if( accumulator%tot_nbins .gt. 0 ) then
        read(iounit,30,end=100,err=100) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)
        do i=1,NumOfABPCVs
            read(iounit,40,end=100,err=100) (accumulator%dpop(i,j),j=1,accumulator%tot_nbins)
        end do
        read(iounit,40,end=100,err=100) (accumulator%pop(i),i=1,accumulator%tot_nbins)
    end if

    return

10  format(A4,1X,A2,1X,I3)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,E20.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[ABP] Unable to read from ABP accumulator!')
102 call pmf_utils_exit(PMF_OUT,1,'[ABF] Unable to read from ABP accumulator v1 - CV section!')

end subroutine abp_accumulator_read

!===============================================================================
! Subroutine:  abp_accumulator_write
!===============================================================================

subroutine abp_accumulator_write(iounit)

    use abp_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer                :: iounit
    ! -----------------------------------------------
    integer                :: i,j,ipotene,ingrps
    !---------------------------------------------------------------------------

    ipotene = 0
    ingrps = 0

    ! write header
    write(iounit,10) 'ABP', 'V1', NumOfABPCVs
    do i=1, NumOfABPCVs
        write(iounit,20) i,trim(ABPCVList(i)%cv%ctype), &
                          accumulator%sizes(i)%min_value,accumulator%sizes(i)%max_value, &
                          accumulator%sizes(i)%nbins
        write(iounit,25) i,trim(ABPCVList(i)%cv%name)
        write(iounit,26) i,pmf_unit_get_rvalue(ABPCVList(i)%cv%unit,1.0d0),ABPCVList(i)%alpha, &
                         trim(pmf_unit_label(ABPCVList(i)%cv%unit))

    end do

    ! write accumulator data
    write(iounit,30) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)
    do i=1,NumOfABPCVs
        write(iounit,40) (accumulator%dpop(i,j),j=1,accumulator%tot_nbins)
    end do
    write(iounit,40) (accumulator%pop(i),i=1,accumulator%tot_nbins)

    return

10  format(A4,1X,A2,1X,I3,1X,I3,1X,I3)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,E20.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

end subroutine abp_accumulator_write

!===============================================================================
! Subroutine:  abp_accumulator_get_la_hramp
!===============================================================================

subroutine abp_accumulator_get_la_hramp

    use abp_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer        :: gi0
    real(PMFDP)    :: sc_ramp
    ! --------------------------------------------------------------------------

    la(:) = 0.0d0

    ! get global index to accumulator for current CV values
    gi0 = abp_accumulator_globalindex(cvvalues)

    if( gi0 .gt. 0 ) then
        sc_ramp = min(1.0d0, real(accumulator%pop(gi0),PMFDP) / real(fhramp,PMFDP) )
        la(:) = sc_ramp * accumulator%dpop(:,gi0) / accumulator%pop(gi0)
    end if

end subroutine abp_accumulator_get_la_hramp

!===============================================================================
! Subroutine:  abp_accumulator_update_direct
!===============================================================================

subroutine abp_accumulator_update_direct

    use abp_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer     :: gi0,gi,i
    real(PMFDP) :: d,w,s,ew
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for current CV values
    gi0 = abp_accumulator_globalindex(cvvalues)
    if( gi0 .le. 0 ) return   ! out of range

    accumulator%nsamples(gi0) = accumulator%nsamples(gi0) + 1
    if( fserver_enabled ) then
        accumulator%nisamples(gi0) = accumulator%nisamples(gi0) + 1
    end if

    ! calculate W factor
    w = exp( cfac )*accumulator%pop(gi0)/accumulator%M

    ! direct mode - for each bin
    do gi=1,accumulator%tot_nbins
        ! update pop
        s = 0.0d0
        do i=1,NumOfABPCVs
            ! calculate difference considering periodicity
            diffvalues(i) = ABPCVList(i)%cv%get_deviation(accumulator%binpos(i,gi),cvvalues(i))
            ! argument for exp
            s = s + diffvalues(i)**2 / ABPCVList(i)%alpha**2
        end do

        ew = exp(-s)*w

        accumulator%pop(gi) = accumulator%pop(gi) + ew
        if( accumulator%pop(gi) .gt. accumulator%M ) then
            accumulator%M = accumulator%pop(gi)
        end if

        if( fserver_enabled ) then
            accumulator%ipop(gi) = accumulator%ipop(gi) + ew
        end if

        ! update dpop
        do i=1,NumOfABPCVs
            d = 2.0d0 * diffvalues(i) / ABPCVList(i)%alpha**2
            accumulator%dpop(i,gi) = accumulator%dpop(i,gi) + d*ew
            if( fserver_enabled ) then
                accumulator%idpop(i,gi) = accumulator%idpop(i,gi) + d*ew
            end if
        end do

    end do

end subroutine abp_accumulator_update_direct

!===============================================================================

end module abp_accumulator

