!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtd_accumulator

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_accumulator_init
!===============================================================================

subroutine mtd_accumulator_init()

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    integer              :: j,k,multibins
    ! --------------------------------------------------------------------------

    ! init dimensions ------------------------------
    allocate(accumulator%sizes(NumOfMTDCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd accumulator!')
    endif

    tot_nbins = 1
    do i=1,NumOfMTDCVs
        accumulator%sizes(i)%min_value  = MTDCVList(i)%min_value
        accumulator%sizes(i)%max_value  = MTDCVList(i)%max_value
        accumulator%sizes(i)%nbins      = MTDCVList(i)%nbins
        accumulator%sizes(i)%width      = abs(accumulator%sizes(i)%max_value - accumulator%sizes(i)%min_value)
        accumulator%sizes(i)%bin_width  = accumulator%sizes(i)%width / accumulator%sizes(i)%nbins
        tot_nbins = tot_nbins * MTDCVList(i)%nbins
    end do

    accumulator%tot_nbins = tot_nbins

    ! MTD potential and force array
    allocate( accumulator%nsamples(accumulator%tot_nbins), &
              accumulator%binpositions(NumOfMTDCVs,accumulator%tot_nbins), &
              accumulator%mtdpotential(accumulator%tot_nbins), &
              accumulator%mtdforce(NumOfMTDCVs,accumulator%tot_nbins), &
              stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd accumulator!')
    endif

    do i=1,NumOfMTDCVs
        multibins = 1

        do j=i+1,NumOfMTDCVs
            multibins = multibins * accumulator%sizes(j)%nbins
        end do

        do j=1,accumulator%tot_nbins
            k = mod(((j-1)/multibins),accumulator%sizes(i)%nbins)
            accumulator%binpositions(i,j) = accumulator%sizes(i)%min_value + &
                                            real(k)*accumulator%sizes(i)%bin_width + accumulator%sizes(i)%bin_width / 2.0d0
        end do
    end do

    call mtd_accumulator_clear()

    return

end subroutine mtd_accumulator_init

!===============================================================================
! Subroutine:  mtd_accumulator_clear
!===============================================================================

subroutine mtd_accumulator_clear()

    use mtd_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    accumulator%nsamples(:)     = 0.0d0
    accumulator%mtdpotential(:) = 0.0d0
    accumulator%mtdforce(:,:)   = 0.0d0

end subroutine mtd_accumulator_clear

!===============================================================================
! Function:  grid_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               gridvalue    ... value that is used to compute the bin index
! compute index for one grid coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function mtd_accumulator_index(idxcoord,gridvalue)

    use mtd_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: gridvalue
    ! --------------------------------------------------------------------------

    if( accumulator%sizes(idxcoord)%bin_width .eq. 0.0d0 ) then
        mtd_accumulator_index = -1
        return
    end if

    ! we need number from zero - therefore we use floor(x)
    mtd_accumulator_index = floor((gridvalue - accumulator%sizes(idxcoord)%min_value) / &
                               accumulator%sizes(idxcoord)%bin_width)

    if( mtd_accumulator_index .lt. 0 .or. mtd_accumulator_index .ge. accumulator%sizes(idxcoord)%nbins) then
        mtd_accumulator_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function mtd_accumulator_index

!===============================================================================
! Function:  grid_globalindex
! Description:  Compute globalindex for grid, based on gridvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function mtd_accumulator_globalindex(lvalues)

    use mtd_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    mtd_accumulator_globalindex = 0

    do i=1,NumOfMTDCVs
        idx_local = mtd_accumulator_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            mtd_accumulator_globalindex = -1
            return
        end if

        mtd_accumulator_globalindex = mtd_accumulator_globalindex*accumulator%sizes(i)%nbins + idx_local
    end do

    mtd_accumulator_globalindex = mtd_accumulator_globalindex + 1

    return

end function mtd_accumulator_globalindex

!===============================================================================
! Subroutine:  mtd_accumulator_read
!===============================================================================

subroutine mtd_accumulator_read(iounit)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer             :: iounit
    ! -----------------------------------------------
    integer             :: ibitem
    character(len=4)    :: smtd
    character(len=2)    :: sver
    character(len=80)   :: buffer
    ! --------------------------------------------------------------------------

    read(iounit,'(A80)') buffer

    ! read header --------------------------
    read(buffer,*,end=200,err=200) smtd, sver, ibitem

    if( trim(adjustl(smtd)) .ne. 'MTD' ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] Missing MTD key in MTD accumulator header!')
    end if

    if( ibitem .ne. NumOfMTDCVs ) then
        call pmf_utils_exit(PMF_OUT,1,'[MTD] MTD accumulator contains different number of CVs!')
    end if

    if( trim(adjustl(sver)) .eq. 'V6' ) then
        call mtd_accumulator_read_v6(iounit)
        return
    end if

    write(PMF_OUT,*) '[MTD] header ver = [',sver,']'
    call pmf_utils_exit(PMF_OUT,1,'[MTD] Unsupported version key in MTD accumulator header!')

    return

    ! unsupported header - process error
200 continue
    call pmf_utils_exit(PMF_OUT,1,'[MTD] Illegal header in MTD accumulator!')

end subroutine mtd_accumulator_read

!===============================================================================
! Subroutine:  mtd_accumulator_read_v6
!===============================================================================

subroutine mtd_accumulator_read_v6(iounit)

    use mtd_dat
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
        read(iounit,5,end=1000,err=100) key

        select case( trim(key) )
            case('CVS')
                ! read header --------------------------
                do i=1, NumOfMTDCVs
                    ! read CV definition
                    read(iounit,20,end=101,err=101) it, stype, min_value, max_value, nbins
                    ! check CV definition
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] Incorrect item in MTD accumulator v4!')
                    end if
                    if( trim(adjustl(stype)) .ne. trim(MTDCVList(i)%cv%ctype) ) then
                        write(PMF_OUT,*) '[MTD] CV type = [',trim(adjustl(stype)),'] should be [',trim(MTDCVList(i)%cv%ctype),']'
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] CV type was redefined in MTD accumulator!')
                    end if
                    if( abs(min_value-accumulator%sizes(i)%min_value) .gt. abs(accumulator%sizes(i)%min_value/100000.0d0) ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] Minimal value of CV was redefined in MTD accumulator!')
                    end if
                    if( abs(max_value-accumulator%sizes(i)%max_value) .gt. abs(accumulator%sizes(i)%max_value/100000.0d0) ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] Maximal value of CV was redefined in MTD accumulator!')
                    end if
                    if( nbins .ne. accumulator%sizes(i)%nbins ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] Number of CV bins was redefined in MTD accumulator!')
                    end if

                    ! read names
                    read(iounit,25,end=102,err=102) it, sname
                    ! check names
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] Incorrect item in MTD accumulator!')
                    end if
                    if( trim(adjustl(sname)) .ne. trim(MTDCVList(i)%cv%name) ) then
                        write(PMF_OUT,*) '[MTD] CV name = [',trim(adjustl(sname)),'] should be [',trim(MTDCVList(i)%cv%name),']'
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] CV name was redefined in MTD accumulator!')
                    end if

                    ! read names
                    read(iounit,26,end=103,err=103) it, fconv, sunit
                    ! check names
                    if( it .ne. i ) then
                        call pmf_utils_exit(PMF_OUT,1,'[MTD] Incorrect item in MTD accumulator!')
                    end if
                    ! ignore values fconv and sunit
                end do
            case('TEMPERATURE')
                ! read but do not use
                read(iounit,6,end=201,err=201) fconv
            case('ENERGY-UNIT')
                ! read but do not use
                read(iounit,27,end=202,err=202) fconv, sunit

            case('NSAMPLES')
                read(iounit,30,end=500,err=500) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)
            case('MTDPOT')
                read(iounit,30,end=600,err=600) (accumulator%mtdpotential(i),i=1,accumulator%tot_nbins)
            case('MTDFORCE')
                do i=1,NumOfMTDCVs
                    read(iounit,40,end=700,err=700) (accumulator%mtdforce(i,j),j=1,accumulator%tot_nbins)
                end do
            case default
                call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - unrecognized keyword: '//trim(key))
        end select
    end do

1000 return

 5  format(A20)
 6  format(F10.4)
20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
25  format(I2,1X,A55)
26  format(I2,1X,E18.11,1X,A36)
27  format(E18.11,1X,A36)
30  format(8(I9,1X))
40  format(4(E19.11,1X))

100 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - keyword!')
101 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - CV section - def!')
102 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - CV section - name!')
103 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - CV section - unit!')

201 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - temperature!')
202 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - energy unit!')

500 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - data section - bins!')
600 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - data section - smtd!')
700 call pmf_utils_exit(PMF_OUT,1,'[MTD] Unable to read from MTD accumulator v6 - data section - smtd2!')

end subroutine mtd_accumulator_read_v6

!===============================================================================
! Subroutine:  mtd_accumulator_write
!===============================================================================

subroutine mtd_accumulator_write(iounit)

    use mtd_dat
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
    write(iounit,10) 'MTD ', 'V6', NumOfMTDCVs

    key = 'CVS'
    write(iounit,5) adjustl(key)
    do i=1, NumOfMTDCVs
        write(iounit,20) i,trim(MTDCVList(i)%cv%ctype), &
                          accumulator%sizes(i)%min_value,accumulator%sizes(i)%max_value, &
                          accumulator%sizes(i)%nbins
        write(iounit,25) i,trim(MTDCVList(i)%cv%name)
        write(iounit,26) i,pmf_unit_get_rvalue(MTDCVList(i)%cv%unit,1.0d0),trim(pmf_unit_label(MTDCVList(i)%cv%unit))
    end do

    key = 'TEMPERATURE'
    write(iounit,5) adjustl(key)
    write(iounit,6) ftemp

    key = 'ENERGY-UNIT'
    write(iounit,5) adjustl(key)
    write(iounit,27) pmf_unit_get_rvalue(EnergyUnit,1.0d0),trim(pmf_unit_label(EnergyUnit))

    key = 'NSAMPLES'
    write(iounit,5) adjustl(key)
    write(iounit,30) (accumulator%nsamples(i),i=1,accumulator%tot_nbins)

    key = 'MTDPOT'
    write(iounit,5) adjustl(key)
    write(iounit,40) (accumulator%mtdpotential(i),i=1,accumulator%tot_nbins)

    key = 'MTDFORCE'
    write(iounit,5) adjustl(key)
    do i=1,NumOfMTDCVs
        write(iounit,40) (accumulator%mtdforce(i,j),j=1,accumulator%tot_nbins)
    end do

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

end subroutine mtd_accumulator_write

!===============================================================================
! Subroutine:  mtd_accumulator_add_data
!===============================================================================

subroutine mtd_accumulator_add_data(cvs, height, widths)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: height
    real(PMFDP)    :: widths(:)
    ! -----------------------------------------------
    integer        :: n, i, gi0
    real(PMFDP)    :: diff, fexparg, fh
    ! --------------------------------------------------------------------------

    ! get global index to accumulator for average values within the set
    gi0 = mtd_accumulator_globalindex(cvs)
    if( gi0 .le. 0 ) then
        outsidesamples = outsidesamples + 1
    else
        insidesamples = insidesamples + 1
        ! increase number of samples
        accumulator%nsamples(gi0) = accumulator%nsamples(gi0) + 1
    end if

    ! update grid data
    do n=1,accumulator%tot_nbins
        fexparg = 0.0d0
        do i=1,NumOfMTDCVs
            diff = MTDCVList(i)%cv%get_deviation(accumulator%binpositions(i,n), cvs(i))
            fexparg = fexparg + diff**2 / (2.0d0 * widths(i)**2)
        end do
        fh = height * exp(-fexparg)
        accumulator%mtdpotential(n) = accumulator%mtdpotential(n) + fh
        do i=1,NumOfMTDCVs
            diff = MTDCVList(i)%cv%get_deviation(accumulator%binpositions(i,n), cvs(i))
            accumulator%mtdforce(i,n) = accumulator%mtdforce(i,n) + fh * diff / widths(i)**2
        end do
    end do

    return

end subroutine mtd_accumulator_add_data

!===============================================================================
! Subroutine:  mtd_accumulator_get_data
!===============================================================================

subroutine mtd_accumulator_get_data(cvs,potential,forces)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: cvs(:)
    real(PMFDP)    :: potential
    real(PMFDP)    :: forces(:)
    ! -----------------------------------------------
    integer        :: gi0
    ! --------------------------------------------------------------------------

    potential = 0.0d0
    forces(:) = 0.0d0

    ! get global index to grid for cvs
    gi0 = mtd_accumulator_globalindex(cvs)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get potential
    potential = accumulator%mtdpotential(gi0)

    ! get forces
    forces(:) = accumulator%mtdforce(:,gi0)

end subroutine mtd_accumulator_get_data

!===============================================================================

end module mtd_accumulator
