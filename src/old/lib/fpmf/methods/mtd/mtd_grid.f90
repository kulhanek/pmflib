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

module mtd_grid

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  mtd_grid_init
!===============================================================================

subroutine mtd_grid_init()

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    integer              :: i,tot_nbins
    integer              :: alloc_failed
    integer              :: j,k,multibins
    ! --------------------------------------------------------------------------

    if (fmode .ne. 2) return

    ! init dimensions ------------------------------
    allocate(grid%sizes(NumOfMTDCVs), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd grid!')
    endif

    tot_nbins = 1
    do i=1,NumOfMTDCVs
        grid%sizes(i)%min_value  = MTDCVList(i)%min_value
        grid%sizes(i)%max_value  = MTDCVList(i)%max_value
        grid%sizes(i)%nbins      = MTDCVList(i)%nbins
        grid%sizes(i)%width      = abs(grid%sizes(i)%max_value - grid%sizes(i)%min_value)
        grid%sizes(i)%bin_width  = grid%sizes(i)%width / grid%sizes(i)%nbins
        tot_nbins = tot_nbins * MTDCVList(i)%nbins
    end do

    grid%tot_nbins = tot_nbins

    ! MTD potential and force array
    allocate( grid%binpositions(NumOfMTDCVs,grid%tot_nbins), &
              grid%mtdpotential(grid%tot_nbins), &
              grid%mtdforce(NumOfMTDCVs,grid%tot_nbins), &
              stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT, 1,'[MTD] Unable to allocate memory for mtd grid!')
    endif

    do i=1,NumOfMTDCVs
        multibins = 1

        do j=i+1,NumOfMTDCVs
            multibins = multibins * grid%sizes(j)%nbins
        end do

        do j=1,grid%tot_nbins
            k = mod(((j-1)/multibins),grid%sizes(i)%nbins)
            grid%binpositions(i,j) = grid%sizes(i)%min_value + &
                                            real(k)*grid%sizes(i)%bin_width + grid%sizes(i)%bin_width / 2.0d0
        end do
    end do

    call mtd_grid_clear()

    return

end subroutine mtd_grid_init

!===============================================================================
! Subroutine:  mtd_grid_clear
!===============================================================================

subroutine mtd_grid_clear()

    use mtd_dat
    use pmf_dat

    implicit none
    ! --------------------------------------------------------------------------

    if (fmode .ne. 2) return

    grid%mtdpotential(:) = 0.0d0
    grid%mtdforce(:,:) = 0.0d0

end subroutine mtd_grid_clear

!===============================================================================
! Function:  grid_index
! Arguments:
!               idxcoord ... number of ksi coordinate
!               gridvalue    ... value that is used to compute the bin index 
! compute index for one grid coordinate
! Return value:     0,1,2, ..., sizes(idxcoord)%numbins-1
!===============================================================================

integer function mtd_grid_index(idxcoord,gridvalue)

    use mtd_dat
    use pmf_dat

    implicit none
    integer        :: idxcoord
    real(PMFDP)    :: gridvalue
    ! --------------------------------------------------------------------------

    if( grid%sizes(idxcoord)%bin_width .eq. 0.0d0 ) then
        mtd_grid_index = -1
        return
    end if

    ! we need number from zero - therefore we use floor(x)
    mtd_grid_index = floor((gridvalue - grid%sizes(idxcoord)%min_value) / &
                               grid%sizes(idxcoord)%bin_width)

    if( mtd_grid_index .lt. 0 .or. mtd_grid_index .ge. grid%sizes(idxcoord)%nbins) then
        mtd_grid_index = -1
        return
    end if

    ! do not try to include right boundary, since it will include the whole additional bin !

    return

end function mtd_grid_index

!===============================================================================
! Function:  grid_globalindex
! Description:  Compute globalindex for grid, based on gridvalues of all coordinates
! Arguments:    none
! Return value: 1,2, ..., totalbins
!===============================================================================

integer function mtd_grid_globalindex(lvalues)

    use mtd_dat
    use pmf_dat

    implicit none
    real(PMFDP)            :: lvalues(:)
    ! -----------------------------------------------
    integer                :: idx_local,i
    ! --------------------------------------------------------------------------

    mtd_grid_globalindex = 0

    do i=1,NumOfMTDCVs
        idx_local = mtd_grid_index(i,lvalues(i))

        if (idx_local .eq. -1) then
            mtd_grid_globalindex = -1
            return
        end if

        mtd_grid_globalindex = mtd_grid_globalindex*grid%sizes(i)%nbins + idx_local
    end do

    mtd_grid_globalindex = mtd_grid_globalindex + 1

    return

end function mtd_grid_globalindex

!===============================================================================
! Subroutine:  mtd_grid_add_data
!===============================================================================

subroutine mtd_grid_add_data(values, height, widths)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: height
    real(PMFDP)    :: widths(:)
    ! -----------------------------------------------
    integer        :: n, i
    real(PMFDP)    :: diff, fexparg, fh
    ! --------------------------------------------------------------------------

    do n=1,grid%tot_nbins
       fexparg = 0.0d0
       do i=1,NumOfMTDCVs
           diff = MTDCVList(i)%cv%get_deviation(grid%binpositions(i,n), values(i))
           fexparg = fexparg + diff**2 / (2.0d0 * widths(i)**2)
       end do
       fh = height * exp(-fexparg)
       grid%mtdpotential(n) = grid%mtdpotential(n) + fh
       do i=1,NumOfMTDCVs
           diff = MTDCVList(i)%cv%get_deviation(grid%binpositions(i,n), values(i))
           grid%mtdforce(i,n) = grid%mtdforce(i,n) + fh * diff / widths(i)**2
       end do
    end do

    return

end subroutine mtd_grid_add_data

!===============================================================================
! Subroutine:  mtd_grid_get_data
!===============================================================================

subroutine mtd_grid_get_data(values,potential,forces)

    use mtd_dat
    use pmf_dat
    use pmf_utils

    implicit none
    real(PMFDP)    :: values(:)
    real(PMFDP)    :: potential
    real(PMFDP)    :: forces(:)
    ! -----------------------------------------------
    integer        :: gi0
    ! --------------------------------------------------------------------------

    potential = 0.0d0
    forces(:) = 0.0d0

    ! get global index to grid for average values within the set
    gi0 = mtd_grid_globalindex(values)
    if( gi0 .le. 0 ) return ! out of valid area

    ! get potential
    potential = grid%mtdpotential(gi0)

    ! get forces
    forces(:) = grid%mtdforce(:,gi0)

end subroutine mtd_grid_get_data

!===============================================================================

end module mtd_grid
