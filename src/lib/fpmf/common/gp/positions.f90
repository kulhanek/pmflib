!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module positions_mod

implicit none
private

integer, parameter :: dp = 8
type sort_array
    real(dp) :: val
    integer :: indx
end type

public :: positions_wrap
public :: positions_vec
public :: positions_dist
public :: positions_find_first_n_neighbors
public :: positions_find_min_max
public :: positions_entry_min_max

interface positions_wrap
   module procedure positions_wrap_1d, positions_wrap_nd
end interface

contains

!--------------------------------------------------------------------------------------------------
   
subroutine positions_wrap_1d(positions, periodicities)

  implicit none

  real(dp), intent(inout) :: positions(:)
  real(dp), intent(in) :: periodicities(:)

  integer :: ndim
  integer :: i

  ! check if sizes are compatible
! if ( size(positions) .ne. size(periodicities) ) then
!      write(6,*) '>>> ERROR: [positions_wrap_single] ndim in positions and periodicities differs!'
!      stop
! end if

  ndim = size(positions)

  do i=1,ndim
     if ( periodicities(i) .ne. 0 ) then
          positions(i) = positions(i) - nint(positions(i)/periodicities(i))*periodicities(i)
     end if
  end do

end subroutine positions_wrap_1d

!--------------------------------------------------------------------------------------------------

subroutine positions_wrap_nd(positions, periodicities)

  implicit none

  real(dp), intent(inout) :: positions(:,:)
  real(dp), intent(in) :: periodicities(:)

  integer :: ndim, nsample
  integer :: i, j

  ! check if sizes are compatible
! if ( size(positions, 1) .ne. size(periodicities) ) then
!      write(6,'(a)') '>>> ERROR: [positions_wrap_array] ndim in positions and periodicities differs!'
!      stop
! end if

  ndim = size(positions, 1)
  nsample = size(positions, 2)

  do j=1,nsample
     do i=1,ndim
        if ( periodicities(i) .ne. 0 ) then
             positions(i,j) = positions(i,j) - nint(positions(i,j)/periodicities(i))*periodicities(i)
        end if
     end do
  end do

end subroutine positions_wrap_nd

!--------------------------------------------------------------------------------------------------

function positions_vec(positions1, positions2, periodicities)

  implicit none

  real(dp), intent(in) :: positions1(:), positions2(:)
  real(dp), intent(in) :: periodicities(:)

  real(dp) :: positions_vec(size(positions1))

  integer :: ndim
  integer :: i
  real(dp) :: d

  ! check if size are compatible
! if ( size(positions1) .ne. size(positions2) ) then
!      write(6,'(a)') '>>> ERROR: [positions_dist] ndim in positions1 and positions2 differs!'
!      stop
! end if

! if ( size(positions1) .ne. size(periodicities) ) then
!      write(6,'(a)') '>>> ERROR: [positions_dist] ndim in positions and periodicities differs!'
!      stop
! end if

  ndim = size(positions1)

  do i=1,ndim
     d = positions1(i) - positions2(i)
     if ( periodicities(i) .ne. 0 ) then
          d = d - nint(d/periodicities(i))*periodicities(i)
     end if
     positions_vec(i) = d
  end do

end function positions_vec

!--------------------------------------------------------------------------------------------------

function positions_dist(positions1, positions2, periodicities)

  implicit none

  real(dp), intent(in) :: positions1(:), positions2(:)
  real(dp), intent(in) :: periodicities(:)

  real(dp) :: positions_dist

  integer :: ndim
  integer :: i
  real(dp) :: d

  ! check if size are compatible
! if ( size(positions1) .ne. size(positions2) ) then
!      write(6,'(a)') '>>> ERROR: [positions_dist] ndim in positions1 and positions2 differs!'
!      stop
! end if

! if ( size(positions1) .ne. size(periodicities) ) then
!      write(6,'(a)') '>>> ERROR: [positions_dist] ndim in positions and periodicities differs!'
!      stop
! end if

  ndim = size(positions1)

  positions_dist = 0.0_dp
  do i=1,ndim
     d = positions1(i) - positions2(i)
     if ( periodicities(i) .ne. 0 ) then
          d = d - nint(d/periodicities(i))*periodicities(i)
     end if
     positions_dist = positions_dist + d*d
  end do

  positions_dist = sqrt(positions_dist)

end function positions_dist

!--------------------------------------------------------------------------------------------------

subroutine positions_find_first_n_neighbors(positions, periodicities, neighbors)

  implicit none

  real(dp), intent(in) :: positions(:,:)
  real(dp), intent(in) :: periodicities(:)

  integer, intent(inout) :: neighbors(:,:)

  integer :: ndim, nsample, nneigh
  integer :: i, j, k
  type(sort_array), allocatable :: array(:)

  ! check if size are compatible
! if ( size(positions,1) .ne. size(periodicities) ) then
!      write(6,'(a)') '>>> ERROR: [positions_find_neighbors] ndim in positions and periodicities differs!'
!      stop
! end if

! if ( size(positions,2) .ne. size(neighbors,2) ) then
!      write(6,'(a)') '>>> ERROR: [positions_find_neighbors] nsamples in positions and neighbors differs!'
!      stop
! end if

  ndim = size(positions, 1)
  nsample = size(positions, 2)
  nneigh = size(neighbors, 1)

  allocate(array(nsample-1))

  do i=1,nsample
     k = 1
     do j=1,nsample
        if(i .eq. j) cycle
        array(k)%val = positions_dist(positions(:,i), positions(:,j), periodicities(:))
        array(k)%indx = j
        k = k + 1
     end do
     call quickselect_left(array, nneigh)
     neighbors(:,i) = array(1:nneigh)%indx
  end do

end subroutine positions_find_first_n_neighbors

!--------------------------------------------------------------------------------------------------

subroutine positions_find_min_max(positions, periodicities, min_positions, max_positions, range_positions)

  implicit none

  real(dp), intent(in) :: positions(:,:)
  real(dp), intent(in) :: periodicities(:)
  real(dp), intent(out) :: min_positions(:)
  real(dp), intent(out) :: max_positions(:)
  real(dp), intent(out) :: range_positions(:)

  integer :: nsample
  integer :: i

  ! check if size are compatible
  if ( size(positions, 1) .ne. size(periodicities) ) then
       write(6,'(a)') '>>> ERROR: [positions_find_min_max] ndim in positions and periodicities differs!'
       stop
  end if

  if ( size(positions, 1) .ne. size(min_positions) ) then
       write(6,'(a)') '>>> ERROR: [positions_find_min_max] ndim in positions and min_positions differs!'
       stop
  end if

  if ( size(positions, 1) .ne. size(max_positions) ) then
       write(6,'(a)') '>>> ERROR: [positions_find_min_max] ndim in positions and max_positions differs!'
       stop
  end if

  if ( size(positions, 1) .ne. size(range_positions) ) then
       write(6,'(a)') '>>> ERROR: [positions_find_min_max] ndim in positions and range_positions differs!'
       stop
  end if

  nsample = size(positions, 2)

  min_positions(:) = positions(:,1)
  max_positions(:) = positions(:,1)
  range_positions(:) = 0.0d0

  do i=2,nsample
     call positions_entry_min_max(positions(:,i), periodicities, min_positions, max_positions, range_positions)
  end do

end subroutine positions_find_min_max

!--------------------------------------------------------------------------------------------------

subroutine positions_entry_min_max(pos, per, min_pos, max_pos, range_pos)

  implicit none

  real(dp), intent(in) :: pos(:)
  real(dp), intent(in) :: per(:)
  real(dp), intent(inout) :: min_pos(:)
  real(dp), intent(inout) :: max_pos(:)
  real(dp), intent(inout) :: range_pos(:)

  integer :: ndim
  integer :: i
  real(dp) :: pval
  real(dp) :: diffmin, diffmax

  ndim=size(pos)

  do i=1,ndim
     ! cv is not periodic
     if ( per(i) .eq. 0 ) then
          if ( pos(i) < min_pos(i) ) then
               min_pos(i) = pos(i)
               range_pos(i) = max_pos(i) - min_pos(i)
          else if( pos(i) > max_pos(i) ) then
               max_pos(i) = pos(i)
               range_pos(i) = max_pos(i) - min_pos(i)
          end if
     ! cv is periodic
     else if ( range_pos(i) .ne. per(i) ) then
          ! move value into default range
          pval = pos(i) - nint(pos(i)/per(i))*per(i)
          diffmin = abs(pval - min_pos(i))
          diffmax = abs(pval - max_pos(i))
          ! criterion of new min/max
          ! min < max
          if ( min_pos(i) < max_pos(i) ) then
               if ( (pval < min_pos(i)) .or. (pval > max_pos(i)) ) then
                    if( diffmin < diffmax ) then
                        min_pos(i) = pval
                    else
                        max_pos(i) = pval
                    end if
               end if
          ! min > max
          else if ( min_pos(i) > max_pos(i) ) then
               if ( (pval < min_pos(i)) .and. (pval > max_pos(i)) ) then
                    if( diffmin < diffmax ) then
                        min_pos(i) = pval
                    else
                        max_pos(i) = pval
                    end if
               end if
          ! min = max
          else
              ! min/max and val are close to each other
              if( diffmin .le. (per(i) / 2.0d0) ) then
                  if( pval < min_pos(i) ) then
                      min_pos(i) = pval
                  else
                      max_pos(i) = pval
                  end if
              ! min/max and val are far from each other
              else
                  ! min/max is close to right edge, val appears on left
                  if( pval < min_pos(i) ) then
                      max_pos(i) = pval
                  ! min/max is close to left edge, val appears on right
                  else
                      min_pos(i) = pval
                  end if
              end if
          end if
          range_pos(i) = max_pos(i) - min_pos(i)
          if ( range_pos(i) .lt. 0 ) then
               range_pos(i) = range_pos(i) + per(i)
          end if
     end if
  end do

end subroutine positions_entry_min_max

!--------------------------------------------------------------------------------------------------

recursive subroutine quicksort(array)

  implicit none

  type(sort_array), intent(inout) :: array(:)

  integer :: ndim
  integer :: marker
  real(dp) :: random
  type(sort_array) :: pivot
  type(sort_array) :: temp
  integer :: left, right

  ndim = size(array)

  if (ndim <= 1) return

  call random_number(random)
  pivot = array(int(random*real(ndim-1))+1)
  left = 0
  right = ndim+1

  do while (left < right)
      right = right - 1
      do while (array(right)%val > pivot%val)
          right = right - 1
      end do
      left = left + 1
      do while (array(left)%val < pivot%val)
          left = left + 1
      end do
      if (left < right) then
         temp = array(left)
         array(left) = array(right)
         array(right) = temp
      end if
  end do

  if (left == right) then
     marker = left + 1
  else
     marker = left
  end if

  call quicksort(array(:marker-1))
  call quicksort(array(marker:))

end subroutine quicksort

!--------------------------------------------------------------------------------------------------

recursive subroutine quickselect_left(array, n)

  implicit none

  type(sort_array), intent(inout) :: array(:)
  integer, intent(in) :: n

  integer :: ndim
  integer :: marker
  real(dp) :: random
  type(sort_array) :: pivot
  type(sort_array) :: temp
  integer :: left, right

  ndim = size(array)

  if (ndim <= 1) return

  call random_number(random)
  marker = int(random*real(ndim))+1
  pivot = array(marker)
  left = 1
  right = ndim
  
  array(marker) = array(right)
  array(right) = pivot

  marker = 1

  do while (left < right)
     if (array(left)%val < pivot%val) then
        temp = array(left)
        array(left) = array(marker)
        array(marker) = temp
        marker = marker + 1
     end if
     left = left + 1
  end do

  array(right) = array(marker)
  array(marker) = pivot

  if ( n == marker ) then
     return
  else if (n < marker) then
     call quickselect_left(array(:marker-1), n)
  else
     call quickselect_left(array(marker+1:), n)
  end if

end subroutine quickselect_left

!--------------------------------------------------------------------------------------------------

end module positions_mod
