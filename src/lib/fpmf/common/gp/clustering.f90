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

module clustering_mod

use positions_mod

implicit none

private

   integer, parameter :: dp = 8

   public :: random_clustering_pick, k_meanspp_clustering_pick, k_means_clustering_pick

   interface random_clustering_pick
      module procedure random_clustering_pick_1d, random_clustering_pick_nd
   end interface

   interface k_meanspp_clustering_pick
      module procedure k_meanspp_clustering_pick_1d, k_meanspp_clustering_pick_nd
   end interface

   interface k_means_clustering_pick
      module procedure k_means_clustering_pick_1d, k_means_clustering_pick_nd
   end interface

contains

!--------------------------------------------------------------------------------------------------

subroutine random_clustering_pick_1d(x, cluster_indices, preselected_indices)

   real(dp), intent(in) :: x(:)
   integer, intent(out) :: cluster_indices(:)
   integer, optional, intent(in) :: preselected_indices(:)

   integer :: n, l

   n = size(x)

   l = 0
   if(present(preselected_indices)) then
      l = size(preselected_indices)
   end if

   ! random initial guess
   if(l /= 0) then
      call random_selection(n, cluster_indices, preselected_indices)
   else
      call random_selection(n, cluster_indices)
   end if

end subroutine random_clustering_pick_1d

!--------------------------------------------------------------------------------------------------

subroutine random_clustering_pick_nd(x, cluster_indices, preselected_indices)

   real(dp), intent(in) :: x(:,:)
   integer, intent(out) :: cluster_indices(:)
   integer, optional, intent(in) :: preselected_indices(:)

   integer :: n, l

   n = size(x,2)

   l = 0
   if(present(preselected_indices)) then
      l = size(preselected_indices)
   end if

   ! random initial guess
   if(l /= 0) then
      call random_selection(n, cluster_indices, preselected_indices)
   else
      call random_selection(n, cluster_indices)
   end if

end subroutine random_clustering_pick_nd

!--------------------------------------------------------------------------------------------------

subroutine k_meanspp_clustering_pick_1d(x, periodicity, cluster_indices, preselected_indices)
   real(dp), intent(in) :: x(:)
   real(dp), intent(in) :: periodicity
   integer, intent(out) :: cluster_indices(:)
   integer, optional, intent(in) :: preselected_indices(:)

   integer :: i, j, k, l, m, n
   real(dp) :: r, p, rrv, rv
   integer, allocatable :: indices(:)
   real(dp), allocatable :: distances(:), probs(:)

   k = size(cluster_indices)
   n = size(x)

   l = 0
   if(present(preselected_indices)) then
      l = size(preselected_indices)
      if(l > k) then
         print *, "size of preselected_indices=", size(preselected_indices), &
                  " is larger than size of cluster_indices=", size(cluster_indices)
         stop
      end if
   end if

   cluster_indices = 0
   if(l /= 0) then
      do i=1, l
         cluster_indices(i) = preselected_indices(i)
      end do
      if(l == k) return
   else
      call random_number(rrv)
      rv = 1+floor(rrv*n)
      cluster_indices(1) = rv
      l = 1
   end if

   do while(l < k)
      if(.not. allocated(indices)) then
         allocate(indices(n-l))
         allocate(distances(n-l))
         allocate(probs(n-l))
         m = 1
         do i=1, n
            if(any(i == cluster_indices)) cycle
            indices(m) = i
            distances(m) = positions_dist((/x(i)/), (/x(cluster_indices(1))/), (/periodicity/))
            do j=2, l
               r = positions_dist((/x(i)/), (/x(cluster_indices(j))/), (/periodicity/))
               if(r > distances(m)) distances(m) = r
            end do
            probs(m) = distances(m)*distances(m)
            m = m + 1
            if(m > n-l) exit
        end do
      else
        do i=1, n-l
           r = positions_dist((/x(i)/), (/x(cluster_indices(l))/), (/periodicity/))
           if(r > distances(i)) then
              distances(i) = r
              probs(i) = distances(i)*distances(i)
           end if
        end do
      end if
      probs(1:n-l) = probs(1:n-l)/sum(probs(1:n-l))

      call random_number(rrv)
      p = 0.0_dp
      do i=1, n-l
         p = p + probs(i)
         if(i == n-l) p = 1.0_dp
         if(p > rrv) then
            cluster_indices(l+1) = indices(i)
            do j=i, n-l-1
               indices(j) = indices(j+1)
               distances(j) = distances(j+1)
            end do
            l = l + 1
            exit
         end if
      end do
   end do

   deallocate(indices)
   deallocate(distances)
   deallocate(probs)

end subroutine k_meanspp_clustering_pick_1d

!--------------------------------------------------------------------------------------------------

subroutine k_meanspp_clustering_pick_nd(x, periodicity, cluster_indices, preselected_indices)
   real(dp), intent(in) :: x(:,:)
   real(dp), intent(in) :: periodicity(:)
   integer, intent(out) :: cluster_indices(:)
   integer, optional, intent(in) :: preselected_indices(:)

   integer :: i, j, k, l, m, n
   integer :: nd
   real(dp) :: r, p, rrv, rv
   integer, allocatable :: indices(:)
   real(dp), allocatable :: distances(:), probs(:)

   k = size(cluster_indices)
   n = size(x,2)
   nd = size(x,1)

   l = 0
   if(present(preselected_indices)) then
      l = size(preselected_indices)
      if(l > k) then
         print *, "size of preselected_indices=", size(preselected_indices), &
                  " is larger than size of cluster_indices=", size(cluster_indices)
         stop
      end if
   end if

   cluster_indices = 0
   if(l /= 0) then
      do i=1, l
         cluster_indices(i) = preselected_indices(i)
      end do
      if(l == k) return
   else
      call random_number(rrv)
      rv = 1+floor(rrv*n)
      cluster_indices(1) = rv
      l = 1
   end if

   do while(l < k)
      if(.not. allocated(indices)) then
         allocate(indices(n-l))
         allocate(distances(n-l))
         allocate(probs(n-l))
         m = 1
         do i=1, n
            if(any(i == cluster_indices)) cycle
            indices(m) = i
            distances(m) = positions_dist(x(:,i), x(:,cluster_indices(1)), periodicity(:))
            do j=2, l
               r = positions_dist(x(:,i), x(:,cluster_indices(j)), periodicity(:))
               if(r > distances(m)) distances(m) = r
            end do
            probs(m) = distances(m)*distances(m)
            m = m + 1
            if(m > n-l) exit
        end do
      else
        do i=1, n-l
           r = positions_dist(x(:,indices(i)), x(:,cluster_indices(l)), periodicity(:))
           if(r > distances(i)) then
              distances(i) = r
              probs(i) = distances(i)*distances(i)
           end if
        end do
      end if
      probs(1:n-l) = probs(1:n-l)/sum(probs(1:n-l))

      call random_number(rrv)
      p = 0.0_dp
      do i=1, n-l
         p = p + probs(i)
         if(i == n-l) p = 1.0_dp
         if(p > rrv) then
            cluster_indices(l+1) = indices(i)
            do j=i, n-l-1
               indices(j) = indices(j+1)
               distances(j) = distances(j+1)
            end do
            l = l + 1
            exit
         end if
      end do
   end do

   deallocate(indices)
   deallocate(distances)
   deallocate(probs)

end subroutine k_meanspp_clustering_pick_nd

!--------------------------------------------------------------------------------------------------
subroutine k_means_clustering_pick_1d(x, periodicity, cluster_indices, starting_indices, thresholds, verbose)
   real(dp), intent(in) :: x(:)
   real(dp), intent(in) :: periodicity
   integer, intent(out) :: cluster_indices(:)
   integer, optional, intent(in) :: starting_indices(:)
   real(dp), optional, intent(in) :: thresholds
   logical, optional, intent(in) :: verbose

   real(dp) :: thr
   logical :: do_verbose
   integer :: k, l
   real(dp), allocatable :: k_means(:)
   integer :: n, rv, i, j, closest_i, closest_j
   real(dp) :: rrv
   real(dp) :: closest_r, r, t_sum
   integer, allocatable :: a(:), prev_a(:)
   real(dp) :: half_periodicity, shift
   integer :: iter

   if(present(verbose)) then
     do_verbose = verbose
   else
     do_verbose = .false.
   end if

   k = size(cluster_indices)
   n = size(x)
   allocate(k_means(k))

   ! calculate the half of the periodicity
   half_periodicity = periodicity / 2.0_dp

   l = 0
   if(present(starting_indices)) then
      l = size(starting_indices)
   end if

   thr = 0.0_dp
   if(present(thresholds)) then
      thr = thresholds
   end if

   ! random initial guess
   if(l /= 0) then
      call random_selection(n, cluster_indices, starting_indices)
   else
      call random_selection(n, cluster_indices)
   end if

   k_means = x(cluster_indices)

   allocate(a(n), prev_a(n))

   ! evaluate initial assignments
   do j=1, n
      closest_i = 0
      closest_r = HUGE(1.0_dp)
      do i=1, k
         r = positions_dist((/x(j)/), (/k_means(i)/),  (/periodicity/))
         if (r < closest_r) then
            closest_r = r
            closest_i = i
         endif
      end do
      a(j) = closest_i
   end do

   prev_a = 0
   iter=0
   do while (any(a /= prev_a))
      prev_a = a

      ! update positions
      do i=1, k
         t_sum = 0.0_dp
         shift = 0.0_dp
         if (abs(k_means(i)-x(j)) > half_periodicity) shift = sign(periodicity, k_means(i))
         do j=1, n
            if (a(j) == i) t_sum = t_sum + x(j) + shift
         end do
         t_sum = t_sum / real(count(a == i),dp)
         if( positions_dist((/k_means(i)/), (/t_sum/), (/periodicity/)) .gt. thr ) then
            k_means(i) = t_sum
            if(abs(k_means(i)) > half_periodicity) k_means(i) = k_means(i) - sign(periodicity, k_means(i))
         end if
      end do

      ! update assignments
      do j=1, n
         closest_i = 0
         closest_r = HUGE(1.0_dp)
         do i=1, k
            r = positions_dist((/x(j)/), (/k_means(i)/),  (/periodicity/))
            if (r < closest_r) then
               closest_r = r
               closest_i = i
            endif
         end do
         a(j) = closest_i
      end do
      iter=iter+1
      if(do_verbose) write(*,*) 'Iteration: ', iter
   end do

   do i=1, k
      closest_i = 0
      closest_r = HUGE(1.0_dp)
      do j=1, n
         r = positions_dist((/x(j)/), (/k_means(i)/),  (/periodicity/))
         if (r < closest_r) then
            closest_r = r
            closest_j = j
         endif
      end do
      cluster_indices(i) = closest_j
   end do

   deallocate(a, prev_a, k_means)

end subroutine k_means_clustering_pick_1d

!--------------------------------------------------------------------------------------------------

subroutine k_means_clustering_pick_nd(x, periodicity, cluster_indices, starting_indices, thresholds, verbose)
   real(dp), intent(in) :: x(:,:)
   real(dp), intent(in) :: periodicity(:)
   integer, intent(out) :: cluster_indices(:)
   integer, optional, intent(in) :: starting_indices(:)
   real(dp), optional, intent(in) :: thresholds(:)
   logical, optional, intent(in) :: verbose

   real(dp), allocatable :: thr(:)
   logical :: do_verbose
   integer :: k, l
   real(dp), allocatable :: k_means(:,:)
   integer :: n, i, j, closest_i, closest_j, nd, d
   real(dp) :: closest_r, r, t_sum(size(x,1))
   integer, allocatable :: a(:), prev_a(:)
   real(dp), allocatable :: half_periodicity(:), shift(:)
   integer :: iter
   logical :: conv

   if(present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   k = size(cluster_indices)
   n = size(x,2)
   nd = size(x,1)
   allocate(k_means(nd,k))
   allocate(half_periodicity(nd), shift(nd))

   ! calculate the half of the periodicity
   half_periodicity = periodicity / 2.0_dp

   l = 0 
   if(present(starting_indices)) then
      l = size(starting_indices)
   end if

   allocate(thr(nd))
   thr = 0.0_dp
   if(present(thresholds)) then
      if(size(thresholds) .ne. nd) then
         print *, "Dimension of convergence=", size(thresholds), " is not equal to nd=", nd
         stop
      end if
      thr = thresholds
   end if
   
   ! random initial guess
   if(l /= 0) then
      call random_selection(n, cluster_indices, starting_indices)
   else
      call random_selection(n, cluster_indices)
   end if

   k_means(:,:) = x(:,cluster_indices(:))

   allocate(a(n), prev_a(n))

   ! evaluate initial assignments
   do j=1, n
      closest_i = 0
      closest_r = HUGE(1.0_dp)
      do i=1, k
         r = positions_dist(x(:,j), k_means(:,i), periodicity(:))
         if (r < closest_r) then
            closest_r = r
            closest_i = i
         endif
      end do
      a(j) = closest_i
   end do

   prev_a = 0
   iter=0
   do while (any(a /= prev_a))
      prev_a = a

      ! update positions
      do i=1, k
         t_sum = 0.0_dp
         do j=1, n
            shift = 0.0_dp
            do d=1, nd
               if (abs(k_means(d,i)-x(d,j)) > half_periodicity(d)) shift(d) = sign(periodicity(d), k_means(d,i))
            end do
            if (a(j) == i) t_sum(:) = t_sum(:) + x(:,j) + shift(:)
         end do
         t_sum(:) = t_sum(:) / real(count(a == i),dp)
         conv = .true.
         do d=1, nd
            if( positions_dist((/k_means(d,i)/), (/t_sum(d)/), (/periodicity(d)/)) .gt. thr(d) ) then
                conv = .false.
                exit
            end if
         end do
         if(.not. conv) then
            k_means(:,i) = t_sum(:)
            do d=1, nd
               if(abs(k_means(d,i)) > half_periodicity(d)) k_means(d,i) = k_means(d,i) - sign(periodicity(d), k_means(d,i))
            end do
         end if
      end do

      ! update assignments
      do j=1, n
         closest_i = 0
         closest_r = HUGE(1.0_dp)
         do i=1, k
            r = positions_dist(x(:,j), k_means(:,i), periodicity(:))
            if (r < closest_r) then
               closest_r = r
               closest_i = i
            endif
         end do
         a(j) = closest_i
      end do
      iter=iter+1
      if(do_verbose) write(*,*) 'Iteration: ', iter
   end do

   do i=1, k
      closest_i = 0
      closest_r = HUGE(1.0_dp)
      do j=1, n
         r = positions_dist(x(:,j), k_means(:,i), periodicity(:))
         if (r < closest_r) then
            closest_r = r
            closest_j = j
         endif
      end do
      cluster_indices(i) = closest_j
   end do

   deallocate(a, prev_a, k_means)
   deallocate(half_periodicity, shift)
   deallocate(thr)

end subroutine k_means_clustering_pick_nd

!--------------------------------------------------------------------------------------------------

subroutine random_selection(n, indices, preselected_indices)
   integer, intent(in) :: n
   integer, intent(out) :: indices(:)
   integer, optional, intent(in) :: preselected_indices(:)

   integer :: i, j, k, l
   integer, allocatable :: ind(:)
   real(dp) :: rv, rrv

   k = size(indices)

   l = 0
   if(present(preselected_indices)) then
      if(any(preselected_indices > n)) then
         print *,"there is at least one index in preselected_indices greater than n=", n
         stop
      end if
      if(any(preselected_indices < 1)) then
         print *,"there is at least one index in preselected_indices less than 1"
         stop
      end if
      l = size(preselected_indices)
   end if

   indices = 0

   if(l > k) then
      if(l/k > 1) then
         allocate(ind(k))
         ind = 0
         do i=1, k
            call random_number(rrv)
            rv = 1+floor(rrv*l)
            do while (any(ind == rv))
               call random_number(rrv)
               rv = 1+floor(rrv*l)
            end do
            ind(i) = rv
         end do
      else
         allocate(ind(l-k))
         ind = 0
         do i=1, l-k
            call random_number(rrv)
            rv = 1+floor(rrv*l)
            do while (any(ind == rv))
               call random_number(rrv)
               rv = 1+floor(rrv*l)
            end do
            ind(i) = rv
         end do
         j = 1
         do i=1, l
            if(all(i /= ind)) then
               indices(j) = preselected_indices(i)
               j = j + 1
               if(j > k) exit
            end if
         end do
      end if
      deallocate(ind)
   else if(l == k) then
      do i=1, k
         indices(i) = preselected_indices(i)
      end do
   else if(k >= n) then
      do i=1, n
         indices(i) = i
      end do
   else
      do i=1, l
         indices(i) = preselected_indices(i)
      end do
      if(n/k > 1) then
         do i=l+1, k
            call random_number(rrv)
            rv = 1+floor(rrv*n)
            do while (any(indices == rv))
               call random_number(rrv)
               rv = 1+floor(rrv*n)
            end do
            indices(i) = rv
         end do
      else
         if((l == 0) .or. (l == n)) then
            allocate(ind(n-k))
            ind = 0
            do i=1, n-k
               call random_number(rrv)
               rv = 1+floor(rrv*n)
               do while (any(ind == rv))
                  call random_number(rrv)
                  rv = 1+floor(rrv*n)
               end do
               ind(i) = rv
            end do
            j = 1
            do i=1, n
               if(all(i /= ind)) then
                  indices(j) = i
                  j = j + 1
                  if(j > k) exit
               end if
            end do
         else
            allocate(ind(n-l))
            j = 1
            do i=1, n
               if(all(i /= preselected_indices)) then
                  ind(j) = i
                  j = j + 1
                  if(j > n-l) exit
               end if
            end do
            do i=l+1, k
               call random_number(rrv)
               rv = 1+floor(rrv*(n-l))
               do while (any(indices(l+1:k) == rv))
                  call random_number(rrv)
                  rv = 1+floor(rrv*(n-l))
               end do
               indices(i) = rv
            end do
            do i=l+1, k
               indices(i) = ind(indices(i))
            end do
         end if
         deallocate(ind)
      end if
   end if
      
end subroutine random_selection

!--------------------------------------------------------------------------------------------------

end module clustering_mod
