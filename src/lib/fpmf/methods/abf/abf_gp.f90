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

!% A simple Gaussian process module for n-D functions learned from samples of
!% their values or gradients, with sparsification

module abf_gp

use abf_dat

implicit none

public :: gp_basic_teach, finalise, f_predict, f_predict_var, f_predict_grad, f_predict_var_grad, f_predict_grad_var
public :: gp_basic_change, gp_basic_merge, gp_basic_complete
public :: gp_basic_write, gp_basic_read
public :: SE_kernel_r_rr

!!% initialise (and teach) a gp_basic
!interface initialise
!   module procedure gp_basic_initialise_nd
!end interface initialise

!% finalise and deallocate a gp_basic
interface finalise
   module procedure gp_basic_finalise
end interface finalise

!% predict a function value from a gp
interface f_predict
   module procedure f_predict_r
end interface f_predict

!% predict a gradient of a function value from a gp
interface f_predict_grad
   module procedure f_predict_grad_r
end interface f_predict_grad

!% predict a function variance from a gp
interface f_predict_var
   module procedure f_predict_var_r
end interface f_predict_var

!% predict the variance of a gradient of a function from a gp
interface f_predict_grad_var
   module procedure f_predict_grad_var_r
end interface f_predict_grad_var

!% predict the gradient of a function variance from a gp
interface f_predict_var_grad
   module procedure f_predict_var_grad_r
end interface f_predict_var_grad

interface Matrix_QR_Solve
   module procedure GP_Matrix_QR_Solve_Vector, GP_Matrix_QR_Solve_Matrix
end interface Matrix_QR_Solve

interface assignment(=)
   module procedure GP_Matrix_Assign
end interface

interface assignment(=)
   module procedure gp_basic_assign
end interface

interface operator(+)
   module procedure gp_basic_add
end interface

contains

subroutine gp_basic_teach(self, len_scale, periodicity, f_var, kernel, f_r, f_v, f_n, g_r, g_v, g_n, &
                          f_pseudo_sparse_set, g_pseudo_sparse_set, jitter, partial, verbose)

   implicit none

   type(gp_basic), intent(inout) :: self !% object to store GP
   real(dp), optional, intent(in) :: len_scale(:)
   real(dp), optional, intent(in) :: periodicity(:)
   real(dp), optional, intent(in) :: f_var !% length scale and variance prior for GP
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), optional, intent(in) :: f_r(:,:), f_v(:), f_n(:) !% arrays of function positions, values, noise 
   real(dp), optional, intent(in) :: g_r(:,:), g_v(:,:), g_n(:,:) !% arrays of function gradient positions, values, noise 
   real(dp), optional, intent(in) :: f_pseudo_sparse_set(:,:), g_pseudo_sparse_set(:,:) !% pseudo positions to use for sparsifcation for values and gradients
   real(dp), optional, intent(in) :: jitter !% jitter for sparse calculation
   logical, optional, intent(in) :: partial
   logical, optional, intent(in) :: verbose

   real(dp), allocatable :: Cmat(:,:), Kmn(:,:), y(:), siginvsq_y(:), Kmn_siginvsq_y(:), siginvsq_Knm(:,:), Kmm(:,:)
   integer :: i, n_f, n_g, n_tot, i_glob, ii
   logical :: do_partial, do_verbose

   real(dp), external :: ddot

   if (present(partial)) then
      do_partial = partial
   else
      do_partial = .false.
   end if

   if (present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   ! check input consistency
   i = count( (/ present(f_r), present(f_v), present(f_n) /) )
   if (i /= 0 .and. i /= 3) then
      print *,"got something other than all or none of f_r, f_v, f_n"
      stop
   endif
   i = count( (/ present(g_r), present(g_v), present(g_n) /) )
   if (i /= 0 .and. i /= 3) then
      print *,"got something other than all or none of g_r, g_v, g_n"
      stop
   endif

   ! compare f_r to f_[gn]
   if (present(f_r)) then
      if (size(f_r,2) /= size(f_v) .or. size(f_r,2) /= size(f_n)) then
         print *,"size(f_r,2)=",size(f_r,2)," size(f_v)=",size(f_v)," size(f_n)=",size(f_n)," not all equal"
         stop
      endif
   endif
   ! compare g_r to g_[gn]
   if (present(g_r)) then
      if (any(shape(g_r) /= shape(g_v)) .or. any(shape(g_r) /= shape(g_n))) then
         print *,"shape(g_r)=",shape(g_r)," shape(g_v)=",shape(g_v)," shape(g_n)=",shape(g_n)," not all equal"
         stop
      endif
   endif
   ! compare f_r to g_r
   if (present(f_r) .and. present(g_r)) then
      if (size(f_r,1) /= size(g_r,1)) then
         print *,"size(f_r,1)=",size(f_r,1)," size(g_r,1)=",size(g_r,1)," not equal"
         stop
      endif
   endif

   if (present(f_r)) then
      n_f = size(f_r,2)
   else
      n_f = 0
   end if

   if (present(g_r)) then
      n_g = size(g_r,2)
   else
      n_g = 0
   end if

   if ( (n_f == 0) .and. (n_g == 0)) then
      print *,"no teaching data provided"
      stop
   endif

   if ( (.not. self%initialised) .or. (.not. self%sparsified) ) then
      call finalise(self)

      if (present(f_r)) then
         self%n_dof = size(f_r,1)
      end if
      if (present(g_r)) then
         self%n_dof = size(g_r,1)
      end if

      if (present(f_pseudo_sparse_set)) then
         if (size(f_pseudo_sparse_set,1) /= self%n_dof) then
            print *,"size(f_pseudo_sparse_set,1)=",size(f_pseudo_sparse_set,1)," n_dof=",self%n_dof," not equal"
            stop
         end if
         self%sparsified = .true.
         self%m_f = size(f_pseudo_sparse_set,2)
         allocate(self%f_r(self%n_dof,self%m_f))
         self%f_r(:,:) = f_pseudo_sparse_set(:,:)
      endif
      if (present(g_pseudo_sparse_set)) then
         if (size(g_pseudo_sparse_set,1) /= self%n_dof) then
            print *,"size(g_pseudo_sparse_set,1)=",size(g_pseudo_sparse_set,1)," n_dof=",self%n_dof," not equal"
            stop
         end if
         self%sparsified = .true.
         self%m_g = size(g_pseudo_sparse_set,2)
         allocate(self%g_r(self%n_dof,self%m_g))
         self%g_r(:,:) = g_pseudo_sparse_set(:,:)
      endif

      if ( (self%m_f == 0) .and. (self%m_g == 0) ) then
         if (present(f_r)) then
            self%m_f = n_f
            allocate(self%f_r(self%n_dof,self%m_f))
            self%f_r(:,:) = f_r(:,:)
         end if
         if (present(g_r)) then
            self%m_g = n_g
            allocate(self%g_r(self%n_dof,self%m_g))
            self%g_r(:,:) = g_r(:,:)
         end if
      end if

      if (present(jitter)) self%jitter = jitter

      if (size(len_scale) /= self%n_dof) then
         print *,"size(len_scale)=",size(len_scale)," , n_dof=",self%n_dof
         stop
      endif

      if (size(periodicity) /= self%n_dof) then
         print *,"size(periodicity)=",size(periodicity)," /= n_dof=",self%n_dof
         stop
      endif

      if (f_var <= 0.0_dp) then
         print *,"invalid f_var=",f_var
         stop
      endif

      do i=1,self%n_dof
         if (periodicity(i) < 0.0_dp) then
            print *,"invalid periodicity(", i, ")=",periodicity(i)
            stop
         end if
      end do

      do i=1,self%n_dof
         if (len_scale(i) <= 0.0_dp) then
            print *,"invalid len_scale(", i, ")=",len_scale(i)
            stop
         end if
      end do

      allocate(self%len_scale_sq(self%n_dof))
      allocate(self%periodicity(self%n_dof))
      self%len_scale_sq = len_scale**2
      self%periodicity = periodicity
      self%f_var = f_var

      if (self%m_f + self%m_g == 0) then
         print *,"no sparsified teaching data provided"
         stop
      endif

      self%m_teach = self%m_f + self%m_g*self%n_dof
   end if

   n_tot = n_f + n_g*self%n_dof

   ! everyone needs Kmm - calculate it only if required and not already available
   if ( (.not. self%sparsified) .or. ((.not. do_partial) .and. (.not. self%Kmm%initialised)) ) then
      allocate(Kmm(self%m_teach, self%m_teach))
      if(do_verbose) print *,"Calculation of Kmm matrix has started..."
      call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, &
                      self%m_f, self%f_r, self%m_g, self%g_r, &
                      self%f_var, self%len_scale_sq, self%periodicity, kernel, Kmm)
      if(do_verbose) print *,"Calculation of Kmm matrix has finished..."
   end if

   if (self%sparsified) then
      allocate(Kmn(self%m_teach, n_tot))
      allocate(Kmn_siginvsq_y(self%m_teach))
      allocate(siginvsq_Knm(n_tot, self%m_teach))
      allocate(siginvsq_y(n_tot))
      allocate(Cmat(self%m_teach, self%m_teach))

      ! calculatn Kmn
      if(do_verbose) print *,"Calculation of Kmn matrix has started..."
      call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, & 
                      n_f, f_r, n_g, g_r, &
                      self%f_var, self%len_scale_sq, self%periodicity, kernel, Kmn)
      if(do_verbose) print *,"Calculation of Kmn matrix has finished..."

      if ( (.not. do_partial) .and. (.not. self%Kmm%initialised) ) then
         ! we'll need Kmm^{-1} for variance
         ! "jitter" (ABP e-mail 14 Aug)
         do i=1, size(Kmm,1)
            Kmm(i,i) = Kmm(i,i) + self%jitter
         end do
         call gp_matrix_initialise(self%Kmm, Kmm)
      end if

      if(do_verbose) print *,"Calculation of siginvsq_y has started..."
      do i=1, n_f
         siginvsq_y(i) = f_v(i)/f_n(i)
      end do
      do i=1, n_g
         siginvsq_y(n_f+(i-1)*self%n_dof+1:n_f+i*self%n_dof) = g_v(:,i)/g_n(:,i)
      end do
      if(do_verbose) print *,"Calculation of siginvsq_y has finished..."

      if(do_verbose) print *,"Calculation of Kmn_siginvsq_y has started..."
      call dgemv('N', size(Kmn,1), size(Kmn,2), 1.0_dp, Kmn(1,1), size(Kmn,1), &
         siginvsq_y(1), 1, 0.0_dp, Kmn_siginvsq_y(1), 1)
      if(do_verbose) print *,"Calculation of Kmn_siginvsq_y has finished..."

      if(do_verbose) print *,"Calculation of siginvsq_Knm has started..."
      siginvsq_Knm = transpose(Kmn)
      do i=1, n_f
         siginvsq_Knm(i,:) = siginvsq_Knm(i,:) / f_n(i)
      end do
      do i=1, n_g
         do ii=1, self%n_dof
            siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) = siginvsq_Knm(n_f+(i-1)*self%n_dof+ii,:) / g_n(ii,i)
         end do
      end do
      if(do_verbose) print *,"Calculation of siginvsq_Knm has finished..."

      Cmat = 0.0_dp
      if(do_verbose) print *,"Calculation of Kmn_siginvsq_Knm has started..."
      call dgemm('N', 'N', size(Cmat,1), size(Cmat,2), size(Kmn,2), &
         1.0_dp, Kmn(1,1), size(Kmn,1), siginvsq_Knm(1,1), size(siginvsq_Knm,1), &
         0.0_dp, Cmat(1,1), size(Cmat,1))
      ! Cmat is now K_{mn} \Sigma^{-2} K_{nm}
      if(do_verbose) print *,"Calculation of Kmn_siginvsq_Knm has finished..."

      if (allocated(self%Kmn_siginvsq_y)) then
         self%Kmn_siginvsq_y = self%Kmn_siginvsq_y + Kmn_siginvsq_y
      else
         allocate(self%Kmn_siginvsq_y(self%m_teach))
         self%Kmn_siginvsq_y = Kmn_siginvsq_y
      end if

      if (allocated(self%Kmn_siginvsq_Knm)) then
         self%Kmn_siginvsq_Knm = self%Kmn_siginvsq_Knm + Cmat
         Cmat = self%Kmn_siginvsq_Knm
      else
         allocate(self%Kmn_siginvsq_Knm(self%m_teach, self%m_teach))
         self%Kmn_siginvsq_Knm = Cmat
      end if

      if (do_partial) then
         self%partial = .true.
         if (self%Cmat%initialised) call gp_matrix_finalise(self%Cmat)
      else
         Cmat = self%Kmm%matrix + Cmat
         ! Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})
         ! Cmat contains jitter via Kmm

         call gp_matrix_initialise(self%Cmat, Cmat)

         if(do_verbose) print *,"Calculation of Cmat_inv_v has started..."
         if(allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
         allocate(self%Cmat_inv_v(self%m_teach))
         call Matrix_QR_Solve(self%Cmat, self%Kmn_siginvsq_y, self%Cmat_inv_v)
         if(do_verbose) print *,"Calculation of Cmat_inv_v has finished..."
      end if

      deallocate(Cmat)
      deallocate(Kmn)
      deallocate(Kmn_siginvsq_y)
      deallocate(siginvsq_Knm)
      deallocate(siginvsq_y)
   else 
      ! not sparsified
      allocate(y(n_tot))

      ! we'll need self%noise_Kmm, which is Kmm shifted by noise, for variance
      if(do_verbose) print *,"Calculation of Cmat has started..."
      do i=1, self%m_f
         Kmm(i,i) = Kmm(i,i) + f_n(i)
      end do
      do i=1, self%m_g
         do ii=1, self%n_dof
            i_glob = self%m_f + (i-1)*self%n_dof+ii
            Kmm(i_glob, i_glob)= Kmm(i_glob, i_glob) + g_n(ii,i)
         end do
      end do
      call gp_matrix_initialise(self%noise_Kmm, Kmm)
      if(do_verbose) print *,"Calculation of Cmat has finished..."

      if (n_f > 0) y(1:n_f) = f_v(1:n_f)
      do i=1, n_g
         i_glob = n_f + (i-1)*self%n_dof
         y(i_glob+1:i_glob+self%n_dof) = g_v(1:self%n_dof,i)
      end do

      if(do_verbose) print *,"Calculation of Cmat_inv_v has started..."
      allocate(self%Cmat_inv_v(self%m_teach))
      call Matrix_QR_Solve(self%noise_Kmm, y, self%Cmat_inv_v)
      if(do_verbose) print *,"Calculation of Cmat_inv_v has finished..."

      deallocate(y)
   endif

   if(allocated(Kmm)) deallocate(Kmm)
   if(.not. allocated(self%k)) allocate(self%k(self%m_teach))
   if(.not. allocated(self%k_grad)) allocate(self%k_grad(self%n_dof, self%m_teach))
   if(.not. allocated(self%mat_inv_k)) allocate(self%mat_inv_k(self%m_teach))
   if(.not. allocated(self%mat_inv_k_grad)) allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))

   self%initialised = .true.
end subroutine gp_basic_teach

subroutine gp_basic_finalise(self)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP

   self%initialised = .false.
   self%sparsified = .false.
   self%partial = .false.

   self%n_dof = 0
   self%m_f = 0
   self%m_g = 0
   self%m_teach = 0
   self%jitter = 0.0_dp
   self%f_var = 0.0_dp

   if (allocated(self%len_scale_sq)) deallocate(self%len_scale_sq)
   if (allocated(self%periodicity)) deallocate(self%periodicity)
   if (allocated(self%f_r)) deallocate(self%f_r)
   if (allocated(self%g_r)) deallocate(self%g_r)
   if (allocated(self%k)) deallocate(self%k)
   if (allocated(self%k_grad)) deallocate(self%k_grad)
   if (allocated(self%mat_inv_k)) deallocate(self%mat_inv_k)
   if (allocated(self%mat_inv_k_grad)) deallocate(self%mat_inv_k_grad)
   if (allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
   if (allocated(self%Kmn_siginvsq_y)) deallocate(self%Kmn_siginvsq_y)
   if (allocated(self%Kmn_siginvsq_Knm)) deallocate(self%Kmn_siginvsq_Knm)
   call gp_matrix_finalise(self%Kmm)
   call gp_matrix_finalise(self%Cmat)
   call gp_matrix_finalise(self%noise_Kmm)

end subroutine gp_basic_finalise

subroutine gp_basic_assign(left_self, right_self)

   implicit none

   type(gp_basic), intent(inout) :: left_self
   type(gp_basic), intent(in) :: right_self

   call gp_basic_finalise(left_self)

   left_self%initialised = right_self%initialised
   if (.not. left_self%initialised) return

   left_self%sparsified = right_self%sparsified
   left_self%partial = right_self%partial
       
   left_self%n_dof = right_self%n_dof
   left_self%m_f = right_self%m_f
   left_self%m_g = right_self%m_g
   left_self%m_teach = right_self%m_teach
   left_self%jitter = right_self%jitter
   left_self%f_var = right_self%f_var

   if (allocated(right_self%len_scale_sq)) then
      allocate(left_self%len_scale_sq(left_self%n_dof))
      left_self%len_scale_sq = right_self%len_scale_sq
   end if
   if (allocated(right_self%periodicity)) then
      allocate(left_self%periodicity(left_self%n_dof))
      left_self%periodicity = right_self%periodicity
   end if
   if (allocated(right_self%f_r)) then
      allocate(left_self%f_r(left_self%n_dof,left_self%m_f))
      left_self%f_r = right_self%f_r
   end if
   if (allocated(right_self%g_r)) then
      allocate(left_self%g_r(left_self%n_dof,left_self%m_g))
      left_self%g_r = right_self%g_r
   end if
   if (allocated(right_self%k)) then
      allocate(left_self%k(left_self%m_teach))
      left_self%k = right_self%k
   end if
   if (allocated(right_self%k_grad)) then
      allocate(left_self%k_grad(left_self%n_dof,left_self%m_teach))
      left_self%k_grad = right_self%k_grad
   end if
   if (allocated(right_self%mat_inv_k)) then
      allocate(left_self%mat_inv_k(left_self%m_teach))
      left_self%mat_inv_k = right_self%mat_inv_k
   end if
   if (allocated(right_self%mat_inv_k_grad)) then
      allocate(left_self%mat_inv_k_grad(left_self%m_teach,left_self%n_dof))
      left_self%mat_inv_k_grad = right_self%mat_inv_k_grad
   end if
   if (allocated(right_self%Cmat_inv_v)) then
      allocate(left_self%Cmat_inv_v(left_self%m_teach))
      left_self%Cmat_inv_v = right_self%Cmat_inv_v
   end if
   if (allocated(right_self%Kmn_siginvsq_y)) then
      allocate(left_self%Kmn_siginvsq_y(left_self%m_teach))
      left_self%Kmn_siginvsq_y = right_self%Kmn_siginvsq_y
   end if
   if (allocated(right_self%Kmn_siginvsq_Knm)) then
      allocate(left_self%Kmn_siginvsq_Knm(left_self%m_teach,left_self%m_teach))
      left_self%Kmn_siginvsq_Knm = right_self%Kmn_siginvsq_Knm
   end if

   left_self%Kmm = right_self%Kmm
   left_self%noise_Kmm = right_self%noise_Kmm
   left_self%Cmat = right_self%Cmat

end subroutine gp_basic_assign

function gp_basic_add(self1, self2) result(self3)

   implicit none

   type(gp_basic), intent(in) :: self1
   type(gp_basic), intent(in) :: self2
   type(gp_basic) :: self3

   integer :: i, j

   if ( .not. self1%initialised ) then
      print *,"self1 is not initialized!"
      stop
   end if

   if ( .not. self2%initialised ) then
      print *,"self2 is not initialized!"
      stop
   end if

   ! check compatibility of self1 and self2
   if (self1%n_dof .ne. self2%n_dof) then
      print *,"self1%n_dof=", self1%n_dof, " is not equal to self2%n_dof=", self2%n_dof
      stop
   end if

   if (self1%m_f .ne. self2%m_f) then
      print *,"self1%m_f=", self1%m_f, " is not equal to self2%m_f=", self2%m_f
      stop
   end if

   if (self1%m_g .ne. self2%m_g) then
      print *,"self1%m_g=", self1%m_g, " is not equal to self2%m_g=", self2%m_g
      stop
   end if

   if (self1%m_teach .ne. self2%m_teach) then
      print *,"self1%m_teach=", self1%m_teach, " is not equal to self2%m_teach=", self2%m_teach
      stop
   end if

   if (self1%jitter .ne. self2%jitter) then
      print *,"self1%jitter=", self1%jitter, " is not equal to self2%jitter=", self2%jitter
      stop
   end if

   if (self1%f_var .ne. self2%f_var) then
      print *,"self1%f_var=", self1%f_var, " is not equal to self2%f_var=", self2%f_var
      stop
   end if

   if (.not. self1%sparsified) then
      print *,"self1 is not sparsified!"
      stop
   end if

   if (.not. self2%sparsified) then
      print *,"self2 is not sparsified!"
      stop
   end if

   do i=1,self1%n_dof
      if (self1%len_scale_sq(i) .ne. self2%len_scale_sq(i)) then
         print *,"self1%len_scale_sq(",i,")=",self1%len_scale_sq(i)," is not equal to self2%len_scale_sq(",i,")=",self2%len_scale_sq(i)
         stop
      end if
   end do

   do i=1,self1%n_dof
      if (self1%periodicity(i) .ne. self2%periodicity(i)) then
         print *,"self1%periodicity(",i,")=",self1%periodicity(i)," is not equal to self2%periodicity(",i,")=",self2%periodicity(i)
         stop
      end if
   end do

   if (self1%m_f .ne. 0) then
      do i=1,self1%n_dof
         do j=1,self1%m_f
            if (self1%f_r(i,j) .ne. self2%f_r(i,j)) then
               print *,"self1%f_r(",i,j,")=",self1%f_r(i,j)," is not equal to self2%f_r(",i,j,")=",self2%f_r(i,j)
               stop
            end if
         end do
      end do
   end if

   if (self1%m_g .ne. 0) then
      do i=1,self1%n_dof
         do j=1,self1%m_g
            if (self1%g_r(i,j) .ne. self2%g_r(i,j)) then
               print *,"self1%g_r(",i,j,")=",self1%g_r(i,j)," is not equal to self2%g_r(",i,j,")=",self2%g_r(i,j)
               stop
            end if
         end do
      end do
   end if

   call finalise(self3)

   self3%sparsified = self1%sparsified
   self3%partial = .true.

   self3%n_dof = self1%n_dof
   self3%m_f = self1%m_f
   self3%m_g = self1%m_g
   self3%m_teach = self1%m_teach
   self3%jitter = self1%jitter
   self3%f_var = self1%f_var

   allocate(self3%len_scale_sq(self3%n_dof))
   self3%len_scale_sq = self1%len_scale_sq
   allocate(self3%periodicity(self3%n_dof))
   self3%periodicity = self1%periodicity
   if (self3%m_f .ne. 0) then
      allocate(self3%f_r(self3%n_dof,self3%m_f))
      self3%f_r = self1%f_r
   end if
   if (self3%m_g .ne. 0) then
      allocate(self3%g_r(self3%n_dof,self3%m_g))
      self3%g_r = self1%g_r
   end if
   allocate(self3%Kmn_siginvsq_y(self3%m_teach))
   self3%Kmn_siginvsq_y = self1%Kmn_siginvsq_y + self2%Kmn_siginvsq_y
   allocate(self3%Kmn_siginvsq_Knm(self3%m_teach,self3%m_teach))
   self3%Kmn_siginvsq_Knm = self1%Kmn_siginvsq_Knm + self2%Kmn_siginvsq_Knm

end function gp_basic_add

subroutine gp_basic_write(self,filename)

   implicit none

   type(gp_basic), intent(in) :: self !% object for GP
   character(len=*), intent(in) :: filename !% filename to write GP out

   if (.not. self%initialised) then
      return
   endif

   open(8,file=filename)

   write(8,*) '# self%n_dof'
   write(8,*) self%n_dof
   write(8,*) '# self%m_f'
   write(8,*) self%m_f
   write(8,*) '# self%m_g'
   write(8,*) self%m_g
   write(8,*) '# self%m_teach'
   write(8,*) self%m_teach
   write(8,*) '# self%f_var'
   write(8,*) self%f_var
   write(8,*) '# self%sparsified'
   write(8,*) self%sparsified
   write(8,*) '# self%len_scale_sq(self%n_dof)'
   write(8,*) self%len_scale_sq
   write(8,*) '# self%periodicity(self%n_dof)'
   write(8,*) self%periodicity
   if(self%m_f .ne. 0) then
      write(8,*) '# self%f_r(self%n_dof,self%m_f)'
      write(8,*) self%f_r
   end if
   if(self%m_g .ne. 0) then
      write(8,*) '# self%g_r(self%n_dof,self%m_g)'
      write(8,*) self%g_r
   end if
   if(self%sparsified) then
      write(8,*) '# self%jitter'
      write(8,*) self%jitter
      write(8,*) '# self%partial'
      write(8,*) self%partial
      write(8,*) '# self%Kmn_siginvsq_y(self%m_teach)'
      write(8,*) self%Kmn_siginvsq_y
      write(8,*) '# self%Kmn_siginvsq_Knm(self%m_teach,self%m_teach)'
      write(8,*) self%Kmn_siginvsq_Knm
      if(.not. self%partial) then
         write(8,*) '# self%Cmat_inv_v(self%m_teach)'
         write(8,*) self%Cmat_inv_v
         write(8,*) '# self%Kmm%n'
         write(8,*) self%Kmm%n
         write(8,*) '# self%Kmm%m'
         write(8,*) self%Kmm%m
         write(8,*) '# self%Kmm%matrix(self%Kmm%n,self%Kmm%m)'
         write(8,*) self%Kmm%matrix
         write(8,*) '# self%Cmat%n'
         write(8,*) self%Cmat%n
         write(8,*) '# self%Cmat%m'
         write(8,*) self%Cmat%m
         write(8,*) '# self%Cmat%matrix(self%Cmat%n,self%Cmat%m)'
         write(8,*) self%Cmat%matrix
      end if
   else
      write(8,*) '# self%Cmat_inv_v(self%m_teach)'
      write(8,*) self%Cmat_inv_v
!     write(8,*) '# self%noise_Kmm%n'
!     write(8,*) self%noise_Kmm%n
!     write(8,*) '# self%noise_Kmm%m'
!     write(8,*) self%noise_Kmm%m
!     write(8,*) '# self%noise_Kmm%matrix(self%noise_Kmm%n,self%noise_Kmm%m)'
!     write(8,*) self%noise_Kmm%matrix
   end if

   close(8)

end subroutine gp_basic_write

subroutine gp_basic_read(self, filename)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   character(len=*), intent(in) :: filename !% filename to read GP from

   integer :: n, m
   real(dp), allocatable :: matrix(:,:)

   open(5,file=filename)

   call finalise(self)

   self%initialised=.true.

   read(5,*)
   read(5,*) self%n_dof
   read(5,*)
   read(5,*) self%m_f 
   read(5,*)
   read(5,*) self%m_g
   read(5,*)
   read(5,*) self%m_teach
   read(5,*)
   read(5,*) self%f_var
   read(5,*)
   read(5,*) self%sparsified
   read(5,*)
   allocate(self%len_scale_sq(self%n_dof))
   read(5,*) self%len_scale_sq
   read(5,*)
   allocate(self%periodicity(self%n_dof))
   read(5,*) self%periodicity
   if(self%m_f .ne. 0) then
      allocate(self%f_r(self%n_dof,self%m_f))
      read(5,*)
      read(5,*) self%f_r
   endif
   if(self%m_g .ne. 0) then
      allocate(self%g_r(self%n_dof,self%m_g))
      read(5,*)
      read(5,*) self%g_r
   endif
   if(self%sparsified) then
      read(5,*)
      read(5,*) self%jitter
      read(5,*)
      read(5,*) self%partial
         allocate(self%Kmn_siginvsq_y(self%m_teach))
         read(5,*)
         read(5,*) self%Kmn_siginvsq_y
         allocate(self%Kmn_siginvsq_Knm(self%m_teach,self%m_teach))
         read(5,*)
         read(5,*) self%Kmn_siginvsq_Knm
      if(.not. self%partial) then
         allocate(self%Cmat_inv_v(self%m_teach))
         read(5,*) 
         read(5,*) self%Cmat_inv_v
         allocate(self%k(self%m_teach))
         allocate(self%k_grad(self%n_dof, self%m_teach))
         allocate(self%mat_inv_k(self%m_teach))
         allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))
         read(5,*)
         read(5,*) n
         read(5,*)
         read(5,*) m
         allocate(matrix(n,m)) 
         read(5,*)
         read(5,*) matrix
         call gp_matrix_initialise(self%Kmm, matrix)
         deallocate(matrix)
         read(5,*)
         read(5,*) n
         read(5,*)
         read(5,*) m
         allocate(matrix(n,m))
         read(5,*)
         read(5,*) matrix
         call gp_matrix_initialise(self%Cmat, matrix)
         deallocate(matrix)
      end if
   else
      allocate(self%Cmat_inv_v(self%m_teach))
      read(5,*)
      read(5,*) self%Cmat_inv_v
!     read(5,*)
!     read(5,*) n
!     read(5,*) 
!     read(5,*) m
!     allocate(matrix(n,m))
!     read(5,*)
!     read(5,*) matrix
!     call gp_matrix_initialise(self%noise_Kmm, matrix)
!     deallocate(matrix)
      allocate(self%k(self%m_teach))
      allocate(self%k_grad(self%n_dof, self%m_teach))
      allocate(self%mat_inv_k(self%m_teach))
      allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))
   end if

   close(5)

end subroutine gp_basic_read

! this subroutine changes the position of sparse points
subroutine gp_basic_change(self, kernel, f_set, g_set, verbose)

   implicit none

   type(gp_basic), intent(inout) :: self
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), optional, intent(in) :: f_set(:,:)
   real(dp), optional, intent(in) :: g_set(:,:)
   logical, optional, intent(in) :: verbose

   integer :: i
   type(gp_basic) :: newself
   logical :: do_verbose
   real(dp), allocatable :: Kll(:,:), Klm(:,:), Qmm(:,:)
   real(dp), allocatable :: Kmm_inv_Kml(:,:), Qmm_Kmm_inv_Kml(:,:)

   real(dp), external :: ddot

   if ( .not. self%initialised ) then
      return
   end if

   if ( .not. self%sparsified ) then
      print *,"self is not sparsified!"
      stop
   end if

   if (present(f_set)) then
      if (size(f_set,1) .ne. self%n_dof) then
         print *,"size(f_set,1)=",size(f_set,1)," is not equal to self%n_dof=",self%n_dof
      end if
   end if

   if (present(g_set)) then
      if (size(g_set,1) .ne. self%n_dof) then
         print *,"size(g_set,1)=",size(g_set,1)," is not equal to self%n_dof=",self%n_dof
      end if
   end if

   if (present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   call finalise(newself)

   newself%initialised = .true.
   newself%sparsified = .true.
   newself%partial = .true.

   newself%n_dof = self%n_dof
   if(present(f_set)) newself%m_f = size(f_set,2)
   if(present(g_set)) newself%m_g = size(g_set,2)
   if (newself%m_f + newself%m_g == 0) then
      print *,"no new set of sparse points provided"
      stop
   end if
   newself%m_teach = newself%m_f + newself%m_g*newself%n_dof

   newself%jitter = self%jitter
   newself%f_var = self%f_var

   allocate(newself%len_scale_sq(newself%n_dof))
   newself%len_scale_sq = self%len_scale_sq
   allocate(newself%periodicity(newself%n_dof))
   newself%periodicity = self%periodicity

   if (newself%m_f .ne. 0) then
      allocate(newself%f_r(newself%n_dof,newself%m_f))
      newself%f_r = f_set
   end if
   if (newself%m_g .ne. 0) then
      allocate(newself%g_r(newself%n_dof,newself%m_g))
      newself%g_r = g_set
   end if

   ! complete self if necessary
   if(do_verbose) print *,"Completion of original teching has started..."
   call gp_basic_complete(self, kernel, do_verbose)
   if(do_verbose) print *,"Completion of original teching has finished..."

   ! calculate Kll
   allocate(Kll(newself%m_teach, newself%m_teach))
   if(do_verbose) print *,"Calculation of Kll matrix has started..."
   call kernel_mat(newself%m_f, newself%f_r, newself%m_g, newself%g_r, &
                   newself%m_f, newself%f_r, newself%m_g, newself%g_r, &
                   newself%f_var, newself%len_scale_sq, newself%periodicity, kernel, Kll)
   if(do_verbose) print *,"Calculation of Kll matrix has finished..."
   ! "jitter" (ABP e-mail 14 Aug)
   do i=1, size(Kll,1)
      Kll(i,i) = Kll(i,i) + newself%jitter
   end do
   call gp_matrix_initialise(newself%Kmm, Kll)

   ! calculate Klm
   allocate(Klm(newself%m_teach, self%m_teach))
   if(do_verbose) print *,"Calculation of Klm matrix has started..."
   call kernel_mat(newself%m_f, newself%f_r, newself%m_g, newself%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, &
                   newself%f_var, newself%len_scale_sq, newself%periodicity, kernel, Klm)
   if(do_verbose) print *,"Calculation of Klm matrix has finished..."

   ! calculate Kmm_inv_Kml
   allocate(Kmm_inv_Kml(self%m_teach, newself%m_teach))
   if(do_verbose) print *,"Calculation of Kmm_inv_Kml matrix has started..."
   call Matrix_QR_Solve(self%Kmm, transpose(Klm), Kmm_inv_Kml)
   if(do_verbose) print *,"Calculation of Kmm_inv_Kml matrix has finished..."

   ! calculation of new Kmn_siginvsq_y
   allocate(newself%Kmn_siginvsq_y(newself%m_teach))
   newself%Kmn_siginvsq_y = 0.0_dp
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_y has started..."
   call dgemv('T', size(Kmm_inv_Kml,1), size(Kmm_inv_Kml,2), 1.0_dp, Kmm_inv_Kml(1,1), size(Kmm_inv_Kml,1), &
              self%Kmn_siginvsq_y(1), 1, 0.0_dp, newself%Kmn_siginvsq_y(1), 1)
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_y has finished..."

   ! calculate Qmm
   allocate(Qmm(self%m_teach,self%m_teach))
!  Qmm = self%Kmm%matrix + self%Kmn_siginvsq_Knm LAM !!! why?
   Qmm = self%Kmn_siginvsq_Knm

   ! calculate Qmm_Kmm_inv_Kml
   allocate(Qmm_Kmm_inv_Kml(self%m_teach, newself%m_teach))
   Qmm_Kmm_inv_Kml = 0.0_dp
   if(do_verbose) print *,"Calculation of Qmm_Kmm_inv_Kml matrix has started..."
   call dgemm('N', 'N', size(Qmm,1), size(Kmm_inv_Kml,2), size(Qmm,2), 1.0_dp, &
              Qmm(1,1), size(Qmm,1), Kmm_inv_Kml(1,1), size(Qmm,2), 0.0_dp, &
              Qmm_Kmm_inv_Kml(1,1), size(Qmm,1))
   if(do_verbose) print *,"Calculation of Qmm_Kmm_inv_Kml matrix has finished..."

   ! calculation of new Kmn_siginvsq_Knm
   allocate(newself%Kmn_siginvsq_Knm(newself%m_teach, newself%m_teach))
   newself%Kmn_siginvsq_Knm = 0.0_dp
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_Knm has started..."
   call dgemm('T', 'N', size(Kmm_inv_Kml,2), size(Qmm_Kmm_inv_Kml,2), size(Kmm_inv_Kml,1), 1.0_dp, &
              Kmm_inv_Kml(1,1), size(Kmm_inv_Kml,1), Qmm_Kmm_inv_Kml(1,1), size(Kmm_inv_Kml,1), 0.0_dp, &
              newself%Kmn_siginvsq_Knm(1,1), size(Kmm_inv_Kml,2))
   if(do_verbose) print *,"Calculation of new Kmn_siginvsq_Knm has finished..."

   deallocate(Kll)
   deallocate(Klm)
   deallocate(Kmm_inv_Kml)
   deallocate(Qmm)
   deallocate(Qmm_Kmm_inv_Kml)

   self = newself

   call finalise(newself)

end subroutine gp_basic_change

! this subroutine merges two teachings that have the same teaching points
! if self3 is specified then the result is written into self3, otherwise into self2 
subroutine gp_basic_merge(self1, self2, self3)

   implicit none

   type(gp_basic), intent(inout) :: self1, self2
   type(gp_basic), optional, intent(inout) :: self3

   integer :: i, j

   if ( .not. self1%initialised ) then
      print *,"self1 is not initialized!"
      stop
   end if

   if ( .not. self2%initialised ) then
      print *,"self2 is not initialized!"
      stop
   end if

   ! check compatibility of self1 and self2
   if (self1%n_dof .ne. self2%n_dof) then
      print *,"self1%n_dof=", self1%n_dof, " is not equal to self2%n_dof=", self2%n_dof
      stop
   end if

   if (self1%m_f .ne. self2%m_f) then
      print *,"self1%m_f=", self1%m_f, " is not equal to self2%m_f=", self2%m_f
      stop
   end if

   if (self1%m_g .ne. self2%m_g) then
      print *,"self1%m_g=", self1%m_g, " is not equal to self2%m_g=", self2%m_g
      stop
   end if

   if (self1%m_teach .ne. self2%m_teach) then
      print *,"self1%m_teach=", self1%m_teach, " is not equal to self2%m_teach=", self2%m_teach
      stop
   end if

   if (self1%f_var .ne. self2%f_var) then
      print *,"self1%f_var=", self1%f_var, " is not equal to self2%f_var=", self2%f_var
      stop
   end if

   if (.not. self1%sparsified) then
      print *,"self1 is not sparsified!"
      stop
   end if

   if (.not. self2%sparsified) then
      print *,"self2 is not sparsified!"
      stop
   end if

   do i=1,self1%n_dof
      if (self1%len_scale_sq(i) .ne. self2%len_scale_sq(i)) then
         print *,"self1%len_scale_sq(",i,")=",self1%len_scale_sq(i)," is not equal to self2%len_scale_sq(",i,")=",self2%len_scale_sq(i)
         stop  
      end if
   end do

   do i=1,self1%n_dof
      if (self1%periodicity(i) .ne. self2%periodicity(i)) then
         print *,"self1%periodicity(",i,")=",self1%periodicity(i)," is not equal to self2%periodicity(",i,")=",self2%periodicity(i)
         stop
      end if
   end do

   if (self1%m_f .ne. 0) then
      do i=1,self1%n_dof
         do j=1,self1%m_f
            if (self1%f_r(i,j) .ne. self2%f_r(i,j)) then
               print *,"self1%f_r(",i,j,")=",self1%f_r(i,j)," is not equal to self2%f_r(",i,j,")=",self2%f_r(i,j)
               stop
            end if
         end do
      end do
   end if

   if (self1%m_g .ne. 0) then
      do i=1,self1%n_dof
         do j=1,self1%m_g
            if (self1%g_r(i,j) .ne. self2%g_r(i,j)) then
               print *,"self1%g_r(",i,j,")=",self1%g_r(i,j)," is not equal to self2%g_r(",i,j,")=",self2%g_r(i,j)
               stop
            end if
         end do
      end do
   end if

   if ( present(self3) ) then
      call finalise(self3)

      self3%n_dof = self1%n_dof
      self3%m_f = self1%m_f
      self3%m_g = self1%m_g 
      self3%m_teach = self1%m_teach
      self3%f_var = self1%f_var
      self3%sparsified = self1%sparsified
      allocate(self3%len_scale_sq(self3%n_dof))
      self3%len_scale_sq = self1%len_scale_sq
      allocate(self3%periodicity(self3%n_dof))
      self3%periodicity = self1%periodicity
      if (self3%m_f .ne. 0) then
         allocate(self3%f_r(self3%n_dof,self3%m_f))
         self3%f_r = self1%f_r
      end if
      if (self3%m_g .ne. 0) then
         allocate(self3%g_r(self3%n_dof,self3%m_g))
         self3%g_r = self1%g_r
      end if

      self3%partial = .true.
      allocate(self3%Kmn_siginvsq_y(self3%m_teach))
      self3%Kmn_siginvsq_y = self1%Kmn_siginvsq_y + self2%Kmn_siginvsq_y
      allocate(self3%Kmn_siginvsq_Knm(self3%m_teach,self3%m_teach))
      self3%Kmn_siginvsq_Knm = self1%Kmn_siginvsq_Knm + self2%Kmn_siginvsq_Knm

      self3%initialised = .true.
   else
      self2%partial = .true.
      self2%Kmn_siginvsq_y = self1%Kmn_siginvsq_y + self2%Kmn_siginvsq_y
      self2%Kmn_siginvsq_Knm = self1%Kmn_siginvsq_Knm + self2%Kmn_siginvsq_Knm
   end if
    
end subroutine gp_basic_merge 

! this subroutine completes a teaching (i.e. calculates Kmm and Cmat) if gp_basic is not complete
subroutine gp_basic_complete(self, kernel, verbose)

   implicit none

   type(gp_basic), intent(inout) :: self
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   logical, optional, intent(in) :: verbose

   logical :: do_verbose 
   integer :: i
   real(dp), allocatable :: Kmm(:,:), Cmat(:,:)

   if (.not. self%initialised) then
      print *,"self is not initialised!"
      stop
   end if

   if (.not. self%partial) return

   if (present(verbose)) then
      do_verbose = verbose
   else
      do_verbose = .false.
   end if

   allocate(Kmm(self%m_teach, self%m_teach))
   if(do_verbose) print *,"Calculation of Kmm has started..."
   call kernel_mat(self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%m_f, self%f_r, self%m_g, self%g_r, &
                   self%f_var, self%len_scale_sq, self%periodicity, kernel, Kmm)
   if(do_verbose) print *,"Calculation of Kmm has finished..."
   ! "jitter" (ABP e-mail 14 Aug)
   do i=1, size(Kmm,1)
      Kmm(i,i) = Kmm(i,i) + self%jitter
   end do
   call gp_matrix_initialise(self%Kmm, Kmm)
   
   allocate(Cmat(self%m_teach, self%m_teach))

   Cmat = self%Kmm%matrix + self%Kmn_siginvsq_Knm
   ! Cmat is now (K_{mm} + K_{mn} \Sigma^{-2} K_{nm})
   ! Cmat contains jitter via Kmm

   call gp_matrix_initialise(self%Cmat, Cmat)

   if(do_verbose) print *,"Calculation of Cmat_inv_v has started..."
   if(allocated(self%Cmat_inv_v)) deallocate(self%Cmat_inv_v)
   allocate(self%Cmat_inv_v(self%m_teach))
   call Matrix_QR_Solve(self%Cmat, self%Kmn_siginvsq_y, self%Cmat_inv_v)
   if(do_verbose) print *,"Calculation of Cmat_inv_v has finished..."

   self%partial = .false.

   deallocate(Kmm)
   deallocate(Cmat)

   if(.not. allocated(self%k)) allocate(self%k(self%m_teach))
   if(.not. allocated(self%k_grad)) allocate(self%k_grad(self%n_dof, self%m_teach))
   if(.not. allocated(self%mat_inv_k)) allocate(self%mat_inv_k(self%m_teach))
   if(.not. allocated(self%mat_inv_k_grad)) allocate(self%mat_inv_k_grad(self%m_teach, self%n_dof))

end subroutine gp_basic_complete

function f_predict_r(self, r, kernel)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_r

   real(dp), external :: ddot

   if (.not. self%initialised) then
      f_predict_r = 0.0_dp
      return
   endif

   if (self%partial) then
      call gp_basic_complete(self, kernel)
   end if

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, kernel, self%k)
   f_predict_r = ddot(size(self%k), self%k(1), 1, self%Cmat_inv_v(1), 1)
end function f_predict_r

function f_predict_grad_r(self, r, kernel)

   implicit none

   type(gp_basic), intent(inout) :: self !% object for GP
   real(dp), intent(in) :: r(:) !% position at which to f_predict value
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_r(size(r))

   if (.not. self%initialised) then
      f_predict_grad_r = 0.0_dp
      return
   endif

   if (self%partial) then
      call gp_basic_complete(self, kernel)
   end if

   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, &
      self%f_var, self%len_scale_sq, self%periodicity, kernel, self%k_grad)
   call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), 1.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
      self%Cmat_inv_v(1), 1, 0.0_dp, f_predict_grad_r(1), 1)
end function f_predict_grad_r

function f_predict_var_r(self, r, kernel)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_r

   real(dp), external :: ddot

   if (.not. self%initialised) then
      f_predict_var_r = 0.0_dp
      return
   endif

   if (self%partial) then
      call gp_basic_complete(self, kernel)
   end if

   call f_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, self%f_var, self%len_scale_sq, self%periodicity, kernel, self%k)

   f_predict_var_r = self%f_var
   if (self%sparsified) then
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - ddot(size(self%k), self%k(1), 1, self%mat_inv_k(1), 1)
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r + ddot(size(self%k), self%k(1), 1, self%mat_inv_k(1), 1)
   else
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      f_predict_var_r = f_predict_var_r - ddot(size(self%k), self%k(1), 1, self%mat_inv_k(1), 1)
   endif
end function f_predict_var_r

function f_predict_grad_var_r(self, r, kernel)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_grad_var_r(size(r))

   integer :: ii

   real(dp), external :: ddot

   if (.not. self%initialised) then
      f_predict_grad_var_r = 0.0_dp
      return
   endif

   if (self%partial) then
      call gp_basic_complete(self, kernel)
   end if

   ! From Lockwood and Anitescu preprint ANL/MCS-P1808-1110 "Gradient-Enhanced Universal Kriging for Uncertainty Propagation"
   ! Eq. 2.22, not including last term which is relevant only for underlying polynomial basis (which we don't have)
   call g_kernel_vec(r, self%m_f, self%f_r, self%m_g, self%g_r, &
      self%f_var, self%len_scale_sq, self%periodicity, kernel, self%k_grad)

   ! first term should be second derivative of covariance, but it's hard wired for now
   f_predict_grad_var_r = self%f_var/self%len_scale_sq 
   if (self%sparsified) then
      ! -k' K_{mm}^{-1} k'
      call Matrix_QR_Solve(self%Kmm, transpose(self%k_grad), self%mat_inv_k_grad)
      ! first term should be second derivative of covariance, but it's hard wired for now
      do ii=1, self%n_dof
         f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) - ddot(size(self%k_grad,2), self%k_grad(ii,1), size(self%k_grad,1), &
                                                                                         self%mat_inv_k_grad(1,ii), 1)
      end do
      ! +k' C^{-1} k'
      call Matrix_QR_Solve(self%Cmat, transpose(self%k_grad), self%mat_inv_k_grad)
      do ii=1, self%n_dof
         f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) + ddot(size(self%k_grad,2), self%k_grad(ii,1), size(self%k_grad,1), &
                                                                                         self%mat_inv_k_grad(1,ii), 1)
      end do
   else
      ! -k' C^{-1} k'
      call Matrix_QR_Solve(self%noise_Kmm, transpose(self%k_grad), self%mat_inv_k_grad)
      do ii=1, self%n_dof
         f_predict_grad_var_r(ii) = f_predict_grad_var_r(ii) - ddot(size(self%k_grad,2), self%k_grad(ii,1), size(self%k_grad,1), &
                                                                                         self%mat_inv_k_grad(1,ii), 1)
      end do
   endif
end function f_predict_grad_var_r

function f_predict_var_grad_r(self, r, kernel)

   implicit none

   type(gp_basic), intent(inout) :: self
   real(dp), intent(in) :: r(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp) :: f_predict_var_grad_r(size(r))

   if (.not. self%initialised) then
      f_predict_var_grad_r = 0.0_dp
      return
   endif

   if (self%partial) then
      call gp_basic_complete(self, kernel)
   end if

   call f_kernel_vec(r,self%m_f,self%f_r,self%m_g,self%g_r,self%f_var,self%len_scale_sq,self%periodicity,kernel,self%k)
   call g_kernel_vec(r,self%m_f,self%f_r,self%m_g,self%g_r,self%f_var,self%len_scale_sq,self%periodicity,kernel,self%k_grad)

   f_predict_var_grad_r = 0.0_dp
   if (self%sparsified) then
      call Matrix_QR_Solve(self%Kmm, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), -2.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
         self%mat_inv_k(1), 1, 1.0_dp, f_predict_var_grad_r(1), 1)
      call Matrix_QR_Solve(self%Cmat, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), 2.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
         self%mat_inv_k(1), 1, 1.0_dp, f_predict_var_grad_r(1), 1)
   else
      call Matrix_QR_Solve(self%noise_Kmm, self%k, self%mat_inv_k)
      call dgemv('N', size(self%k_grad,1), size(self%k_grad,2), -2.0_dp, self%k_grad(1,1), size(self%k_grad,1), &
         self%mat_inv_k(1), 1, 1.0_dp, f_predict_var_grad_r(1), 1)
   endif
end function f_predict_var_grad_r

subroutine kernel_mat(l_n_f, l_f_r, l_n_g, l_g_r, r_n_f, r_f_r, r_n_g, r_g_r, f_var, len_scale_sq, periodicity, kernel, mat)

   implicit none

   integer, intent(in) :: l_n_f, l_n_g
   real(dp), optional, intent(in) :: l_f_r(:,:), l_g_r(:,:)
   integer, intent(in) :: r_n_f, r_n_g
   real(dp), optional, intent(in) :: r_f_r(:,:), r_g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: mat(:,:)

   integer :: i, i_glob, n_dof

   if (l_n_f > 0) n_dof=size(l_f_r,1)
   if (l_n_g > 0) n_dof=size(l_g_r,1)

   do i=1, l_n_f
      ! f on f
      if (r_n_f > 0) call kernel(l_f_r(1:n_dof,i), r_f_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                                 f_f=mat(i,1:r_n_f))
      ! f on g
      if (r_n_g > 0) call kernel(l_f_r(1:n_dof,i), r_g_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                                 f_g=mat(i,r_n_f+1:r_n_f+r_n_g*n_dof))
   end do
   do i=1, l_n_g
      i_glob = l_n_f + (i-1)*n_dof + 1
      ! g on f
      if (r_n_f > 0) call kernel(l_g_r(1:n_dof,i), r_f_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                                 g_f=mat(i_glob:i_glob+n_dof-1,1:r_n_f))
      ! g on g
      if (r_n_g > 0) call kernel(l_g_r(1:n_dof,i), r_g_r(1:n_dof,:), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                                 g_g=mat(i_glob:i_glob+n_dof-1,r_n_f+1:r_n_f+r_n_g*n_dof))
   end do
end subroutine

subroutine f_kernel_vec(r, n_f, f_r, n_g, g_r, f_var, len_scale_sq, periodicity, kernel, vec)

   implicit none

   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: vec(:)

   integer :: n_dof

   n_dof=size(r)

   if (n_f > 0) call kernel(r(1:n_dof), f_r(1:n_dof,1:n_f), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                            f_f=vec(1:n_f))
   if (n_g > 0) call kernel(r(1:n_dof), g_r(1:n_dof,1:n_g), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                            f_g=vec(n_f+1:n_f+n_g*n_dof))
end subroutine f_kernel_vec

subroutine g_kernel_vec(r, n_f, f_r, n_g, g_r, f_var, len_scale_sq, periodicity, kernel, vec)

   implicit none

   real(dp), intent(in) :: r(:)
   integer, intent(in) :: n_f
   real(dp), intent(in) :: f_r(:,:)
   integer, intent(in) :: n_g
   real(dp), intent(in) :: g_r(:,:)
   real(dp), intent(in) :: f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   interface
      subroutine kernel(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)
         integer, parameter :: dp = 8
         real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
         real(dp), intent(in) :: periodicity(:)
         real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)
      end subroutine kernel
   end interface
   real(dp), intent(out) :: vec(:,:)

   integer n_dof

   n_dof = size(r)

   if (n_f > 0) call kernel(r(1:n_dof), f_r(1:n_dof,1:n_f), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                            g_f=vec(1:n_dof,1:n_f))
   if (n_g > 0) call kernel(r(1:n_dof), g_r(1:n_dof,1:n_g), f_var, len_scale_sq(1:n_dof), periodicity(1:n_dof), &
                            g_g=vec(1:n_dof,n_f+1:n_f+n_g*n_dof))
end subroutine g_kernel_vec

subroutine SE_kernel_r_rr(x1, x2, f_var, len_scale_sq, periodicity, f_f, g_f, f_g, g_g)

   implicit none

   real(dp), intent(in) :: x1(:), x2(:,:), f_var, len_scale_sq(:)
   real(dp), intent(in) :: periodicity(:)
   real(dp), optional, intent(out) :: f_f(:), g_f(:,:), f_g(:), g_g(:,:)

   real(dp), allocatable :: exp_arg(:), dexp_arg_i(:,:), ddexp_arg_i(:,:), dexp_arg_j(:,:)
   integer :: i, j, nv, n_dof


   n_dof = size(x1)
   nv = size(x2,2)

   allocate(exp_arg(nv))
   allocate(dexp_arg_i(n_dof,nv))
   allocate(ddexp_arg_i(n_dof,nv))
   allocate(dexp_arg_j(n_dof,nv))
   exp_arg = 0.0_dp
   do i=1, n_dof
      if (periodicity(i) /= 0.0_dp) then
         exp_arg(:) = exp_arg(:) - 2.0_dp*sin(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)**2/len_scale_sq(i)
      else
         exp_arg(:) = exp_arg(:) - 0.5_dp*(x2(i,:)-x1(i))**2/len_scale_sq(i)
      endif
   end do

   if (present(f_f)) then
      f_f(:) = f_var * exp(exp_arg(:))
   end if
   if (present(g_f) .or. present(g_g)) then
      do i=1, n_dof
         if (periodicity(i) /= 0.0_dp) then
            dexp_arg_i(i,:) = 2.0_dp*sin(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp) &
                            * cos(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)*TWO_PI/periodicity(i) / len_scale_sq(i)
            if (present(g_g)) ddexp_arg_i(i,:) = (-sin(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)**2 &
                            + cos(TWO_PI/periodicity(i)*(x2(i,:)-x1(i))/2.0_dp)**2)*(TWO_PI/periodicity(i))**2/len_scale_sq(i)
         else
            dexp_arg_i(i,:) = (x2(i,:)-x1(i))/len_scale_sq(i)
            if (present(g_g)) ddexp_arg_i(i,:) = 1.0_dp/len_scale_sq(i)
         endif
      end do
   endif
   if (present(f_g) .or. present(g_g)) then
      do j=1, n_dof
         if (periodicity(j) /= 0.0_dp) then
            dexp_arg_j(j,:) = -2.0_dp*sin(TWO_PI/periodicity(j)*(x2(j,:)-x1(j))/2.0_dp) & 
                            * cos(TWO_PI/periodicity(j)*(x2(j,:)-x1(j))/2.0_dp)*TWO_PI/periodicity(j) / len_scale_sq(j)
         else
            dexp_arg_j(j,:) = -(x2(j,:)-x1(j))/len_scale_sq(j)
         endif
      end do
   endif

   if (present(g_f)) then
      do i=1, n_dof
         g_f(i,1:nv) = f_var * exp(exp_arg(:))*(dexp_arg_i(i,:))
      end do
   endif
   if (present(f_g)) then
      do j=1, n_dof
         f_g(j:nv*n_dof:n_dof) = f_var * exp(exp_arg(:))*(dexp_arg_j(j,:))
      end do
   endif

   if (present(g_g)) then
      do i=1, n_dof
         do j=1, n_dof
            if (i /= j) then
               g_g(i,j:nv*n_dof:n_dof) = f_var * exp(exp_arg(:))*(dexp_arg_i(i,:))*(dexp_arg_j(j,:))
            else
               g_g(i,j:nv*n_dof:n_dof) = f_var * exp(exp_arg(:))*(ddexp_arg_i(i,:)-dexp_arg_i(i,:)**2)
            endif
         end do
      end do
   end if

   deallocate(exp_arg)
   deallocate(dexp_arg_i)
   deallocate(ddexp_arg_i)
   deallocate(dexp_arg_j)

end subroutine SE_kernel_r_rr

subroutine GP_Matrix_Initialise(this,matrix)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix

  if(this%initialised) call gp_matrix_finalise(this)

  this%n = size(matrix,1)
  this%m = size(matrix,2)
  allocate(this%matrix(this%n,this%m), this%factor(this%n,this%m), this%s(this%n), this%tau(this%m) )

  this%matrix = matrix
  this%initialised = .true.

end subroutine GP_Matrix_Initialise

subroutine GP_Matrix_Update(this,matrix)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix

  integer :: factorised

  factorised = this%factorised

  if(this%initialised) then
     if( all(shape(matrix) == (/this%n,this%m/)) ) then
        this%matrix = matrix
     else
        call gp_matrix_initialise(this,matrix)
     endif
  else
     call gp_matrix_initialise(this,matrix)
  endif

  select case(factorised)
  case(QR)
     call GP_Matrix_QR_Factorise(this)
  endselect

end subroutine GP_Matrix_Update

subroutine GP_Matrix_Finalise(this)

  implicit none

  type(GP_Matrix), intent(inout) :: this

  if(.not. this%initialised) return

  this%initialised = .false.
  this%equilibrated = .false.
  this%factorised = NOT_FACTORISED

  this%n = 0
  this%m = 0
  if(allocated(this%matrix) ) deallocate(this%matrix)
  if(allocated(this%factor) ) deallocate(this%factor)
  if(allocated(this%s) ) deallocate(this%s)
  if(allocated(this%tau) ) deallocate(this%tau)

end subroutine GP_Matrix_Finalise

subroutine GP_Matrix_Assign(left_this, right_this)

  implicit none

  type(GP_Matrix), intent(inout) :: left_this
  type(GP_Matrix), intent(in) :: right_this

  call gp_matrix_finalise(left_this)

  left_this%initialised = right_this%initialised
  if( .not. left_this%initialised) return

  left_this%n = right_this%n
  left_this%m = right_this%m

  allocate(left_this%matrix(left_this%n,left_this%m))
  allocate(left_this%factor(left_this%n,left_this%m))
  allocate(left_this%s(left_this%n))
  allocate(left_this%tau(left_this%m))

  left_this%matrix = right_this%matrix
  left_this%factorised = right_this%factorised

  if(left_this%factorised) then
     left_this%factor = right_this%factor
     left_this%s = right_this%s
     left_this%tau = right_this%tau
  end if

end subroutine GP_Matrix_Assign

subroutine GP_Matrix_QR_Factorise(this,q,r)

  implicit none

  type(GP_Matrix), intent(inout) :: this         
  real(dp), dimension(:,:), intent(out), optional :: q, r

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: info

  this%factor = this%matrix

  allocate(work(1))
  lwork = -1
  ! enquiry about optimal lwork for factorization
  call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  ! do the QR factorization
  call dgeqrf(this%n, this%m, this%factor, this%n, this%tau, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_Factorise: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  this%factorised = QR

  if( present(q) .or. present(r) ) call GP_Matrix_GetQR(this,q,r)

end subroutine GP_Matrix_QR_Factorise

subroutine GP_Matrix_GetQR(this,q,r)

  implicit none

  type(GP_Matrix), intent(inout) :: this         
  real(dp), dimension(:,:), intent(out), optional :: q, r

  integer :: lwork
  real(dp), allocatable :: work(:)
  integer :: j, info

  if( this%factorised /= QR ) then
     print *,'GP_Matrix_GetQR: not QR-factorised, call GP_Matrix_QR_Factorise first.'
     stop
  endif

  if(present(q)) then
     if (size(q,1) /= this%n .or. size(q,2) /= this%m) then
        print *, "GT_Matrix_GetQR: shape(q) ",shape(q),"does not match",this%n,this%m
     endif
     q = this%factor

     allocate(work(1))
     lwork = -1
     call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
     lwork = nint(work(1))
     deallocate(work)

     allocate(work(lwork))
     call dorgqr(this%n, this%m, this%m, q, this%n, this%tau, work, lwork, info)
     deallocate(work)
  endif

  if(present(r)) then
     if (size(r,1) /= this%n .or. size(r,2) /= this%m) then
        print *, "GP_Matrix_GetQR: shape(r) ",shape(r),"does not match",this%n,this%m
     endif
     r = this%factor(1:this%m,1:this%m)
     do j = 1, this%m - 1
        r(j+1:this%m,j) = 0.0_dp
     enddo
  endif

  if( info /= 0 ) then
     print *,'GP_Matrix_GetQR: ',(info),'-th parameter had an illegal value.'
     stop
  endif

end subroutine GP_Matrix_GetQR

! computes A^{-1}B matrix product using QR factorization of A: A^{-1}B=R^{-1}Q^{T}B
! this contains the factorization of A, matrix is B
subroutine GP_Matrix_QR_Solve_Matrix(this,matrix,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:,:), intent(in) :: matrix
  real(dp), dimension(:,:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: info, i, j, n, o

  integer :: lwork
  real(dp), allocatable :: work(:)

  if(this%factorised == NOT_FACTORISED) then
     call GP_Matrix_QR_Factorise(this)
  elseif(this%factorised /= QR) then
     print *,'GP_Matrix_QR_Solve_Matrix: matrix not QR-factorised'
     stop
  endif

  n = size(matrix,1)
  o = size(matrix,2)
  if (size(result,1) /= this%m .or. size(result,2) /= o) then
     print *, "GP_Matrix_QR_Solve_matrix: shape(result) ",shape(result),"does not match",this%m,o
  endif

  if( n /= this%n ) then
     print *,'GP_Matrix_QR_Solve_Matrix: dimensions of Q and matrix do not match.'
     stop
  endif

  allocate(my_result(n,o))
  my_result = matrix
  lwork = -1
  allocate(work(1))
  ! enquiry about the optimal lwork for calculating Q^{T}B
  call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
  lwork = nint(work(1))
  deallocate(work)

  allocate(work(lwork))
  ! calculate Q^{T}B
  call dormqr('L', 'T', this%n, o, this%m, this%factor, this%n, this%tau, my_result, this%n, work, lwork, info)
  deallocate(work)

  if( info /= 0 ) then
     print *,'GP_Matrix_QR_Solve_Matrix: ',(-info),'-th parameter had an illegal value.'
     stop
  endif

  ! calculate R^{-1}Q^{T}B using back substitution
  do i = 1, o
     do j = this%m, 2, -1
        my_result(j,i) = my_result(j,i)/this%factor(j,j)
        my_result(1:j-1,i) = my_result(1:j-1,i) - my_result(j,i)*this%factor(1:j-1,j)
     enddo
     my_result(1,i) = my_result(1,i) / this%factor(1,1)
  enddo

  result = my_result(1:this%m,:)
  deallocate(my_result)

end subroutine GP_Matrix_QR_Solve_Matrix

subroutine GP_Matrix_QR_Solve_Vector(this,vector,result)

  implicit none

  type(GP_Matrix), intent(inout) :: this
  real(dp), dimension(:), intent(in) :: vector
  real(dp), dimension(:), intent(out) :: result

  real(dp), dimension(:,:), allocatable :: my_result
  integer :: n, m

  n = size(vector)
  m = size(result)

  allocate(my_result(m,1))

  call GP_Matrix_QR_Solve_Matrix(this,reshape(vector,(/n,1/)),my_result)
  result = my_result(:,1)

  deallocate(my_result)

end subroutine GP_Matrix_QR_Solve_Vector

end module abf_gp 
