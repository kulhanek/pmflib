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

module abf_gp_dat
implicit none

integer, parameter :: dp = 8

real(dp), parameter :: TWO_PI = 2.0_dp*3.14159265358979_dp

integer, parameter :: NOT_FACTORISED = 0
integer, parameter :: QR             = 2

type GP_Matrix
    logical :: initialised = .false.
    logical :: equilibrated = .false.
    integer :: factorised = NOT_FACTORISED
    integer :: n
    integer :: m
    real(dp), dimension(:,:), allocatable :: matrix
    real(dp), dimension(:,:), allocatable :: factor
    real(dp), dimension(:), allocatable :: s
    real(dp), dimension(:), allocatable :: tau
end type GP_Matrix

type gp_basic
    logical :: initialised = .false.
    logical :: sparsified = .false.
    logical :: partial = .false.
    integer :: n_dof
    integer :: m_f
    integer :: m_g
    integer :: m_teach
    real(dp) :: jitter
    real(dp) :: f_var
    real(dp), allocatable :: len_scale_sq(:)
    real(dp), allocatable :: periodicity(:)
    real(dp), allocatable :: f_r(:,:)
    real(dp), allocatable :: g_r(:,:)
    real(dp), allocatable :: k(:)
    real(dp), allocatable :: k_grad(:,:)
    real(dp), allocatable :: mat_inv_k(:)
    real(dp), allocatable :: mat_inv_k_grad(:,:)
    real(dp), allocatable :: Cmat_inv_v(:)
    real(dp), allocatable :: Kmn_siginvsq_y(:)
    real(dp), allocatable :: Kmn_siginvsq_Knm(:,:)
    type(GP_Matrix) :: Kmm
    type(GP_Matrix) :: noise_Kmm
    type(GP_Matrix) :: Cmat
end type gp_basic

end module abf_gp_dat
