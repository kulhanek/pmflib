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

module cv_math

use pmf_sizes
use pmf_constants

implicit none

type :: SImpStrData
    real(PMFDP)     :: ingr
    real(PMFDP)     :: xs(3)
    real(PMFDP)     :: xr(3)
    real(PMFDP)     :: eivals(4)
    real(PMFDP)     :: eivecs(4,4)
    real(PMFDP)     :: u(3,3)
    real(PMFDP)     :: o(3)
end type SImpStrData

contains

!===============================================================================
! Subroutine:  get_vlen
! Input vectors:        a
! Output vector size:    v
!===============================================================================

subroutine get_vlen(a,v)

    implicit none
    real(PMFDP)  :: a(3)
    real(PMFDP) :: v
    ! --------------------------------------------
    real(PMFDP) :: dp
    ! --------------------------------------------------------------------------

    dp = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    v = sqrt(dp)

end subroutine get_vlen

!===============================================================================
! Subroutine:  get_vlen_der
!===============================================================================

subroutine get_vlen_der(a,d_v,d_a)

    implicit none
    real(PMFDP)      :: a(3),d_v
    real(PMFDP)   :: d_a(3)
    ! --------------------------------------------
    real(PMFDP) :: dp,v,f
    ! --------------------------------------------------------------------------

    dp = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    v = sqrt(dp)

    if( v .lt. PMF_MEPS ) then
        v = PMF_MEPS
    end if

    f = 1.0d0 / v

    d_a(:) = d_a(:) + d_v*a(:)*f

end subroutine get_vlen_der

!===============================================================================
! Subroutine:  get_vangle
! Input vectors:    a,b
! Output angle:     v
!===============================================================================

subroutine get_vangle(a,b,v)

    implicit none
    real(PMFDP)  :: a(3),b(3)
    real(PMFDP) :: v
    ! --------------------------------------------
    real(PMFDP) :: n_a(3),n_b(3)
    ! --------------------------------------------------------------------------

    call norm_vec(a,n_a)
    call norm_vec(b,n_b)

    call get_nvangle(n_a,n_b,v)

end subroutine get_vangle

!===============================================================================
! Subroutine:  get_vangle_der
!===============================================================================

subroutine get_vangle_der(a,b,d_v,d_a,d_b)

    implicit none
    real(PMFDP)      :: a(3),b(3),d_v
    real(PMFDP)   :: d_a(3),d_b(3)
    ! --------------------------------------------
    real(PMFDP) :: n_a(3),n_b(3),d_na(3),d_nb(3)
    ! --------------------------------------------------------------------------

!    call norm_vec(a,n_a)
!    call norm_vec(b,n_b)
!
!    call get_nvangle(n_a,n_b,v)

    call norm_vec(a,n_a)
    call norm_vec(b,n_b)

    d_na(:) = 0.0d0
    d_nb(:) = 0.0d0
    call get_nvangle_der(n_a,n_b,d_v,d_na,d_nb)

    call norm_vec_der(a,d_na,d_a)
    call norm_vec_der(b,d_nb,d_b)

end subroutine get_vangle_der

!===============================================================================
! Subroutine:  get_nvangle
! Input vectors:    a,b (normalized vecotrs)
! Output angle:     v
!===============================================================================

subroutine get_nvangle(a,b,v)

    implicit none
    real(PMFDP)  :: a(3),b(3)
    real(PMFDP) :: v
    ! --------------------------------------------
    real(PMFDP) :: dp
    ! --------------------------------------------------------------------------

    dp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

    if ( dp .gt.  1.0d0 ) then
        dp =  1.0d0
    else if ( dp .lt. -1.0d0 ) then
        dp = -1.0d0
    end if

    v = acos(dp)

end subroutine get_nvangle

!===============================================================================
! Subroutine:  get_nvangle_der
!===============================================================================

subroutine get_nvangle_der(a,b,d_v,d_a,d_b)

    implicit none
    real(PMFDP)      :: a(3),b(3),d_v
    real(PMFDP)   :: d_a(3),d_b(3)
    ! --------------------------------------------
    real(PMFDP) :: dp,df,f
    ! --------------------------------------------------------------------------

    dp  = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

    if ( dp .gt.  1.0d0 ) then
        dp =  1.0d0
    else if ( dp .lt. -1.0d0 ) then
        dp = -1.0d0
    end if

    ! v = acos(dp)

    df = sqrt( 1.0d0 - dp**2 )
    if( df .lt. PMF_MEPS ) then
        df = PMF_MEPS
    end if

    f = - 1.0d0 / df

    d_a(:) = d_a(:) + d_v*b(:)*f
    d_b(:) = d_b(:) + d_v*a(:)*f

end subroutine get_nvangle_der

!===============================================================================
! Subroutine:  get_vtors
! Input vectors:        a,b,c ! c - is the rotational axis !!!!!
! Output torsion angle: v
!===============================================================================

subroutine get_vtors(a,b,c,v)

    implicit none
    real(PMFDP)  :: a(3),b(3),c(3)
    real(PMFDP) :: v
    ! --------------------------------------------
    real(PMFDP) :: sc
    real(PMFDP) :: c_a(3)
    real(PMFDP) :: c_b(3)
    ! --------------------------------------------------------------------------

    ! get normal vectors
    call get_cross_product(c,a,c_a)
    call get_cross_product(c,b,c_b)

    ! get angle
    call get_vangle(c_a,c_b,v)

    ! determine sign
    call get_vtors_sign(a,b,c,sc)
    v = v * sc

end subroutine get_vtors

!===============================================================================
! Subroutine:  get_vtors_der
!===============================================================================

subroutine get_vtors_der(a,b,c,d_v,d_a,d_b,d_c)

    implicit none
    real(PMFDP)      :: a(3),b(3),c(3),d_v
    real(PMFDP)   :: d_a(3),d_b(3),d_c(3)
    ! --------------------------------------------
    real(PMFDP) :: sc
    real(PMFDP) :: c_a(3),c_b(3),ta(3),tb(3)
    ! --------------------------------------------------------------------------

   ! get normal vectors
    call get_cross_product(c,a,c_a)
    call get_cross_product(c,b,c_b)
!
!    ! get angle
!    call get_vangle(c_a,c_b,v)

    ! determine sign
    call get_vtors_sign(a,b,c,sc)

    ta(:) = 0.0d0
    tb(:) = 0.0d0

    call get_vangle_der(c_a,c_b,d_v*sc,ta,tb)

    call get_cross_product_der(c,a,ta,d_c,d_a)
    call get_cross_product_der(c,b,tb,d_c,d_b)

end subroutine get_vtors_der

!===============================================================================
! Subroutine:  get_vtors_sign
! Input vectors:        a,b,c
! Output torsion angle: v
!===============================================================================

subroutine get_vtors_sign(a,b,c,sc)

    implicit none
    real(PMFDP)  :: a(3),b(3),c(3)
    real(PMFDP) :: sc
    ! --------------------------------------------
    real(PMFDP) :: c_ab(3), dt
    ! --------------------------------------------------------------------------

    ! get normal vectors
    call get_cross_product(a,b,c_ab)
    dt = c_ab(1)*c(1) + c_ab(2)*c(2) + c_ab(3)*c(3)
    sc = sign(1.0d0,dt)

end subroutine get_vtors_sign

!===============================================================================
! Subroutine:  get_cross_product
!===============================================================================

subroutine get_cross_product(a,b,c)

    implicit none
    real(PMFDP)  :: a(3),b(3)
    real(PMFDP) :: c(3)
    ! --------------------------------------------------------------------------

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine get_cross_product

!===============================================================================
! Subroutine:  get_cross_product_der
!===============================================================================

subroutine get_cross_product_der(a,b,d_c,d_a,d_b)

    implicit none
    real(PMFDP)      :: a(3),b(3),d_c(3)
    real(PMFDP)   :: d_a(3),d_b(3)
    ! --------------------------------------------------------------------------

!    c(1) = a(2)*b(3) - a(3)*b(2)
!    c(2) = a(3)*b(1) - a(1)*b(3)
!    c(3) = a(1)*b(2) - a(2)*b(1)

    d_a(1) = d_a(1) + d_c(3)*b(2) - d_c(2)*b(3)
    d_a(2) = d_a(2) + d_c(1)*b(3) - d_c(3)*b(1)
    d_a(3) = d_a(3) + d_c(2)*b(1) - d_c(1)*b(2)

    d_b(1) = d_b(1) + d_c(2)*a(3) - d_c(3)*a(2)
    d_b(2) = d_b(2) + d_c(3)*a(1) - d_c(1)*a(3)
    d_b(3) = d_b(3) + d_c(1)*a(2) - d_c(2)*a(1)

end subroutine get_cross_product_der

!===============================================================================
! Subroutine:  norm_vec
!===============================================================================

subroutine norm_vec(a,na)

    use pmf_utils
    use pmf_constants

    implicit none
    real(PMFDP)  :: a(3)
    real(PMFDP) :: na(3)
    ! --------------------------------------------
    real(PMFDP) :: a2,v,dp
    ! --------------------------------------------------------------------------

    a2 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    v = sqrt(a2)
    if( v .lt. PMF_MEPS ) then
        call pmf_utils_exit(PMF_OUT,1,'norm_vec - vector has nearly zero length!')
    end if
    dp = 1.0d0/v

    na(:) = a(:)*dp

end subroutine norm_vec

!===============================================================================
! Subroutine:  norm_vec_der
!===============================================================================

subroutine norm_vec_der(a,d_na,d_a)

    use pmf_utils
    use pmf_constants

    implicit none
    real(PMFDP)      :: a(3)
    real(PMFDP)      :: d_na(3)
    real(PMFDP)   :: d_a(3)
    ! --------------------------------------------
    real(PMFDP) :: a2,v,dp,ds,pre
    ! --------------------------------------------------------------------------

    a2 = a(1)*a(1) + a(2)*a(2) + a(3)*a(3)
    v = sqrt(a2)
    if( v .lt. PMF_MEPS ) then
        call pmf_utils_exit(PMF_OUT,1,'norm_vec_der - vector has nearly zero length!')
    end if

    dp = 1.0d0 / v
    ds = dp / a2

! na(:) = a(:)*dp
    ! na(1) = a(1)*dp
    ! na(2) = a(2)*dp
    ! na(3) = a(3)*dp

    ! derivatives
    pre = (d_na(1)*a(1) + d_na(2)*a(2) + d_na(3)*a(3))*ds
    d_a(1) = d_a(1) + d_na(1)*dp - pre*a(1)
    d_a(2) = d_a(2) + d_na(2)*dp - pre*a(2)
    d_a(3) = d_a(3) + d_na(3)*dp - pre*a(3)

end subroutine norm_vec_der

!===============================================================================
! Subroutine:  get_mst
!===============================================================================

subroutine get_mst(a,b,c)

    use pmf_constants

    implicit none
    real(PMFDP)  :: a(3,3),b(3,3)
    real(PMFDP) :: c(3,3)
    ! --------------------------------------------
    real(PMFDP) :: l(3,3)
    ! --------------------------------------------------------------------------

    ! calculate average for each axis
    l(:,1) = (a(:,1) + b(:,1))*0.5d0
    l(:,2) = (a(:,2) + b(:,2))*0.5d0
    l(:,3) = (a(:,3) + b(:,3))*0.5d0

    ! normalize
    call norm_vec(l(:,1),c(:,1))
    call norm_vec(l(:,2),c(:,2))
    call norm_vec(l(:,3),c(:,3))

end subroutine get_mst

!===============================================================================
! Subroutine:  get_mst_der
!===============================================================================

subroutine get_mst_der(a,b,d_c,d_a,d_b)

    use pmf_constants

    implicit none
    real(PMFDP)      :: a(3,3),b(3,3),d_c(3,3)
    real(PMFDP)   :: d_a(3,3),d_b(3,3)
    ! --------------------------------------------
    real(PMFDP) :: l(3,3)
    real(PMFDP) :: d_n(3,3)
    ! --------------------------------------------------------------------------

    ! calculate average for each axis
    l(:,1) = (a(:,1) + b(:,1))*0.5d0
    l(:,2) = (a(:,2) + b(:,2))*0.5d0
    l(:,3) = (a(:,3) + b(:,3))*0.5d0

!    ! normalize
!    call norm_vec(c(:,1))
!    call norm_vec(c(:,2))
!    call norm_vec(c(:,3))

    ! normalize
    d_n(:,:) = 0.0d0
    call norm_vec_der(l(:,1),d_c(:,1),d_n(:,1))
    call norm_vec_der(l(:,2),d_c(:,2),d_n(:,2))
    call norm_vec_der(l(:,3),d_c(:,3),d_n(:,3))

!    ! calculate average for each axis
!    c(:,1) = (a(:,1) + b(:,1))*0.5d0
!    c(:,2) = (a(:,2) + b(:,2))*0.5d0
!    c(:,3) = (a(:,3) + b(:,3))*0.5d0

    d_a(:,1) = d_a(:,1) + 0.5d0*d_n(:,1)
    d_a(:,2) = d_a(:,2) + 0.5d0*d_n(:,2)
    d_a(:,3) = d_a(:,3) + 0.5d0*d_n(:,3)

    d_b(:,1) = d_b(:,1) + 0.5d0*d_n(:,1)
    d_b(:,2) = d_b(:,2) + 0.5d0*d_n(:,2)
    d_b(:,3) = d_b(:,3) + 0.5d0*d_n(:,3)

end subroutine get_mst_der

!===============================================================================
! Subroutine:  get_rotmat
!===============================================================================

subroutine get_rotmat(h,g,r)

    use pmf_constants

    implicit none
    real(PMFDP)  :: h(3)    ! this must be a normal vector
    real(PMFDP)  :: g
    real(PMFDP) :: r(3,3)
    ! --------------------------------------------
    real(PMFDP) :: c, s, dc
    ! --------------------------------------------------------------------------

    ! https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    c = cos(g)
    s = sin(g)
    dc = 1.0d0 - c

    r(1,1) = c + dc * h(1) * h(1)
    r(1,2) = h(1) * h(2) * dc - h(3) * s
    r(1,3) = h(1) * h(3) * dc + h(2) * s
    r(2,1) = h(1) * h(2) * dc + h(3) * s
    r(2,2) = c + dc * h(2) * h(2)
    r(2,3) = h(2) * h(3) * dc - h(1) * s
    r(3,1) = h(1) * h(3) * dc - h(2) * s
    r(3,2) = h(2) * h(3) * dc + h(1) * s
    r(3,3) = c + dc * h(3) * h(3)

end subroutine get_rotmat

!===============================================================================
! Subroutine:  get_rotmat_der
!===============================================================================

subroutine get_rotmat_der(h,g,d_r,d_h,d_g)

    use pmf_constants

    implicit none
    real(PMFDP)      :: h(3)    ! this must be a normal vector
    real(PMFDP)      :: g
    real(PMFDP)      :: d_r(3,3)
    real(PMFDP)   :: d_h(3)
    real(PMFDP)   :: d_g
    ! --------------------------------------------
    real(PMFDP) :: c, s, dc, d_c, d_s, d_dc
    ! --------------------------------------------------------------------------

    ! https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    c = cos(g)
    s = sin(g)
    dc = 1.0d0 - c

!    r(1,1) = c + dc * h(1) * h(1)
!    r(1,2) = h(1) * h(2) * dc - h(3) * s
!    r(1,3) = h(1) * h(3) * dc + h(2) * s

!    r(2,1) = h(1) * h(2) * dc + h(3) * s
!    r(2,2) = c + dc * h(2) * h(2)
!    r(2,3) = h(2) * h(3) * dc - h(1) * s

!    r(3,1) = h(1) * h(3) * dc - h(2) * s
!    r(3,2) = h(2) * h(3) * dc + h(1) * s
!    r(3,3) = c + dc * h(3) * h(3)

    d_h(1) = d_h(1) + 2.0d0*d_r(1,1)*h(1)*dc       + d_r(1,2)*h(2)*dc        + d_r(1,3)*h(3)*dc        &
                          + d_r(2,1)*h(2)*dc                                 - d_r(2,3)*s              &
                          + d_r(3,1)*h(3)*dc       + d_r(3,2)*s

    d_h(2) = d_h(2)                                + d_r(1,2)*h(1)*dc        + d_r(1,3)*s              &
                          + d_r(2,1)*h(1)*dc + 2.0d0*d_r(2,2)*h(2)*dc        + d_r(2,3)*h(3)*dc        &
                          - d_r(3,1)*s             + d_r(3,2)*h(3)*dc

    d_h(3) = d_h(3)                                - d_r(1,2)*s              + d_r(1,3)*h(1)*dc        &
                          + d_r(2,1)*s                                       + d_r(2,3)*h(2)*dc        &
                          + d_r(3,1)*h(1)*dc       + d_r(3,2)*h(2)*dc  + 2.0d0*d_r(3,3)*h(3)*dc

    d_c  =   d_r(1,1) + d_r(2,2) + d_r(3,3)
    d_s  = - d_r(1,2)*h(3) + d_r(1,3)*h(2) + d_r(2,1)*h(3) - d_r(2,3)*h(1) - d_r(3,1)*h(2) + d_r(3,2)*h(1)
    d_dc = + d_r(1,1)*h(1)*h(1) + d_r(1,2)*h(1)*h(2) + d_r(1,3)*h(1)*h(3) &
           + d_r(2,1)*h(2)*h(1) + d_r(2,2)*h(2)*h(2) + d_r(2,3)*h(2)*h(3) &
           + d_r(3,1)*h(3)*h(1) + d_r(3,2)*h(3)*h(2) + d_r(3,3)*h(3)*h(3)
    d_c  = d_c - d_dc

    d_g = d_g - d_c*sin(g) + d_s*cos(g)


end subroutine get_rotmat_der

!===============================================================================
! Subroutine:  rotate_ux
!===============================================================================

subroutine rotate_vec(r,a,b)

    use pmf_constants

    implicit none
    real(PMFDP) :: r(3,3)
    real(PMFDP) :: a(3)
    real(PMFDP):: b(3)
    ! --------------------------------------------------------------------------

    b(1) = r(1,1)*a(1) + r(1,2)*a(2) + r(1,3)*a(3)
    b(2) = r(2,1)*a(1) + r(2,2)*a(2) + r(2,3)*a(3)
    b(3) = r(3,1)*a(1) + r(3,2)*a(2) + r(3,3)*a(3)

end subroutine rotate_vec

!===============================================================================
! Subroutine:  rotate_ux
!===============================================================================

subroutine rotate_vec_der(r,a,d_b,d_r,d_a)

    use pmf_constants

    implicit none
    real(PMFDP)      :: r(3,3),a(3),d_b(3)
    real(PMFDP)   :: d_r(3,3),d_a(3)
    ! --------------------------------------------------------------------------

!    b(1) = r(1,1)*a(1) + r(1,2)*a(2) + r(1,3)*a(3)
!    b(2) = r(2,1)*a(1) + r(2,2)*a(2) + r(2,3)*a(3)
!    b(3) = r(3,1)*a(1) + r(3,2)*a(2) + r(3,3)*a(3)

    d_a(1) = d_a(1) + d_b(1)*r(1,1) + d_b(2)*r(2,1) + d_b(3)*r(3,1)
    d_a(2) = d_a(2) + d_b(1)*r(1,2) + d_b(2)*r(2,2) + d_b(3)*r(3,2)
    d_a(3) = d_a(3) + d_b(1)*r(1,3) + d_b(2)*r(2,3) + d_b(3)*r(3,3)

    d_r(1,1) = d_r(1,1) + d_b(1)*a(1)
    d_r(1,2) = d_r(1,2) + d_b(1)*a(2)
    d_r(1,3) = d_r(1,3) + d_b(1)*a(3)

    d_r(2,1) = d_r(2,1) + d_b(2)*a(1)
    d_r(2,2) = d_r(2,2) + d_b(2)*a(2)
    d_r(2,3) = d_r(2,3) + d_b(2)*a(3)

    d_r(3,1) = d_r(3,1) + d_b(3)*a(1)
    d_r(3,2) = d_r(3,2) + d_b(3)*a(2)
    d_r(3,3) = d_r(3,3) + d_b(3)*a(3)

end subroutine rotate_vec_der

!===============================================================================
! Subroutine:  rotate_ux
!===============================================================================

subroutine rotate_ux(h,g,u,ru)

    implicit none
    real(PMFDP)  :: h(3)
    real(PMFDP)  :: g
    real(PMFDP)  :: u(3,3)
    real(PMFDP) :: ru(3,3)
    ! --------------------------------------------
    real(PMFDP) :: rm(3,3)
    ! --------------------------------------------------------------------------

    call get_rotmat(h,g,rm)

    call rotate_vec(rm,u(:,1),ru(:,1))
    call rotate_vec(rm,u(:,2),ru(:,2))
    call rotate_vec(rm,u(:,3),ru(:,3))

end subroutine rotate_ux

!===============================================================================
! Subroutine:  rotate_ux_der
!===============================================================================

subroutine rotate_ux_der(h,g,u,d_ru,d_h,d_g,d_u)

    implicit none
    real(PMFDP)      :: h(3)
    real(PMFDP)      :: g
    real(PMFDP)      :: u(3,3)
    real(PMFDP)      :: d_ru(3,3)
    real(PMFDP)   :: d_h(3)
    real(PMFDP)   :: d_g
    real(PMFDP)   :: d_u(3,3)
    ! --------------------------------------------
    real(PMFDP) :: rm(3,3), d_rm(3,3)
    ! --------------------------------------------------------------------------

    call get_rotmat(h,g,rm)

!    call rotate_vec(rm,u(:,1),ru(:,1))
!    call rotate_vec(rm,u(:,2),ru(:,2))
!    call rotate_vec(rm,u(:,3),ru(:,3))

    d_rm(:,:) = 0.0d0

    call rotate_vec_der(rm,u(:,1),d_ru(:,1),d_rm,d_u(:,1))
    call rotate_vec_der(rm,u(:,2),d_ru(:,2),d_rm,d_u(:,2))
    call rotate_vec_der(rm,u(:,3),d_ru(:,3),d_rm,d_u(:,3))

    call get_rotmat_der(h,g,d_rm,d_h,d_g)

end subroutine rotate_ux_der

!===============================================================================
! Subroutine:  superimpose_str
!===============================================================================

subroutine superimpose_str(cv_item,gi,gj,x,str,simpdat,u,o)

    use pmf_dat
    use pmf_utils
    use smf_xyzfile
    use smf_xyzfile_type

    implicit none
    class(CVType)     :: cv_item
    integer              :: gi,gj
    real(PMFDP)          :: x(:,:)
    type(XYZFILE_TYPE)   :: str
    class(SImpStrData)  :: simpdat
    real(PMFDP)         :: u(3,3)
    real(PMFDP)         :: o(3)
    ! -----------------------------------------------
    integer             :: i,ai,info,best
    real(PMFDP)         :: work(26*4),tmp1(3)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    ! --------------------------------------------------------------------------

    ! reset data
    simpdat%ingr        = 0.0d0
    simpdat%xs(:)       = 0.0d0
    simpdat%xr(:)       = 0.0d0
    simpdat%eivals(:)   = 0.0d0
    simpdat%eivecs(:,:) = 0.0d0
    simpdat%u(:,:)      = 0.0d0
    simpdat%o(:)        = 0.0d0

    ! inverse number of atoms
    simpdat%ingr = 1.0d0 / (gj-gi)

    ! calculate geometrical centres (source and target) -------------------

    do  i = gi+1,gj
        ai = cv_item%lindexes(i)
        ! source
        simpdat%xs(:) = simpdat%xs(:) + x(:,ai)

        ! reference
        simpdat%xr(:) = simpdat%xr(:) + str%cvs(:,i-gi)
    end do

    simpdat%xs(:) = simpdat%xs(:) * simpdat%ingr
    simpdat%xr(:) = simpdat%xr(:) * simpdat%ingr

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    do i = gi+1,gj
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - simpdat%xs(1))*(str%cvs(1,i-gi) - simpdat%xr(1))
        r12 = r12 + (x(1,ai) - simpdat%xs(1))*(str%cvs(2,i-gi) - simpdat%xr(2))
        r13 = r13 + (x(1,ai) - simpdat%xs(1))*(str%cvs(3,i-gi) - simpdat%xr(3))

        r21 = r21 + (x(2,ai) - simpdat%xs(2))*(str%cvs(1,i-gi) - simpdat%xr(1))
        r22 = r22 + (x(2,ai) - simpdat%xs(2))*(str%cvs(2,i-gi) - simpdat%xr(2))
        r23 = r23 + (x(2,ai) - simpdat%xs(2))*(str%cvs(3,i-gi) - simpdat%xr(3))

        r31 = r31 + (x(3,ai) - simpdat%xs(3))*(str%cvs(1,i-gi) - simpdat%xr(1))
        r32 = r32 + (x(3,ai) - simpdat%xs(3))*(str%cvs(2,i-gi) - simpdat%xr(2))
        r33 = r33 + (x(3,ai) - simpdat%xs(3))*(str%cvs(3,i-gi) - simpdat%xr(3))
    end do

    r11 = r11 * simpdat%ingr
    r12 = r12 * simpdat%ingr
    r13 = r13 * simpdat%ingr

    r21 = r21 * simpdat%ingr
    r22 = r22 * simpdat%ingr
    r23 = r23 * simpdat%ingr

    r31 = r31 * simpdat%ingr
    r32 = r32 * simpdat%ingr
    r33 = r33 * simpdat%ingr

    ! construct matrix for quaterion fitting
    simpdat%eivecs(1,1) =  r11 + r22 + r33
    simpdat%eivecs(1,2) =  r23 - r32
    simpdat%eivecs(1,3) =  r31 - r13
    simpdat%eivecs(1,4) =  r12 - r21

    simpdat%eivecs(2,1) =  r23 - r32
    simpdat%eivecs(2,2) =  r11 - r22 - r33
    simpdat%eivecs(2,3) =  r12 + r21
    simpdat%eivecs(2,4) =  r13 + r31

    simpdat%eivecs(3,1) =  r31 - r13
    simpdat%eivecs(3,2) =  r12 + r21
    simpdat%eivecs(3,3) = -r11 + r22 - r33
    simpdat%eivecs(3,4) =  r23 + r32

    simpdat%eivecs(4,1) =  r12 - r21
    simpdat%eivecs(4,2) =  r13 + r31
    simpdat%eivecs(4,3) =  r23 + r32
    simpdat%eivecs(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    simpdat%eivals(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, simpdat%eivecs, 4, simpdat%eivals, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in superimpose_str!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    simpdat%u(1,1) = simpdat%eivecs(1,best)**2 + simpdat%eivecs(2,best)**2 - simpdat%eivecs(3,best)**2 - simpdat%eivecs(4,best)**2
    simpdat%u(2,1) = 2.0d0*( simpdat%eivecs(2,best)*simpdat%eivecs(3,best) - simpdat%eivecs(1,best)*simpdat%eivecs(4,best) )
    simpdat%u(3,1) = 2.0d0*( simpdat%eivecs(2,best)*simpdat%eivecs(4,best) + simpdat%eivecs(1,best)*simpdat%eivecs(3,best) )

    simpdat%u(1,2) = 2.0d0*( simpdat%eivecs(2,best)*simpdat%eivecs(3,best) + simpdat%eivecs(1,best)*simpdat%eivecs(4,best) )
    simpdat%u(2,2) = simpdat%eivecs(1,best)**2 - simpdat%eivecs(2,best)**2 + simpdat%eivecs(3,best)**2 - simpdat%eivecs(4,best)**2
    simpdat%u(3,2) = 2.0d0*( simpdat%eivecs(3,best)*simpdat%eivecs(4,best) - simpdat%eivecs(1,best)*simpdat%eivecs(2,best) )

    simpdat%u(1,3) = 2.0d0*( simpdat%eivecs(2,best)*simpdat%eivecs(4,best) - simpdat%eivecs(1,best)*simpdat%eivecs(3,best) )
    simpdat%u(2,3) = 2.0d0*( simpdat%eivecs(3,best)*simpdat%eivecs(4,best) + simpdat%eivecs(1,best)*simpdat%eivecs(2,best) )
    simpdat%u(3,3) = simpdat%eivecs(1,best)**2 - simpdat%eivecs(2,best)**2 - simpdat%eivecs(3,best)**2 + simpdat%eivecs(4,best)**2

    ! resulting rotational matrix
    u(:,:) = simpdat%u(:,:)

    ! get origins of bases
    simpdat%o(:) = 0.0d0
    ! move reference point to origin
    simpdat%o(:) = simpdat%o(:) - simpdat%xr(:)
    ! rotate
    tmp1(1) = simpdat%u(1,1)*simpdat%o(1) + simpdat%u(1,2)*simpdat%o(2) + simpdat%u(1,3)*simpdat%o(3)
    tmp1(2) = simpdat%u(2,1)*simpdat%o(1) + simpdat%u(2,2)*simpdat%o(2) + simpdat%u(2,3)*simpdat%o(3)
    tmp1(3) = simpdat%u(3,1)*simpdat%o(1) + simpdat%u(3,2)*simpdat%o(2) + simpdat%u(3,3)*simpdat%o(3)
    ! move origin to new reference point (experimental structure)
    simpdat%o(:) = tmp1(:) + simpdat%xs(:)

    o(:) = simpdat%o(:)

end subroutine superimpose_str

!===============================================================================
! Subroutine:  superimpose_str_der
!===============================================================================

subroutine superimpose_str_der(cv_item,gi,gj,ctx,str,simpdat,a_tu,a_o)

    use pmf_dat
    use pmf_utils
    use smf_xyzfile
    use smf_xyzfile_type

    implicit none
    class(CVType)            :: cv_item
    integer                  :: gi,gj
    type(CVContextType)   :: ctx
    type(XYZFILE_TYPE)       :: str
    class(SImpStrData)       :: simpdat
    real(PMFDP)              :: a_tu(3,3)
    real(PMFDP)              :: a_o(3)
    ! -----------------------------------------------
    integer             :: i,best,mi,mj,ai
    real(PMFDP)         :: a_fa(4),a_rij(4,4)
    real(PMFDP)         :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4)
    real(PMFDP)         :: l_u(3,3),a_xs(3)
    ! --------------------------------------------------------------------------

    l_u(:,:) = a_tu(:,:)

! origin
!    ! get origins of bases
!    o(:) = 0.0d0
!    ! move reference point to origin
!    o(:) = o(:) - simpdat%xr(:)
!    ! rotate
!    tmp1(1) = simpdat%u(1,1)*o(1) + simpdat%u(1,2)*o(2) + simpdat%u(1,3)*o(3)
!    tmp1(2) = simpdat%u(2,1)*o(1) + simpdat%u(2,2)*o(2) + simpdat%u(2,3)*o(3)
!    tmp1(3) = simpdat%u(3,1)*o(1) + simpdat%u(3,2)*o(2) + simpdat%u(3,3)*o(3)
!    ! move origin to new reference point (experimental structure)
!    o(:) = tmp1(:) + simpdat%xs(:)

    l_u(1,1) = l_u(1,1) - simpdat%xr(1)*a_o(1)
    l_u(2,1) = l_u(2,1) - simpdat%xr(1)*a_o(2)
    l_u(3,1) = l_u(3,1) - simpdat%xr(1)*a_o(3)
    l_u(1,2) = l_u(1,2) - simpdat%xr(2)*a_o(1)
    l_u(2,2) = l_u(2,2) - simpdat%xr(2)*a_o(2)
    l_u(3,2) = l_u(3,2) - simpdat%xr(2)*a_o(3)
    l_u(1,3) = l_u(1,3) - simpdat%xr(3)*a_o(1)
    l_u(2,3) = l_u(2,3) - simpdat%xr(3)*a_o(2)
    l_u(3,3) = l_u(3,3) - simpdat%xr(3)*a_o(3)

    a_xs(:) = a_o(:)*simpdat%ingr

    best = 4

! rotation matrix a ------------------------------
!     ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
!     ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
!     ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    a_fa(1) = 2.0d0*( simpdat%eivecs(1,best)*l_u(1,1) - simpdat%eivecs(4,best)*l_u(2,1) + simpdat%eivecs(3,best)*l_u(3,1))
    a_fa(2) = 2.0d0*( simpdat%eivecs(2,best)*l_u(1,1) + simpdat%eivecs(3,best)*l_u(2,1) + simpdat%eivecs(4,best)*l_u(3,1))
    a_fa(3) = 2.0d0*(-simpdat%eivecs(3,best)*l_u(1,1) + simpdat%eivecs(2,best)*l_u(2,1) + simpdat%eivecs(1,best)*l_u(3,1))
    a_fa(4) = 2.0d0*(-simpdat%eivecs(4,best)*l_u(1,1) - simpdat%eivecs(1,best)*l_u(2,1) + simpdat%eivecs(2,best)*l_u(3,1))

!     ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
!     ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
!     ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    a_fa(1) = a_fa(1) + 2.0d0*(simpdat%eivecs(4,best)*l_u(1,2) + simpdat%eivecs(1,best)*l_u(2,2) - simpdat%eivecs(2,best)*l_u(3,2))
    a_fa(2) = a_fa(2) + 2.0d0*(simpdat%eivecs(3,best)*l_u(1,2) - simpdat%eivecs(2,best)*l_u(2,2) - simpdat%eivecs(1,best)*l_u(3,2))
    a_fa(3) = a_fa(3) + 2.0d0*(simpdat%eivecs(2,best)*l_u(1,2) + simpdat%eivecs(3,best)*l_u(2,2) + simpdat%eivecs(4,best)*l_u(3,2))
    a_fa(4) = a_fa(4) + 2.0d0*(simpdat%eivecs(1,best)*l_u(1,2) - simpdat%eivecs(4,best)*l_u(2,2) + simpdat%eivecs(3,best)*l_u(3,2))

!     ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
!     ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
!     ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    a_fa(1) = a_fa(1) + 2.0d0*(-simpdat%eivecs(3,best)*l_u(1,3) + simpdat%eivecs(2,best)*l_u(2,3) + simpdat%eivecs(1,best)*l_u(3,3))
    a_fa(2) = a_fa(2) + 2.0d0*( simpdat%eivecs(4,best)*l_u(1,3) + simpdat%eivecs(1,best)*l_u(2,3) - simpdat%eivecs(2,best)*l_u(3,3))
    a_fa(3) = a_fa(3) + 2.0d0*(-simpdat%eivecs(1,best)*l_u(1,3) + simpdat%eivecs(4,best)*l_u(2,3) - simpdat%eivecs(3,best)*l_u(3,3))
    a_fa(4) = a_fa(4) + 2.0d0*( simpdat%eivecs(2,best)*l_u(1,3) + simpdat%eivecs(3,best)*l_u(2,3) + simpdat%eivecs(4,best)*l_u(3,3))

! derivatives of fa with respect to matrix elements
    v(:,:) = simpdat%eivecs(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(simpdat%eivals(i) - simpdat%eivals(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + simpdat%eivecs(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

! merge xij with a_fa, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fa(1)*xij(1,mi,mj)+a_fa(2)*xij(2,mi,mj)+a_fa(3)*xij(3,mi,mj)+a_fa(4)*xij(4,mi,mj))*simpdat%ingr
        end do
    end do

! finally gradients for group_a
    do i = gi+1, gj

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(str%cvs(1,i-gi) - simpdat%xr(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(str%cvs(2,i-gi) - simpdat%xr(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(str%cvs(3,i-gi) - simpdat%xr(3)) &
                + a_xs(1)

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(str%cvs(1,i-gi) - simpdat%xr(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(str%cvs(2,i-gi) - simpdat%xr(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(str%cvs(3,i-gi) - simpdat%xr(3)) &
                + a_xs(2)

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(str%cvs(1,i-gi) - simpdat%xr(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(str%cvs(2,i-gi) - simpdat%xr(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(str%cvs(3,i-gi) - simpdat%xr(3)) &
                + a_xs(3)
    end do

end subroutine superimpose_str_der

! ==============================================================================

subroutine get_mst_morg(ua,oa,ub,ob,y1,y2,mst,morg)

    implicit none
    real(PMFDP)      :: ua(3,3),ub(3,3)
    real(PMFDP)      :: oa(3),ob(3)
    real(PMFDP)      :: y1(3),y2(3)
    real(PMFDP)     :: mst(3,3)
    real(PMFDP)     :: morg(3)
    ! ---------------------------------------------------
    real(PMFDP)     :: xaxis(3),yaxis(3),zaxis(3)
    real(PMFDP)     :: yaxisr(3),zaxisr(3)
    real(PMFDP)     :: y0axis(3),zsc
    ! --------------------------------------------------------------------------

! z-axis =========================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    call norm_vec(zaxisr,zaxis)

! y-axis =========================================
    y0axis(:) = y1(:) - y2(:)
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    call norm_vec(yaxisr,yaxis)

! x-axis =========================================
    call get_cross_product(yaxis,zaxis,xaxis)

! get origin and complete mst ====================
    morg(:) = 0.5d0*(oa(:) + ob(:))
    mst(:,1) = xaxis(:)
    mst(:,2) = yaxis(:)
    mst(:,3) = zaxis(:)

end subroutine get_mst_morg

! ==============================================================================

subroutine get_mst_morg_der(ua,ub,y1,y2,a_mst,a_morg,a_ua,a_oa,a_ub,a_ob,a_y1,a_y2)

    implicit none
    real(PMFDP)      :: ua(3,3),ub(3,3)
    real(PMFDP)      :: y1(3),y2(3)
    real(PMFDP)      :: a_mst(3,3)
    real(PMFDP)      :: a_morg(3)
    real(PMFDP)   :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)   :: a_oa(3),a_ob(3)
    real(PMFDP)   :: a_y1(3),a_y2(3)
    ! ---------------------------------------------------
    real(PMFDP)     :: yaxis(3),zaxis(3)
    real(PMFDP)     :: yaxisr(3),zaxisr(3)
    real(PMFDP)     :: y0axis(3),zsc
    real(PMFDP)     :: a_xaxis(3),a_yaxis(3),a_zaxis(3)
    real(PMFDP)     :: a_zaxisr(3),a_yaxisr(3),a_y0axis(3)
    real(PMFDP)     :: t1
    ! --------------------------------------------------------------------------

! z-axis =========================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    call norm_vec(zaxisr,zaxis)

! y-axis =========================================
    y0axis(:) = y1(:) - y2(:)
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    call norm_vec(yaxisr,yaxis)
!
!! x-axis =========================================
!    call get_cross_product(yaxis,zaxis,xaxis)
!
!! get origin and complete mst ====================
!    morg(:) = 0.5d0*(oa(:) + ob(:))
!    mst(:,1) = xaxis(:)
!    mst(:,2) = yaxis(:)
!    mst(:,3) = zaxis(:)

! final derivatives
    a_oa(:)    = a_oa(:) + 0.5d0*a_morg(:)
    a_ob(:)    = a_ob(:) + 0.5d0*a_morg(:)

    a_xaxis(:) = a_mst(:,1)
    a_yaxis(:) = a_mst(:,2)
    a_zaxis(:) = a_mst(:,3)

! x-axis
    call get_cross_product_der(yaxis,zaxis,a_xaxis,a_yaxis,a_zaxis)

! y-axis
    call norm_vec_der(yaxisr,a_yaxis,a_yaxisr)

!    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
!    ! remove projections to z-axis
!    yaxisr(1) = y0axis(1) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(1)
!    yaxisr(2) = y0axis(2) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(2)
!    yaxisr(3) = y0axis(3) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(3)

    ! with respect to y0axis
    t1 = zaxis(1)*a_yaxisr(1) + zaxis(2)*a_yaxisr(2) + zaxis(3)*a_yaxisr(3)
    a_y0axis(1) = a_yaxisr(1) - zaxis(1)*t1
    a_y0axis(2) = a_yaxisr(2) - zaxis(2)*t1
    a_y0axis(3) = a_yaxisr(3) - zaxis(3)*t1

    ! with respect to zaxis
    a_zaxis(1) = a_zaxis(1) - a_yaxisr(1) * ( 2.0d0*y0axis(1)*zaxis(1) + y0axis(2)*zaxis(2) + y0axis(3)*zaxis(3) ) &
                            - a_yaxisr(2) * y0axis(1)*zaxis(2) &
                            - a_yaxisr(3) * y0axis(1)*zaxis(3)

    a_zaxis(2) = a_zaxis(2) - a_yaxisr(1) * y0axis(2)*zaxis(1) &
                            - a_yaxisr(2) * ( y0axis(1)*zaxis(1) + 2.0d0*y0axis(2)*zaxis(2) + y0axis(3)*zaxis(3) ) &
                            - a_yaxisr(3) * y0axis(2)*zaxis(3)

    a_zaxis(3) = a_zaxis(3) - a_yaxisr(1) * y0axis(3)*zaxis(1) &
                            - a_yaxisr(2) * y0axis(3)*zaxis(2) &
                            - a_yaxisr(3) * ( y0axis(1)*zaxis(1) + y0axis(2)*zaxis(2) + 2.0d0*y0axis(3)*zaxis(3) )

! z-axis
    call norm_vec_der(zaxisr,a_zaxis,a_zaxisr)
    a_ua(:,3) = a_ua(:,3) + 0.5d0*a_zaxisr(:)
    a_ub(:,3) = a_ub(:,3) + 0.5d0*zsc*a_zaxisr(:)

    a_y1(:) = a_y1(:) + a_y0axis(:)
    a_y2(:) = a_y2(:) - a_y0axis(:)

end subroutine get_mst_morg_der

!===============================================================================

end module cv_math

