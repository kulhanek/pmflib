

program cv_math_test

    use cv_math
    implicit none
    ! --------------------------------------------------------------------------

    call test_norm_vec_der
    call test_get_vlen_der
    call test_get_nvangle
    call test_get_vangle
    call test_get_cross_product
    call test_get_vtors
    call test_get_rotmat

contains

! ==============================================================================

subroutine test_norm_vec_der

    implicit none
    real(PMFDP) :: a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: d_v(3)
    real(PMFDP) :: h,v1(3),v2(3)
    integer     :: i,j
    ! ---------------------------------------------------

    a(1) =  0.5
    a(2) =  0.6
    a(3) =  0.7

    d_v(1) =  0.8d0
    d_v(2) =  1.3d0
    d_v(3) = -0.6d0

    d_a(:) = 0.0d0

    call norm_vec_der(a,d_v,d_a)

    h = 1e-5

    nd_a(:) = 0.0
    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call norm_vec(as,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call norm_vec(as,v2)
        do j=1,3
            nd_a(i) = nd_a(i) + d_v(j)*(v2(j) - v1(j)) / ( 2.0d0 * h )
        end do
    end do

    write(*,*) 'test_norm_vec_der'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a

end subroutine test_norm_vec_der

! ==============================================================================

subroutine test_get_vlen_der

    implicit none
    real(PMFDP) :: a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: d_v
    real(PMFDP) :: h,v1,v2
    integer     :: i
    ! ---------------------------------------------------

    a(1) =  0.5
    a(2) =  0.6
    a(3) =  0.7

    d_v    = 1.0d0
    d_a(:) = 0.0d0

    call get_vlen_der(a,d_v,d_a)

    h = 1e-5

    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call get_vlen(as,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call get_vlen(as,v2)
        nd_a(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    write(*,*) 'test_get_vlen_der'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a

end subroutine test_get_vlen_der

! ==============================================================================

subroutine test_get_nvangle

    implicit none
    real(PMFDP) :: ia(3),a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: ib(3),b(3),bs(3),d_b(3),nd_b(3)
    real(PMFDP) :: d_v
    real(PMFDP) :: h,v1,v2
    integer     :: i
    ! ---------------------------------------------------

    ia(1) =  0.5
    ia(2) =  0.6
    ia(3) =  0.7
    call norm_vec(ia,a)

    ib(1) =  0.3
    ib(2) =  0.2
    ib(3) = -0.1
    call norm_vec(ib,b)

    d_v    = 1.0d0
    d_a(:) = 0.0d0
    d_b(:) = 0.0d0

    call get_nvangle_der(a,b,d_v,d_a,d_b)

    h = 1e-5

    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call get_nvangle(as,b,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call get_nvangle(as,b,v2)
        nd_a(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    do i=1,3
        bs(:) = b(:)
        bs(i) = bs(i) - h
        call get_nvangle(a,bs,v1)
        bs(:) = b(:)
        bs(i) = bs(i) + h
        call get_nvangle(a,bs,v2)
        nd_b(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    write(*,*) 'test_get_nvangle'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a
    write(*,'(3(F10.5,1X))') d_b, nd_b, d_b-nd_b

end subroutine test_get_nvangle

! ==============================================================================

subroutine test_get_vangle

    implicit none
    real(PMFDP) :: a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: b(3),bs(3),d_b(3),nd_b(3)
    real(PMFDP) :: d_v
    real(PMFDP) :: h,v1,v2
    integer     :: i
    ! ---------------------------------------------------

    a(1) =  0.5
    a(2) =  0.6
    a(3) =  0.7

    b(1) =  0.3
    b(2) =  0.2
    b(3) = -0.1

    d_v    = 1.0d0
    d_a(:) = 0.0d0
    d_b(:) = 0.0d0

    call get_vangle_der(a,b,d_v,d_a,d_b)

    h = 1e-5

    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call get_vangle(as,b,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call get_vangle(as,b,v2)
        nd_a(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    do i=1,3
        bs(:) = b(:)
        bs(i) = bs(i) - h
        call get_vangle(a,bs,v1)
        bs(:) = b(:)
        bs(i) = bs(i) + h
        call get_vangle(a,bs,v2)
        nd_b(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    write(*,*) 'test_get_vangle'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a
    write(*,'(3(F10.5,1X))') d_b, nd_b, d_b-nd_b

end subroutine test_get_vangle

! ==============================================================================

subroutine test_get_cross_product

    implicit none
    real(PMFDP) :: a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: b(3),bs(3),d_b(3),nd_b(3)
    real(PMFDP) :: d_v(3)
    real(PMFDP) :: h,v1(3),v2(3)
    integer     :: i,j
    ! ---------------------------------------------------

    a(1) =  0.5
    a(2) =  0.6
    a(3) =  0.7

    b(1) =  0.3
    b(2) =  0.2
    b(3) = -0.1

    d_v(1) = -1.2d0
    d_v(2) =  2.3d0
    d_v(3) = -0.3d0

    d_a(:) = 0.0d0
    d_b(:) = 0.0d0

    call get_cross_product_der(a,b,d_v,d_a,d_b)

    h = 1e-5

    nd_a(:) = 0.0d0
    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call get_cross_product(as,b,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call get_cross_product(as,b,v2)
        do j=1,3
            nd_a(i) = nd_a(i) + d_v(j)*(v2(j) - v1(j)) / ( 2.0d0 * h )
        end do
    end do

    nd_b(:) = 0.0d0
    do i=1,3
        bs(:) = b(:)
        bs(i) = bs(i) - h
        call get_cross_product(a,bs,v1)
        bs(:) = b(:)
        bs(i) = bs(i) + h
        call get_cross_product(a,bs,v2)
        do j=1,3
            nd_b(i) = nd_b(i) + d_v(j)*(v2(j) - v1(j)) / ( 2.0d0 * h )
        end do
    end do

    write(*,*) 'test_get_cross_product'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a
    write(*,'(3(F10.5,1X))') d_b, nd_b, d_b-nd_b

end subroutine test_get_cross_product

! ==============================================================================

subroutine test_get_vtors

    implicit none
    real(PMFDP) :: a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: b(3),bs(3),d_b(3),nd_b(3)
    real(PMFDP) :: c(3),cs(3),d_c(3),nd_c(3)
    real(PMFDP) :: d_v
    real(PMFDP) :: h,v1,v2
    integer     :: i
    ! ---------------------------------------------------

    a(1) =  0.5
    a(2) =  0.6
    a(3) =  0.7

    b(1) =  0.3
    b(2) =  0.2
    b(3) = -0.1

    d_v    = 1.0d0
    d_a(:) = 0.0d0
    d_b(:) = 0.0d0
    d_c(:) = 0.0d0

    call get_vtors_der(a,b,c,d_v,d_a,d_b,d_c)

    h = 1e-5

    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call get_vtors(as,b,c,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call get_vtors(as,b,c,v2)
        nd_a(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    do i=1,3
        bs(:) = b(:)
        bs(i) = bs(i) - h
        call get_vtors(a,bs,c,v1)
        bs(:) = b(:)
        bs(i) = bs(i) + h
        call get_vtors(a,bs,c,v2)
        nd_b(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    do i=1,3
        cs(:) = c(:)
        cs(i) = cs(i) - h
        call get_vtors(a,b,cs,v1)
        cs(:) = c(:)
        cs(i) = cs(i) + h
        call get_vtors(a,b,cs,v2)
        nd_c(i) = (v2 - v1) / ( 2.0d0 * h )
    end do

    write(*,*) 'test_get_vtors'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a
    write(*,'(3(F10.5,1X))') d_b, nd_b, d_b-nd_b
    write(*,'(3(F10.5,1X))') d_c, nd_c, d_c-nd_c

end subroutine test_get_vtors

! ==============================================================================

subroutine test_get_rotmat

    implicit none
    real(PMFDP) :: ia(3),a(3),as(3),d_a(3),nd_a(3)
    real(PMFDP) :: g,gs,d_g,nd_g
    real(PMFDP) :: d_v(3,3)
    real(PMFDP) :: h,v1(3,3),v2(3,3)
    integer     :: i,j,k
    ! ---------------------------------------------------

    ia(1) =  0.5
    ia(2) =  0.6
    ia(3) =  0.7
    call norm_vec(ia,a)

    g = 0.23

    d_v(3,3)    = 1.0d0
    d_a(:)      = 0.0d0
    d_g         = 0.0d0

    call get_rotmat_der(a,g,d_v,d_a,d_g)

    h = 1e-5

    nd_a(:) = 0.0d0
    do i=1,3
        as(:) = a(:)
        as(i) = as(i) - h
        call get_rotmat(as,g,v1)
        as(:) = a(:)
        as(i) = as(i) + h
        call get_rotmat(as,g,v2)
        do j=1,3
            do k=1,3
                nd_a(i) = nd_a(i) + d_v(j,k)*(v2(j,k) - v1(j,k)) / ( 2.0d0 * h )
            end do
        end do
    end do

    nd_g = 0.0d0
    gs = g - h
    call get_rotmat(a,gs,v1)
    gs = g + h
    call get_rotmat(a,gs,v2)
    do j=1,3
        do k=1,3
           nd_g = nd_g + d_v(j,k)*(v2(j,k) - v1(j,k)) / ( 2.0d0 * h )
        end do
    end do

    write(*,*) 'test_get_rotmat'
    write(*,'(3(F10.5,1X))') d_a, nd_a, d_a-nd_a
    write(*,'((F10.5,1X))') d_g, nd_g, d_g-nd_g

end subroutine test_get_rotmat

! ==============================================================================

end program
