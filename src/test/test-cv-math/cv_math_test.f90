

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
    call test_get_mst_morg_simple

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

    g = 2.98

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

subroutine test_get_mst_morg_simple

    implicit none
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: oa(3),ob(3)
    real(PMFDP)         :: y1(3),y2(3)
    real(PMFDP)         :: a_ua(3,3),a_ub(3,3)
    real(PMFDP)         :: a_oa(3),a_ob(3)
    real(PMFDP)         :: a_y1(3),a_y2(3)
    real(PMFDP)         :: r,h
    real(PMFDP)         :: a_mst(3,3),a_morg(3)
    real(PMFDP)         :: mst1(3,3),morg1(3)
    real(PMFDP)         :: mst2(3,3),morg2(3)
    real(PMFDP)         :: na_oa(3),na_y1(3),na_y2(3),na_ua(3,3),na_ub(3,3),na_ob(3)
    real(PMFDP)         :: oas(3),y1s(3),y2s(3),uas(3,3),ubs(3,3),obs(3)
    integer             :: i,j,k,l
    ! ---------------------------------------------------

    do i=1,3
        do j=1,3
            call RANDOM_NUMBER(r)
            ua(i,j) = r
            call RANDOM_NUMBER(r)
            ub(i,j) = r
            call RANDOM_NUMBER(r)
            a_mst(i,j)  = r
        end do
        call RANDOM_NUMBER(r)
        y1(i) = r
        call RANDOM_NUMBER(r)
        y2(i) = r
        call RANDOM_NUMBER(r)
        oa(i) = r
        call RANDOM_NUMBER(r)
        ob(i) = r
        call RANDOM_NUMBER(r)
        a_morg(i)  = r
    end do

    a_ua(:,:) = 0.0d0
    a_oa(:) = 0.0d0
    a_ub(:,:) = 0.0d0
    a_ob(:) = 0.0d0
    a_y1(:) = 0.0d0
    a_y2(:) = 0.0d0

    call get_mst_morg_simple_der(ua,ub,y1,y2,a_mst,a_morg,a_ua,a_oa,a_ub,a_ob,a_y1,a_y2)

    h = 1e-5

    na_oa(:) = 0.0d0
    do i=1,3
        oas(:) = oa(:)
        oas(i) = oas(i) - h
        call get_mst_morg_simple(ua,oas,ub,ob,y1,y2,mst1,morg1)
        oas(:) = oa(:)
        oas(i) = oas(i) + h
        call get_mst_morg_simple(ua,oas,ub,ob,y1,y2,mst2,morg2)
        do j=1,3
            do k=1,3
                na_oa(i) = na_oa(i) + a_mst(j,k)*(mst2(j,k) - mst1(j,k)) / ( 2.0d0 * h )
            end do
        end do
        do j=1,3
            na_oa(i) = na_oa(i) + a_morg(j)*(morg2(j) - morg1(j)) / ( 2.0d0 * h )
        end do
    end do

    na_ob(:) = 0.0d0
    do i=1,3
        obs(:) = ob(:)
        obs(i) = obs(i) - h
        call get_mst_morg_simple(ua,oa,ub,obs,y1,y2,mst1,morg1)
        obs(:) = ob(:)
        obs(i) = obs(i) + h
        call get_mst_morg_simple(ua,oa,ub,obs,y1,y2,mst2,morg2)
        do j=1,3
            do k=1,3
                na_ob(i) = na_ob(i) + a_mst(j,k)*(mst2(j,k) - mst1(j,k)) / ( 2.0d0 * h )
            end do
        end do
        do j=1,3
            na_ob(i) = na_ob(i) + a_morg(j)*(morg2(j) - morg1(j)) / ( 2.0d0 * h )
        end do
    end do

    na_y1(:) = 0.0d0
    do i=1,3
        y1s(:) = y1(:)
        y1s(i) = y1s(i) - h
        call get_mst_morg_simple(ua,oa,ub,ob,y1s,y2,mst1,morg1)
        y1s(:) = y1(:)
        y1s(i) = y1s(i) + h
        call get_mst_morg_simple(ua,oa,ub,ob,y1s,y2,mst2,morg2)
        do j=1,3
            do k=1,3
                na_y1(i) = na_y1(i) + a_mst(j,k)*(mst2(j,k) - mst1(j,k)) / ( 2.0d0 * h )
            end do
        end do
        do j=1,3
            na_y1(i) = na_y1(i) + a_morg(j)*(morg2(j) - morg1(j)) / ( 2.0d0 * h )
        end do
    end do

    na_y2(:) = 0.0d0
    do i=1,3
        y2s(:) = y2(:)
        y2s(i) = y2s(i) - h
        call get_mst_morg_simple(ua,oa,ub,ob,y1,y2s,mst1,morg1)
        y2s(:) = y2(:)
        y2s(i) = y2s(i) + h
        call get_mst_morg_simple(ua,oa,ub,ob,y1,y2s,mst2,morg2)
        do j=1,3
            do k=1,3
                na_y2(i) = na_y2(i) + a_mst(j,k)*(mst2(j,k) - mst1(j,k)) / ( 2.0d0 * h )
            end do
        end do
        do j=1,3
            na_y2(i) = na_y2(i) + a_morg(j)*(morg2(j) - morg1(j)) / ( 2.0d0 * h )
        end do
    end do

    na_ua(:,:) = 0.0d0
    do i=1,3
        do l=1,3
            uas(:,:) = ua(:,:)
            uas(i,l) = uas(i,l) - h
            call get_mst_morg_simple(uas,oa,ub,ob,y1,y2,mst1,morg1)
            uas(:,:) = ua(:,:)
            uas(i,l) = uas(i,l) + h
            call get_mst_morg_simple(uas,oa,ub,ob,y1,y2,mst2,morg2)
            do j=1,3
                do k=1,3
                    na_ua(i,l) = na_ua(i,l) + a_mst(j,k)*(mst2(j,k) - mst1(j,k)) / ( 2.0d0 * h )
                end do
            end do
            do j=1,3
                na_ua(i,l) = na_ua(i,l) + a_morg(j)*(morg2(j) - morg1(j)) / ( 2.0d0 * h )
            end do
        end do
    end do

    na_ub(:,:) = 0.0d0
    do i=1,3
        do l=1,3
            ubs(:,:) = ub(:,:)
            ubs(i,l) = ubs(i,l) - h
            call get_mst_morg_simple(ua,oa,ubs,ob,y1,y2,mst1,morg1)
            ubs(:,:) = ub(:,:)
            ubs(i,l) = ubs(i,l) + h
            call get_mst_morg_simple(ua,oa,ubs,ob,y1,y2,mst2,morg2)
            do j=1,3
                do k=1,3
                    na_ub(i,l) = na_ub(i,l) + a_mst(j,k)*(mst2(j,k) - mst1(j,k)) / ( 2.0d0 * h )
                end do
            end do
            do j=1,3
                na_ub(i,l) = na_ub(i,l) + a_morg(j)*(morg2(j) - morg1(j)) / ( 2.0d0 * h )
            end do
        end do
    end do

    write(*,*) '>>>> test_get_mst_morg_simple'
    write(*,*) '> oa'
    write(*,'(3(F10.5,1X))') a_oa, na_oa, a_oa-na_oa
    write(*,*) '> ob'
    write(*,'(3(F10.5,1X))') a_ob, na_ob, a_ob-na_ob
    write(*,*) '> y1'
    write(*,'(3(F10.5,1X))') a_y1, na_y1, a_y1-na_y1
    write(*,*) '> y2'
    write(*,'(3(F10.5,1X))') a_y2, na_y2, a_y2-na_y2
    write(*,*) '> ua'
    write(*,'(3(F10.5,1X))') a_ua, na_ua, a_ua-na_ua
    write(*,*) '> ub'
    write(*,'(3(F10.5,1X))') a_ub, na_ub, a_ub-na_ub

end subroutine test_get_mst_morg_simple


! ==============================================================================

end program
