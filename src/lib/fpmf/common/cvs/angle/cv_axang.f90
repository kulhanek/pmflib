!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@kulhanek@chemi.muni.cz
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

! angle between two principal axes - AXANG

module cv_axang

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeAXANG

    integer :: direction_a
    integer :: head_a
    integer :: direction_b
    integer :: head_b

    contains
        procedure :: load_cv        => load_axang
        procedure :: calculate_cv   => calculate_axang
end type CVTypeAXANG

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_axang
!===============================================================================

subroutine load_axang(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common

    implicit none
    class(CVTypeAXANG)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: m
    character(1)                        :: cdir
    logical                             :: found
    character(len=PRMFILE_MAX_LINE)     :: mask
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'AXANG'
    cv_item%unit          = AngleUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_init_groups(cv_item,prm_fin)

    ! load axang specific data  ----------------------

    ! load axis a -----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    if( .not. prmfile_get_string_by_key(prm_fin,'direction_a',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'direction_a is not specified!')
    end if

    write(PMF_OUT,60) cdir

    select case(cdir)
        case('x')
            cv_item%direction_a = 1
        case('y')
            cv_item%direction_a = 2
        case('z')
            cv_item%direction_a = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,'direction_a has to be x, y, or z!')
    end select

    ! load cv_item%head_a -----------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'head_a',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'head_a atom is not specified!')
    end if
    write(PMF_OUT,65) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%head_a)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=1,cv_item%grps(1)
        if( cv_item%rindexes(m) .eq. cv_item%head_a ) then
            cv_item%head_a = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'head_a atom has to be part of group_a!')
    end if

    ! load axis b -----------------------------------
    write(PMF_OUT,80)
    call cv_common_read_group(cv_item,prm_fin,2)

    if( .not. prmfile_get_string_by_key(prm_fin,'direction_b',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'direction_b is not specified!')
    end if

    write(PMF_OUT,90) cdir

    select case(cdir)
        case('x')
            cv_item%direction_b = 1
        case('y')
            cv_item%direction_b = 2
        case('z')
            cv_item%direction_b = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,'direction_b has to be x, y, or z!')
    end select

    ! load cv_item%head_a -----------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'head_b',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'head_b atom is not specified!')
    end if
    write(PMF_OUT,95) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%head_b)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=cv_item%grps(1)+1,cv_item%grps(2)
        if( cv_item%rindexes(m) .eq. cv_item%head_b ) then
            cv_item%head_b = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'head_b atom has to be part of group_b!')
    end if

50 format('   == Axis A =====================================')
60 format('   ** Moment of gyration : ',A)
65 format('   ** Axis head          : ',A)
80 format('   == Axis B =====================================')
90 format('   ** Moment of gyration : ',A)
95 format('   ** Axis head          : ',A)

end subroutine load_axang

!===============================================================================
! Subroutine:  calculate_axang
!===============================================================================

subroutine calculate_axang(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeAXANG)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,mi,mj,ai,m,info
    real(PMFDP)    :: amass,arg,sc,s1,s2
    real(PMFDP)    :: t11,t22,t33,t12,t13,t23,work(26*3)
    real(PMFDP)    :: totmass1,itotmass1,totmass2,itotmass2
    real(PMFDP)    :: com1(3),com2(3)
    real(PMFDP)    :: t1(3,3),eigenvalues1(3),t2(3,3),eigenvalues2(3)
    real(PMFDP)    :: v(3,3),api(3,3),cij(3),a_xij(3,3,3),b_xij(3,3,3),bint(3,3)
    ! --------------------------------------------------------------------------

    ! calculate centres of masses --------------------
    totmass1 = 0.0d0
    com1(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        com1(:) = com1(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_axang!')
    end if
    itotmass1 = 1.0d0 / totmass1
    com1(:) = com1(:) * itotmass1

    totmass2 = 0.0d0
    com2(:) = 0.0
    do m = cv_item%grps(1)+1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        com2(:) = com2(:) + x(:,ai)*amass
        totmass2 = totmass2 + amass
    end do
    if( totmass2 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_axang!')
    end if
    itotmass2 = 1.0d0 / totmass2
    com2(:) = com2(:) * itotmass2

    ! calculate tensor of gyration
    t11 = 0.0d0
    t22 = 0.0d0
    t33 = 0.0d0
    t12 = 0.0d0
    t13 = 0.0d0
    t23 = 0.0d0

    ! construct matrix
    do m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        t11 = t11 + amass*(x(1,ai) - com1(1))*(x(1,ai) - com1(1))
        t22 = t22 + amass*(x(2,ai) - com1(2))*(x(2,ai) - com1(2))
        t33 = t33 + amass*(x(3,ai) - com1(3))*(x(3,ai) - com1(3))
        t12 = t12 + amass*(x(1,ai) - com1(1))*(x(2,ai) - com1(2))
        t13 = t13 + amass*(x(1,ai) - com1(1))*(x(3,ai) - com1(3))
        t23 = t23 + amass*(x(2,ai) - com1(2))*(x(3,ai) - com1(3))
    end do

    t1(1,1) = t11 * itotmass1
    t1(2,2) = t22 * itotmass1
    t1(3,3) = t33 * itotmass1

    t1(1,2) = t12 * itotmass1
    t1(1,3) = t13 * itotmass1
    t1(2,3) = t23 * itotmass1

    t1(2,1) = t12 * itotmass1
    t1(3,1) = t13 * itotmass1
    t1(3,2) = t23 * itotmass1

    ! calculate eignevalues and eigenvectors of matrix
    eigenvalues1(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 3, t1, 3, eigenvalues1, work, 26*3, info)

    if( info .ne. 0 ) then
    call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize T1 matrix in calculate_axang!')
    end if

    ! calculate tensor of gyration
    t11 = 0.0d0
    t22 = 0.0d0
    t33 = 0.0d0
    t12 = 0.0d0
    t13 = 0.0d0
    t23 = 0.0d0

    ! construct matrix
    do m = cv_item%grps(1) + 1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        t11 = t11 + amass*(x(1,ai) - com2(1))*(x(1,ai) - com2(1))
        t22 = t22 + amass*(x(2,ai) - com2(2))*(x(2,ai) - com2(2))
        t33 = t33 + amass*(x(3,ai) - com2(3))*(x(3,ai) - com2(3))
        t12 = t12 + amass*(x(1,ai) - com2(1))*(x(2,ai) - com2(2))
        t13 = t13 + amass*(x(1,ai) - com2(1))*(x(3,ai) - com2(3))
        t23 = t23 + amass*(x(2,ai) - com2(2))*(x(3,ai) - com2(3))
    end do

    t2(1,1) = t11 * itotmass2
    t2(2,2) = t22 * itotmass2
    t2(3,3) = t33 * itotmass2

    t2(1,2) = t12 * itotmass2
    t2(1,3) = t13 * itotmass2
    t2(2,3) = t23 * itotmass2

    t2(2,1) = t12 * itotmass2
    t2(3,1) = t13 * itotmass2
    t2(3,2) = t23 * itotmass2

    ! calculate eignevalues and eigenvectors of matrix
    eigenvalues2(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 3, t2, 3, eigenvalues2, work, 26*3, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize T2 matrix in calculate_axang!')
    end if

    ! disable axis flipping using axis head atoms
    s1 = sign(1.0d0,dot_product(x(:,cv_item%lindexes(cv_item%head_a))-com1(:),t1(:,cv_item%direction_a)))
    s2 = sign(1.0d0,dot_product(x(:,cv_item%lindexes(cv_item%head_b))-com2(:),t2(:,cv_item%direction_b)))

    arg = dot_product(s1*t1(:,cv_item%direction_a),s2*t2(:,cv_item%direction_b))

    ! write(*,'(F10.4,1X,F4.1,1X,F4.1,1X,6(F10.6,1X),1X,6(F10.6,1X))') arg, s1, s2, com1, &
    !        com1+10.0d0*s1*t1(:,cv_item%direction_a),com2, com2+10.0d0*t2(:,cv_item%direction_b)

    ! coordinate ctx%CVsValues(cv_item%idx) ------------------------------------------------------------

    if ( arg .gt.  1.0 ) then
        arg =  1.0
    else if ( arg .lt. -1.0 ) then
        arg = -1.0
    end if

    ctx%CVsValues(cv_item%idx) = acos(arg)
    sc = sin(ctx%CVsValues(cv_item%idx))

    ! derivatives -----------------------------------------------------------------
    if( abs(sc) .lt. 1.e-12 ) then
        ! avoid division by zero
        sc = -1.e12
    else
        sc = -1.0d0 / sc
    end if

    sc = sc*s1*s2

    ! construct pseudoinverse matrix of A, api
    v(:,:) = t1(:,:)
    api(:,:) = 0.0d0
    do i=1,3
        if( i .ne. cv_item%direction_a ) api(i,i) = 1.0d0/(eigenvalues1(i) - eigenvalues1(cv_item%direction_a))
    end do
    call dgemm('N','N',3,3,3,1.0d0,v,3,api,3,0.0d0,bint,3)
    call dgemm('N','T',3,3,3,1.0d0,bint,3,v,3,0.0d0,api,3)

    ! and solve system of equations
    a_xij(:,:,:) = 0.0d0
    do mi=1,3
        do mj=1,3
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + t1(mj,cv_item%direction_a)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',3,3,-1.0d0,api,3,cij,1,0.0d0,a_xij(:,mi,mj),1)

            a_xij(:,mi,mj) = a_xij(:,mi,mj)*t2(:,cv_item%direction_b)
        end do
    end do

    ! construct pseudoinverse matrix of B, api
    v(:,:) = t2(:,:)
    api(:,:) = 0.0d0
    do i=1,3
        if( i .ne. cv_item%direction_b ) api(i,i) = 1.0d0/(eigenvalues2(i) - eigenvalues2(cv_item%direction_b))
    end do
    call dgemm('N','N',3,3,3,1.0d0,v,3,api,3,0.0d0,bint,3)
    call dgemm('N','T',3,3,3,1.0d0,bint,3,v,3,0.0d0,api,3)

    ! and solve system of equations
    b_xij(:,:,:) = 0.0d0
    do mi=1,3
        do mj=1,3
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + t2(mj,cv_item%direction_b)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',3,3,-1.0d0,api,3,cij,1,0.0d0,b_xij(:,mi,mj),1)

            b_xij(:,mi,mj) = b_xij(:,mi,mj)*t1(:,cv_item%direction_a)
        end do
    end do

    ! and finaly gradients --------------------------
    do m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + sc*amass*itotmass1*( &
                              2.0d0*(x(1,ai) - com1(1))*(a_xij(1,1,1) + a_xij(2,1,1) + a_xij(3,1,1)) &
                            +       (x(2,ai) - com1(2))*(a_xij(1,1,2) + a_xij(2,1,2) + a_xij(3,1,2)) &
                            +       (x(2,ai) - com1(2))*(a_xij(1,2,1) + a_xij(2,2,1) + a_xij(3,2,1)) &
                            +       (x(3,ai) - com1(3))*(a_xij(1,1,3) + a_xij(2,1,3) + a_xij(3,1,3)) &
                            +       (x(3,ai) - com1(3))*(a_xij(1,3,1) + a_xij(2,3,1) + a_xij(3,3,1)))

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + sc*amass*itotmass1*( &
                                    (x(1,ai) - com1(1))*(a_xij(1,1,2) + a_xij(2,1,2) + a_xij(3,1,2)) &
                            +       (x(1,ai) - com1(1))*(a_xij(1,2,1) + a_xij(2,2,1) + a_xij(3,2,1)) &
                            + 2.0d0*(x(2,ai) - com1(2))*(a_xij(1,2,2) + a_xij(2,2,2) + a_xij(3,2,2)) &
                            +       (x(3,ai) - com1(3))*(a_xij(1,2,3) + a_xij(2,2,3) + a_xij(3,2,3)) &
                            +       (x(3,ai) - com1(3))*(a_xij(1,3,2) + a_xij(2,3,2) + a_xij(3,3,2)))

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + sc*amass*itotmass1*( &
                                    (x(1,ai) - com1(1))*(a_xij(1,1,3) + a_xij(2,1,3) + a_xij(3,1,3)) &
                            +       (x(1,ai) - com1(1))*(a_xij(1,3,1) + a_xij(2,3,1) + a_xij(3,3,1)) &
                            +       (x(2,ai) - com1(2))*(a_xij(1,2,3) + a_xij(2,2,3) + a_xij(3,2,3)) &
                            +       (x(2,ai) - com1(2))*(a_xij(1,3,2) + a_xij(2,3,2) + a_xij(3,3,2)) &
                            + 2.0d0*(x(3,ai) - com1(3))*(a_xij(1,3,3) + a_xij(2,3,3) + a_xij(3,3,3)))
    end do

    do m = cv_item%grps(1) + 1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + sc*amass*itotmass2*( &
                              2.0d0*(x(1,ai) - com2(1))*(b_xij(1,1,1) + b_xij(2,1,1) + b_xij(3,1,1)) &
                            +       (x(2,ai) - com2(2))*(b_xij(1,1,2) + b_xij(2,1,2) + b_xij(3,1,2)) &
                            +       (x(2,ai) - com2(2))*(b_xij(1,2,1) + b_xij(2,2,1) + b_xij(3,2,1)) &
                            +       (x(3,ai) - com2(3))*(b_xij(1,1,3) + b_xij(2,1,3) + b_xij(3,1,3)) &
                            +       (x(3,ai) - com2(3))*(b_xij(1,3,1) + b_xij(2,3,1) + b_xij(3,3,1)))

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + sc*amass*itotmass2*( &
                                    (x(1,ai) - com2(1))*(b_xij(1,1,2) + b_xij(2,1,2) + b_xij(3,1,2)) &
                            +       (x(1,ai) - com2(1))*(b_xij(1,2,1) + b_xij(2,2,1) + b_xij(3,2,1)) &
                            + 2.0d0*(x(2,ai) - com2(2))*(b_xij(1,2,2) + b_xij(2,2,2) + b_xij(3,2,2)) &
                            +       (x(3,ai) - com2(3))*(b_xij(1,2,3) + b_xij(2,2,3) + b_xij(3,2,3)) &
                            +       (x(3,ai) - com2(3))*(b_xij(1,3,2) + b_xij(2,3,2) + b_xij(3,3,2)))

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + sc*amass*itotmass2*( &
                                    (x(1,ai) - com2(1))*(b_xij(1,1,3) + b_xij(2,1,3) + b_xij(3,1,3)) &
                            +       (x(1,ai) - com2(1))*(b_xij(1,3,1) + b_xij(2,3,1) + b_xij(3,3,1)) &
                            +       (x(2,ai) - com2(2))*(b_xij(1,2,3) + b_xij(2,2,3) + b_xij(3,2,3)) &
                            +       (x(2,ai) - com2(2))*(b_xij(1,3,2) + b_xij(2,3,2) + b_xij(3,3,2)) &
                            + 2.0d0*(x(3,ai) - com2(3))*(b_xij(1,3,3) + b_xij(2,3,3) + b_xij(3,3,3)))
    end do

    return

end subroutine calculate_axang

!===============================================================================

end module cv_axang

