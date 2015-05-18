!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@enzim.hu
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

module cv_pvang

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePVANG

    integer :: x_direction
    integer :: y_direction

    contains
        procedure :: load_cv        => load_pvang
        procedure :: calculate_cv   => calculate_pvang
end type CVTypePVANG

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_pvang
!===============================================================================

subroutine load_pvang(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common

    implicit none
    class(CVTypePVANG)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: m
    logical                             :: found
    character(len=PRMFILE_MAX_LINE)     :: mask
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'PVANG'
    cv_item%unit          = AngleUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 3
    call cv_common_init_groups(cv_item,prm_fin)

    ! read group a ----------------------------------
    write(PMF_OUT,60)
    call cv_common_read_group(cv_item,prm_fin,1)

    ! load x_direction_a ----------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'x_direction',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'x_direction atom is not specified!')
    end if
    write(PMF_OUT,80) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%x_direction)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=1,cv_item%grps(1)
        if( cv_item%rindexes(m) .eq. cv_item%x_direction ) then
            cv_item%x_direction = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'x_direction atom has to be part of group_a!')
    end if

    ! load y_direction_a ----------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'y_direction',mask) ) then
        call pmf_utils_exit(PMF_OUT,1,'y_direction atom is not specified!')
    end if
    write(PMF_OUT,80) trim(mask)
    call cv_common_set_atom_from_mask(mask,cv_item%y_direction)

    ! find atom in cv_item%rindexes
    found = .false.
    do m=1,cv_item%grps(1)
        if( cv_item%rindexes(m) .eq. cv_item%y_direction ) then
            cv_item%y_direction = m
            found = .true.
            exit
        end if
    end do

    if( .not. found ) then
        call pmf_utils_exit(PMF_OUT,1,'y_direction atom has to be part of group_a!')
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,85)
    call cv_common_read_group(cv_item,prm_fin,2)
    call cv_common_read_group(cv_item,prm_fin,3)


    60 format('   == Plane A ====================================')
    70 format('   ** x-direction atom   : ',A)
    80 format('   ** y-direction atom   : ',A)
    85 format('   == Vector B ===================================')

end subroutine load_pvang

!===============================================================================
! Subroutine:  calculate_pvang
!===============================================================================

subroutine calculate_pvang(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypePVANG)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,ai,m,info,orient,mi,mj
    real(PMFDP)    :: d1(3),d2(3),d3(3),dx(3),dzx(3),dzy(3),dzz(3),cang
    real(PMFDP)    :: a(3,3),a11,a22,a33,a12,a13,a23,st
    real(PMFDP)    :: eigenvalues(3)
    real(PMFDP)    :: totmass1,totmass2,totmass3,amass,ac,dzz_s,dzz_s2
    real(PMFDP)    :: work(26*3),sc
    real(PMFDP)    :: v(3,3),api(3,3),cij(3),xij(3,3,3),bint(3,3)
    ! --------------------------------------------------------------------------

    orient = 1   ! best plane orientation for group points A

    ! calculate centres of mases --------------------
    totmass1 = 0.0d0
    d1(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d1(:) = d1(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_pvang!')
    end if
    d1(:) = d1(:) / totmass1

    ! calculate parameters of plane A -----------------
    a11 = 0.0d0
    a22 = 0.0d0
    a33 = 0.0d0
    a12 = 0.0d0
    a13 = 0.0d0
    a23 = 0.0d0

    ! construct matrix
    do m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        a11 = a11 + amass*(x(1,ai) - d1(1))**2
        a22 = a22 + amass*(x(2,ai) - d1(2))**2
        a33 = a33 + amass*(x(3,ai) - d1(3))**2
        a12 = a12 + amass*(x(1,ai) - d1(1))*(x(2,ai) - d1(2))
        a13 = a13 + amass*(x(1,ai) - d1(1))*(x(3,ai) - d1(3))
        a23 = a23 + amass*(x(2,ai) - d1(2))*(x(3,ai) - d1(3))
    end do

    a(1,1) = a11
    a(2,2) = a22
    a(3,3) = a33

    a(1,2) = a12
    a(1,3) = a13
    a(2,3) = a23

    a(2,1) = a12
    a(3,1) = a13
    a(3,2) = a23

    ! calculate eignevalues and eigenvectors of matrix
    eigenvalues(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 3, a, 3, eigenvalues, work, 26*3, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_pvang (plane A)!')
    end if

    ! determine z-axis orientation
    dzx(:) = x(:,cv_item%lindexes(cv_item%x_direction)) - d1(:)
    dzy(:) = x(:,cv_item%lindexes(cv_item%y_direction)) - d1(:)

    ! cross product
    dzz(1) = dzx(2)*dzy(3) - dzx(3)*dzy(2)
    dzz(2) = dzx(3)*dzy(1) - dzx(1)*dzy(3)
    dzz(3) = dzx(1)*dzy(2) - dzx(2)*dzy(1)

    ! norm dzz (this is not neccessary, but it is keep for debug)
    dzz_s = sqrt(dzz(1)**2 + dzz(2)**2 + dzz(3)**2)
    dzz(:) = dzz(:)/dzz_s

    ! angle between z-axis and plane vector
    ac = dzz(1)*a(1,orient) + dzz(2)*a(2,orient) + dzz(3)*a(3,orient)

    ! correct vector direction
    a(:,orient) = a(:,orient)*sign(1.0d0,ac)

    ! calculate parameters of vector B -----------------

    totmass2 = 0.0d0
    d2(:) = 0.0d0
    do  m = cv_item%grps(1) + 1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d2(:) = d2(:) + x(:,ai)*amass
        totmass2 = totmass2 + amass
    end do
    if( totmass2 .gt. 0 ) d2(:) = d2(:) / totmass2

    totmass3 = 0.0d0
    d3(:) = 0.0d0
    do  m = cv_item%grps(2) + 1, cv_item%grps(3)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d3(:) = d3(:) + x(:,ai)*amass
        totmass3 = totmass3 + amass
    end do
    if( totmass3 .gt. 0 ) d3(:) = d3(:) / totmass3

    ! vector
    dx(:) = d3(:) - d2(:)

    ! normalize vector
    dzz_s2 = dx(1)**2 + dx(2)**2 + dx(3)**2
    dzz_s = sqrt(dzz_s2)
    if( dzz_s .gt. 0d0 ) dx(:) = dx(:) / dzz_s

    ! calculate the angle between plane and vector -----------------------------

    cang = a(1,orient)*dx(1) + a(2,orient)*dx(2) + a(3,orient)*dx(3)

    if ( cang .gt.  1.0 ) then
        cang =  1.0
    else if ( cang .lt. -1.0 ) then
        cang = -1.0
    end if

    ctx%CVsValues(cv_item%idx) = acos(cang)

    ! derivatives --------------------------------------------------------------

    sc = sin(ctx%CVsValues(cv_item%idx))
    if( abs(sc) .lt. 1.e-12 ) then
        ! avoid division by zero
        sc = -1.e12
    else
        sc = -1.0d0 / sc
    end if

    ! eigenvector derivatives -----------------------

    ! construct pseudoinverse matrix of A, api
    v(:,:) = a(:,:)
    api(:,:) = 0.0d0
    do i=1,3
        if( i .ne. orient ) api(i,i) = 1.0d0/(eigenvalues(i) - eigenvalues(orient))
    end do
    call dgemm('N','N',3,3,3,1.0d0,v,3,api,3,0.0d0,bint,3)
    call dgemm('N','T',3,3,3,1.0d0,bint,3,v,3,0.0d0,api,3)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,3
        do mj=1,3
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + a(mj,orient)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',3,3,-1.0d0,api,3,cij,1,0.0d0,xij(:,mi,mj),1)

            xij(:,mi,mj) = sc*xij(:,mi,mj)*dx
        end do
    end do

    ! and finaly gradients --------------------------
    do m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + amass*( &
                              2.0d0*(x(1,ai) - d1(1))*(xij(1,1,1) + xij(2,1,1) + xij(3,1,1)) &
                            +       (x(2,ai) - d1(2))*(xij(1,1,2) + xij(2,1,2) + xij(3,1,2)) &
                            +       (x(2,ai) - d1(2))*(xij(1,2,1) + xij(2,2,1) + xij(3,2,1)) &
                            +       (x(3,ai) - d1(3))*(xij(1,1,3) + xij(2,1,3) + xij(3,1,3)) &
                            +       (x(3,ai) - d1(3))*(xij(1,3,1) + xij(2,3,1) + xij(3,3,1)))

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + amass*( &
                                    (x(1,ai) - d1(1))*(xij(1,1,2) + xij(2,1,2) + xij(3,1,2)) &
                            +       (x(1,ai) - d1(1))*(xij(1,2,1) + xij(2,2,1) + xij(3,2,1)) &
                            + 2.0d0*(x(2,ai) - d1(2))*(xij(1,2,2) + xij(2,2,2) + xij(3,2,2)) &
                            +       (x(3,ai) - d1(3))*(xij(1,2,3) + xij(2,2,3) + xij(3,2,3)) &
                            +       (x(3,ai) - d1(3))*(xij(1,3,2) + xij(2,3,2) + xij(3,3,2)))

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + amass*( &
                                    (x(1,ai) - d1(1))*(xij(1,1,3) + xij(2,1,3) + xij(3,1,3)) &
                            +       (x(1,ai) - d1(1))*(xij(1,3,1) + xij(2,3,1) + xij(3,3,1)) &
                            +       (x(2,ai) - d1(2))*(xij(1,2,3) + xij(2,2,3) + xij(3,2,3)) &
                            +       (x(2,ai) - d1(2))*(xij(1,3,2) + xij(2,3,2) + xij(3,3,2)) &
                            + 2.0d0*(x(3,ai) - d1(3))*(xij(1,3,3) + xij(2,3,3) + xij(3,3,3)))
    end do

    ! vector derivatives -----------------------

    sc = sc*dzz_s/dzz_s2

    do m = cv_item%grps(1) + 1, cv_item%grps(2)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        st = sc*amass/totmass2
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) - st*(a(1,orient) - dx(1)*cang)
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) - st*(a(2,orient) - dx(2)*cang)
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) - st*(a(3,orient) - dx(3)*cang)
    end do

    do m = cv_item%grps(2) + 1, cv_item%grps(3)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        st = sc*amass/totmass3
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + st*(a(1,orient) - dx(1)*cang)
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + st*(a(2,orient) - dx(2)*cang)
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + st*(a(3,orient) - dx(3)*cang)
    end do

    return

end subroutine calculate_pvang

!===============================================================================

end module cv_pvang

