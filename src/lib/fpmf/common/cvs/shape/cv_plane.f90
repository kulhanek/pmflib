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

! the root mean square deviation to the best fit plane

module cv_plane

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePLANE
    contains
        procedure :: load_cv        => load_plane
        procedure :: calculate_cv   => calculate_plane
end type CVTypePLANE

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_plane
!===============================================================================

subroutine load_plane(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common

    implicit none
    class(CVTypePLANE)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'PLANE'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_init_groups(cv_item,prm_fin)

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    return

50 format('   == Plane ======================================')

end subroutine load_plane

!===============================================================================
! Subroutine:  calculate_plane
!===============================================================================

subroutine calculate_plane(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypePLANE)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,ai,aj,m,n,info,mi,mj,orient
    real(PMFDP)    :: d1(3),dx(3)
    real(PMFDP)    :: a(3,3),a11,a22,a33,a12,a13,a23,eigenvalues(3)
    real(PMFDP)    :: totmass1,amass,sd,sc
    real(PMFDP)    :: work(26*3)
    real(PMFDP)    :: v(3,3),api(3,3),cij(3),xij(3,3,3),bint(3,3),txij(3,3,3)
    ! -----------------------------------------------------------------------------

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
    if( totmass1 .gt. 0 ) d1(:) = d1(:) / totmass1

    ! calculate parameters of plane -----------------
    a11 = 0.0d0
    a22 = 0.0d0
    a33 = 0.0d0
    a12 = 0.0d0
    a13 = 0.0d0
    a23 = 0.0d0

    ! construct matrix:
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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_plane!')
    end if

    ! root mean square distance -------------------
    ctx%CVsValues(cv_item%idx) = 0.0d0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        dx(:) = x(:,ai) - d1(:)
        sd = a(1,orient)*dx(1) + a(2,orient)*dx(2) + a(3,orient)*dx(3)
        ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + amass * sd**2
    end do
    if( totmass1 .gt. 0 ) ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) / totmass1
    ctx%CVsValues(cv_item%idx) = sqrt(ctx%CVsValues(cv_item%idx))

    ! first derivatives --------------------------
    if( ctx%CVsValues(cv_item%idx) .gt. 1.0d-7 ) then
        sc = 1.0d0 / ctx%CVsValues(cv_item%idx)
    else
        sc = 0.0
    end if

    ! first part, e.g. a*dx'
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        dx(:) = x(:,ai) - d1(:)
        sd = a(1,orient)*dx(1) + a(2,orient)*dx(2) + a(3,orient)*dx(3)
        ctx%CVsDrvs(:,ai,cv_item%idx) =  ctx%CVsDrvs(:,ai,cv_item%idx) + sc*sd*amass*a(:,1)/totmass1
    end do

    ! second part, e.g. a'dx

    ! eigenvector derivatives -----------------------

    ! construct pseudoinverse matrix, api
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
        end do
    end do

    txij = xij

    ! and finaly gradients --------------------------
    do m = 1, cv_item%grps(1)
        cij(:) = 0.0d0
        do n = 1, cv_item%grps(1)
            aj = cv_item%lindexes(n)
            amass = mass(aj)
            dx(:) = x(:,aj) - d1(:)
            sd = a(1,orient)*dx(1) + a(2,orient)*dx(2) + a(3,orient)*dx(3)
            cij(:) = cij(:) + sc*sd*amass*dx(:)/totmass1
        end do
        do mi = 1,3
            xij(mi,:,:) = txij(mi,:,:)*cij(mi)
        end do
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                            + amass*(2.0d0*(x(1,ai) - d1(1))*(xij(1,1,1) + xij(2,1,1) + xij(3,1,1)) &
                            +       (x(2,ai) - d1(2))*(xij(1,1,2) + xij(2,1,2) + xij(3,1,2)) &
                            +       (x(2,ai) - d1(2))*(xij(1,2,1) + xij(2,2,1) + xij(3,2,1)) &
                            +       (x(3,ai) - d1(3))*(xij(1,1,3) + xij(2,1,3) + xij(3,1,3)) &
                            +       (x(3,ai) - d1(3))*(xij(1,3,1) + xij(2,3,1) + xij(3,3,1)))

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                            + amass*((x(1,ai) - d1(1))*(xij(1,1,2) + xij(2,1,2) + xij(3,1,2)) &
                            +       (x(1,ai) - d1(1))*(xij(1,2,1) + xij(2,2,1) + xij(3,2,1)) &
                            + 2.0d0*(x(2,ai) - d1(2))*(xij(1,2,2) + xij(2,2,2) + xij(3,2,2)) &
                            +       (x(3,ai) - d1(3))*(xij(1,2,3) + xij(2,2,3) + xij(3,2,3)) &
                            +       (x(3,ai) - d1(3))*(xij(1,3,2) + xij(2,3,2) + xij(3,3,2)))

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                            + amass*((x(1,ai) - d1(1))*(xij(1,1,3) + xij(2,1,3) + xij(3,1,3)) &
                            +       (x(1,ai) - d1(1))*(xij(1,3,1) + xij(2,3,1) + xij(3,3,1)) &
                            +       (x(2,ai) - d1(2))*(xij(1,2,3) + xij(2,2,3) + xij(3,2,3)) &
                            +       (x(2,ai) - d1(2))*(xij(1,3,2) + xij(2,3,2) + xij(3,3,2)) &
                            + 2.0d0*(x(3,ai) - d1(3))*(xij(1,3,3) + xij(2,3,3) + xij(3,3,3)))
    end do

    return

end subroutine calculate_plane

!===============================================================================

end module cv_plane

