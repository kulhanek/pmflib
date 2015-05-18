!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_evec

use pmf_sizes
use pmf_constants
use pmf_dat
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeEVEC

    type(XYZFILE_TYPE)          :: xyz_str      ! reference structure
    double precision, pointer   :: evec(:,:)    ! essential vector

    contains
        procedure :: load_cv        => load_evec
        procedure :: calculate_cv   => calculate_evec
end type CVTypeEVEC

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_evec
!===============================================================================

subroutine load_evec(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use smf_xyzfile_type
    use smf_xyzfile

    implicit none
    class(CVTypeEVEC)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_VALUE)    :: file_name
    integer                             :: alloc_failed, i, nvsize
    character(len=PRMFILE_MAX_VALUE)    :: evec_comment
    ! -----------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'EVEC'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_read_groups(cv_item,prm_fin)

    ! read reference structure -------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'average',file_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of average structure (average) is not specified!')
    end if
    write(PMF_OUT,50) trim(file_name)

    call init_xyz(cv_item%xyz_str)
    call open_xyz(PMF_XYZ,file_name,cv_item%xyz_str,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str)
    write(PMF_OUT,70) trim(cv_item%xyz_str%comment)
    call close_xyz(PMF_XYZ,cv_item%xyz_str)

    if( cv_item%xyz_str%natoms .ne. cv_item%natoms ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group and target structure differ!')
    end if

    ! allocate evec array
    allocate(cv_item%evec(3,cv_item%natoms),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    if( .not. prmfile_get_string_by_key(prm_fin,'evec',file_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of essential vector (evec) is not specified!')
    end if
    write(PMF_OUT,60) trim(file_name)

    open(unit=PMF_EVEC, file=file_name, status='OLD', form='FORMATTED', err=100)

    ! header
    read(PMF_EVEC,*,err=105,end=105) nvsize
    if( nvsize .ne. 3*cv_item%natoms ) then
        call pmf_utils_exit(PMF_OUT,1,'Essential vector does not have correct size!')
    end if

    read(PMF_EVEC,*,err=105,end=105) evec_comment
    write(PMF_OUT,70) trim(evec_comment)

    ! data
    do i=1, cv_item%natoms
        read(PMF_EVEC,*,err=110,end=110) cv_item%evec(1,i), cv_item%evec(2,i), cv_item%evec(3,i)
    end do
    close(PMF_EVEC)

    return

100 call pmf_utils_exit(PMF_OUT,1,'Unable to open essential vector file!')
105 call pmf_utils_exit(PMF_OUT,1,'Unable to read header of essential vector file!')
110 call pmf_utils_exit(PMF_OUT,1,'Unable to read essential vector file!')

50 format('   ** average structure  : ',A)
60 format('   ** essential vector   : ',A)
70 format('      >> ',A)

end subroutine load_evec

!===============================================================================
! Subroutine:  calculate_evec
!===============================================================================

subroutine calculate_evec(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeEVEC)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,j,k,ai,aj,info,best,mi,mj
    real(PMFDP)    :: x1,x2,x3,xr1,xr2,xr3,amass,totmass,itotmass,value,sc
    real(PMFDP)    :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)    :: f(4,4),u(3,3),ut(3,3),duq(4),dqx(4),dqy(4),dqz(4)
    real(PMFDP)    :: eigenvalues(4),work(26*4)
    real(PMFDP)    :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4)
    ! -----------------------------------------------------------------------------

    ! calculate gemetrical centres (reference and source) -------------------
    x1 = 0.0d0
    x2 = 0.0d0
    x3 = 0.0d0
    xr1 = 0.0d0
    xr2 = 0.0d0
    xr3 = 0.0d0
    totmass = 0.0d0

    do  i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        ! source
        x1 = x1 + x(1,ai)*amass
        x2 = x2 + x(2,ai)*amass
        x3 = x3 + x(3,ai)*amass

        ! reference
        xr1 = xr1 + cv_item%xyz_str%cvs(1,i)*amass
        xr2 = xr2 + cv_item%xyz_str%cvs(2,i)*amass
        xr3 = xr3 + cv_item%xyz_str%cvs(3,i)*amass

        totmass = totmass + amass
    end do

    itotmass = 1.0d0 / totmass
    x1 = x1 * itotmass
    x2 = x2 * itotmass
    x3 = x3 * itotmass
    xr1 = xr1 * itotmass
    xr2 = xr2 * itotmass
    xr3 = xr3 * itotmass

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

    do i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)

        r11 = r11 + amass*(x(1,ai) - x1)*(cv_item%xyz_str%cvs(1,i) - xr1)
        r12 = r12 + amass*(x(1,ai) - x1)*(cv_item%xyz_str%cvs(2,i) - xr2)
        r13 = r13 + amass*(x(1,ai) - x1)*(cv_item%xyz_str%cvs(3,i) - xr3)

        r21 = r21 + amass*(x(2,ai) - x2)*(cv_item%xyz_str%cvs(1,i) - xr1)
        r22 = r22 + amass*(x(2,ai) - x2)*(cv_item%xyz_str%cvs(2,i) - xr2)
        r23 = r23 + amass*(x(2,ai) - x2)*(cv_item%xyz_str%cvs(3,i) - xr3)

        r31 = r31 + amass*(x(3,ai) - x3)*(cv_item%xyz_str%cvs(1,i) - xr1)
        r32 = r32 + amass*(x(3,ai) - x3)*(cv_item%xyz_str%cvs(2,i) - xr2)
        r33 = r33 + amass*(x(3,ai) - x3)*(cv_item%xyz_str%cvs(3,i) - xr3)
    end do

    ! construct matrix for quaterion fitting
    f(1,1) =  r11 + r22 + r33
    f(1,2) =  r23 - r32
    f(1,3) =  r31 - r13
    f(1,4) =  r12 - r21

    f(2,1) =  r23 - r32
    f(2,2) =  r11 - r22 - r33
    f(2,3) =  r12 + r21
    f(2,4) =  r13 + r31

    f(3,1) =  r31 - r13
    f(3,2) =  r12 + r21
    f(3,3) = -r11 + r22 - r33
    f(3,4) =  r23 + r32

    f(4,1) =  r12 - r21
    f(4,2) =  r13 + r31
    f(4,3) =  r23 + r32
    f(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvalues(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, f, 4, eigenvalues, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_evec!')
    end if

    best = 4

    ! rotation matrix ------------------------------
    u(1,1) = f(1,best)**2 + f(2,best)**2 - f(3,best)**2 - f(4,best)**2
    u(1,2) = 2.0d0*( f(2,best)*f(3,best) - f(1,best)*f(4,best) )
    u(1,3) = 2.0d0*( f(2,best)*f(4,best) + f(1,best)*f(3,best) )

    u(2,1) = 2.0d0*( f(2,best)*f(3,best) + f(1,best)*f(4,best) )
    u(2,2) = f(1,best)**2 - f(2,best)**2 + f(3,best)**2 - f(4,best)**2
    u(2,3) = 2.0d0*( f(3,best)*f(4,best) - f(1,best)*f(2,best) )

    u(3,1) = 2.0d0*( f(2,best)*f(4,best) - f(1,best)*f(3,best) )
    u(3,2) = 2.0d0*( f(3,best)*f(4,best) + f(1,best)*f(2,best) )
    u(3,3) = f(1,best)**2 - f(2,best)**2 - f(3,best)**2 + f(4,best)**2

    ! transpose u matrix ------
    ut(1,1) = u(1,1)
    ut(1,2) = u(2,1)
    ut(1,3) = u(3,1)

    ut(2,1) = u(1,2)
    ut(2,2) = u(2,2)
    ut(2,3) = u(3,2)

    ut(3,1) = u(1,3)
    ut(3,2) = u(2,3)
    ut(3,3) = u(3,3)

    value = 0.0d0

    ! calculate projections ----
    do  i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)

        r11 = x(1,ai)  - x1
        r12 = x(2,ai)  - x2
        r13 = x(3,ai)  - x3

        r21 = cv_item%xyz_str%cvs(1,i) - xr1
        r22 = cv_item%xyz_str%cvs(2,i) - xr2
        r23 = cv_item%xyz_str%cvs(3,i) - xr3

        ! transform r1 point
        r31 = u(1,1)*r11 + u(1,2)*r12 + u(1,3)*r13
        r32 = u(2,1)*r11 + u(2,2)*r12 + u(2,3)*r13
        r33 = u(3,1)*r11 + u(3,2)*r12 + u(3,3)*r13

        value = value + (r31 - r21)*cv_item%evec(1,i) + (r32 - r22)*cv_item%evec(2,i) + (r33 - r23)*cv_item%evec(3,i)
    end do

    ctx%CVsValues(cv_item%idx) = value

    ! construct pseudoinverse matrix of A, api
    v(:,:) = f(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvalues(i) - eigenvalues(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + f(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

    ! derivative of
    ! uxe

    ! calculate derivatives ----
    do  i = 1, cv_item%natoms

        ! transform evec
        r31 = ut(1,1)*cv_item%evec(1,i) + ut(1,2)*cv_item%evec(2,i) + ut(1,3)*cv_item%evec(3,i)
        r32 = ut(2,1)*cv_item%evec(1,i) + ut(2,2)*cv_item%evec(2,i) + ut(2,3)*cv_item%evec(3,i)
        r33 = ut(3,1)*cv_item%evec(1,i) + ut(3,2)*cv_item%evec(2,i) + ut(3,3)*cv_item%evec(3,i)

        ! get qx derivatives
        r11 = cv_item%xyz_str%cvs(1,i) - xr1
        r12 = cv_item%xyz_str%cvs(2,i) - xr2
        r13 = cv_item%xyz_str%cvs(3,i) - xr3

        do k=1,4
            dqx(k) =   xij(k,1,1)*r11 - xij(k,1,3)*r13 + xij(k,1,4)*r12 &
                     + xij(k,2,2)*r11 + xij(k,2,3)*r12 + xij(k,2,4)*r13 &
                     - xij(k,3,1)*r13 + xij(k,3,2)*r12 - xij(k,3,3)*r11 &
                     + xij(k,4,1)*r12 + xij(k,4,2)*r13 - xij(k,4,4)*r11

            dqy(k) =   xij(k,1,1)*r12 + xij(k,1,2)*r13 - xij(k,1,4)*r11 &
                     + xij(k,2,1)*r13 - xij(k,2,2)*r12 + xij(k,2,3)*r11 &
                     + xij(k,3,2)*r11 + xij(k,3,3)*r12 + xij(k,3,4)*r13 &
                     - xij(k,4,1)*r11 + xij(k,4,3)*r13 - xij(k,4,4)*r12

            dqz(k) =   xij(k,1,1)*r13 - xij(k,1,2)*r12 + xij(k,1,3)*r11 &
                     - xij(k,2,1)*r12 - xij(k,2,2)*r13 + xij(k,2,4)*r11 &
                     + xij(k,3,1)*r11 - xij(k,3,3)*r13 + xij(k,3,4)*r12 &
                     + xij(k,4,2)*r11 + xij(k,4,3)*r12 + xij(k,4,4)*r13
        end do

        ! first part
        do  j = 1, cv_item%natoms
            aj = cv_item%lindexes(j)
            sc = - mass(aj) * itotmass
            if( i .eq. j ) then
                sc = 1.0d0 + sc
            end if
            ctx%CVsDrvs(1,aj,cv_item%idx) = ctx%CVsDrvs(1,aj,cv_item%idx) + sc*r31
            ctx%CVsDrvs(2,aj,cv_item%idx) = ctx%CVsDrvs(2,aj,cv_item%idx) + sc*r32
            ctx%CVsDrvs(3,aj,cv_item%idx) = ctx%CVsDrvs(3,aj,cv_item%idx) + sc*r33
        end do

        ! second part

!    u11(q0,q1,q2,q3) := q0*q0+q1*q1-q2*q2-q3*q3;
!    u12(q0,q1,q2,q3) := 2*(q1*q2-q0*q3);
!    u13(q0,q1,q2,q3) := 2*(q1*q3+q0*q2);
!    u21(q0,q1,q2,q3) := 2*(q1*q2+q0*q3);
!    u22(q0,q1,q2,q3) := q0*q0 -q1*q1 + q2*q2 -q3*q3;
!    u23(q0,q1,q2,q3) := 2*(q2*q3-q0*q1);
!    u31(q0,q1,q2,q3) := 2*(q1*q3-q0*q2);
!    u32(q0,q1,q2,q3) := 2*(q2*q3+q0*q1);
!    u33(q0,q1,q2,q3) := q0*q0 -q1*q1 - q2*q2 +q3*q3;

        ai = cv_item%lindexes(i)
        amass = mass(ai)

        do  j = 1, cv_item%natoms
            aj = cv_item%lindexes(j)
            r21 = x(1,aj) - x1
            r22 = x(2,aj) - x2
            r23 = x(3,aj) - x3
!           (2*e1*q2-2*e2*q1+2*e3*q0)*x3+(-2*e1*q3+2*e3*q1+2*e2*q0)*x2+(2*e2*q3-2*e3*q2+2*e1*q0)*x1
            duq(1) = ( 2*cv_item%evec(1,j)*f(3,best)-2*cv_item%evec(2,j)*f(2,best)+2*cv_item%evec(3,j)*f(1,best))*r23 &
                   + (-2*cv_item%evec(1,j)*f(4,best)+2*cv_item%evec(3,j)*f(2,best)+2*cv_item%evec(2,j)*f(1,best))*r22 &
                   + ( 2*cv_item%evec(2,j)*f(4,best)-2*cv_item%evec(3,j)*f(3,best)+2*cv_item%evec(1,j)*f(1,best))*r21

!           (2*e1*q3-2*e3*q1-2*e2*q0)*x3+(2*e1*q2-2*e2*q1+2*e3*q0)*x2+(2*e3*q3+2*e2*q2+2*e1*q1)*x1
            duq(2) = (2*cv_item%evec(1,j)*f(4,best)-2*cv_item%evec(3,j)*f(2,best)-2*cv_item%evec(2,j)*f(1,best))*r23 &
                   + (2*cv_item%evec(1,j)*f(3,best)-2*cv_item%evec(2,j)*f(2,best)+2*cv_item%evec(3,j)*f(1,best))*r22 &
                   + (2*cv_item%evec(3,j)*f(4,best)+2*cv_item%evec(2,j)*f(3,best)+2*cv_item%evec(1,j)*f(2,best))*r21

!           (2*e2*q3-2*e3*q2+2*e1*q0)*x3+(2*e3*q3+2*e2*q2+2*e1*q1)*x2+(-2*e1*q2+2*e2*q1-2*e3*q0)*x1
            duq(3) = ( 2*cv_item%evec(2,j)*f(4,best)-2*cv_item%evec(3,j)*f(3,best)+2*cv_item%evec(1,j)*f(1,best))*r23 &
                   + ( 2*cv_item%evec(3,j)*f(4,best)+2*cv_item%evec(2,j)*f(3,best)+2*cv_item%evec(1,j)*f(2,best))*r22 &
                   + (-2*cv_item%evec(1,j)*f(3,best)+2*cv_item%evec(2,j)*f(2,best)-2*cv_item%evec(3,j)*f(1,best))*r21

!           (2*e3*q3+2*e2*q2+2*e1*q1)*x3+(-2*e2*q3+2*e3*q2-2*e1*q0)*x2+(-2*e1*q3+2*e3*q1+2*e2*q0)*x1
            duq(4) = ( 2*cv_item%evec(3,j)*f(4,best)+2*cv_item%evec(2,j)*f(3,best)+2*cv_item%evec(1,j)*f(2,best))*r23 &
                   + (-2*cv_item%evec(2,j)*f(4,best)+2*cv_item%evec(3,j)*f(3,best)-2*cv_item%evec(1,j)*f(1,best))*r22 &
                   + (-2*cv_item%evec(1,j)*f(4,best)+2*cv_item%evec(3,j)*f(2,best)+2*cv_item%evec(2,j)*f(1,best))*r21

            ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                                          + amass*(duq(1)*dqx(1)+duq(2)*dqx(2)+duq(3)*dqx(3)+duq(4)*dqx(4))
            ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                                          + amass*(duq(1)*dqy(1)+duq(2)*dqy(2)+duq(3)*dqy(3)+duq(4)*dqy(4))
            ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                                          + amass*(duq(1)*dqz(1)+duq(2)*dqz(2)+duq(3)*dqz(3)+duq(4)*dqz(4))
        end do

    end do

    return

end subroutine calculate_evec

!===============================================================================

end module cv_evec

