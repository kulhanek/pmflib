!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2021 Ivo Durnik,
!    Copyright (C) 2008 Silvia Cereda, sc578@cam.ac.uk &
!                       Petr Kulhanek, kulhanek@enzim.hu
!
!    This library is free software; you can repostribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is postributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module cv_orad2

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeORAD2

    type(XYZFILE_TYPE)  :: xyz_str_a
    integer :: direction

    contains
        procedure :: load_cv        => load_orad2
        procedure :: calculate_cv   => calculate_orad2
end type CVTypeORAD2

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_orad2
!===============================================================================

subroutine load_orad2(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeORAD2)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(1)                        :: cdir

    integer                             :: i,ar
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    logical                             :: skiptest,lresult
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'OPOS2'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! init groups and control array -----------------
    cv_item%ngrps = 2
    call cv_common_init_groups(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    ! read reference structure ----------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_a',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_a) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_a)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a)

    if( cv_item%xyz_str_a%natoms .ne. cv_item%grps(1) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group A and reference structure A differs!')
    end if

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,90)
    call cv_common_read_group(cv_item,prm_fin,2)

    ! load orientation ------------------------------
    write(PMF_OUT,80)
    if( .not. prmfile_get_string_by_key(prm_fin,'direction',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'direction is not specified!')
    else
        write(PMF_OUT,60) cdir
    end if

    select case(cdir)
        case('x')
            cv_item%direction = 1
        case('y')
            cv_item%direction = 2
        case('z')
            cv_item%direction = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,&
                                'direction has to be x, y, or z!')
    end select

    return

 50 format('   == System A ===================================')
 70 format('   ** reference structure: ',A)
 90 format('   == Point B ====================================')
 60 format('   ** direction          : ',A1)
 80 format('   -----------------------------------------------')
100 format('Atom mismatch between group A and reference A atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

end subroutine load_orad2

!===============================================================================
! Subroutine:  calculate_orad2
!===============================================================================

subroutine calculate_orad2(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeORAD2)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,mi,mj,ai,info,best
    real(PMFDP)         :: px(3)
    real(PMFDP)         :: f(4,4), work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: u(3,3)
    real(PMFDP)         :: eigenvalues(4)
    real(PMFDP)         :: xs(3),xr(3),xb(3)
    real(PMFDP)         :: ingra,ingrb
    real(PMFDP)         :: o(3),tmp1(3),vx(3)
    ! --------------------------------------------------------------------------
    real(PMFDP)         :: sc
    real(PMFDP)         :: a_u(3,3)
    real(PMFDP)         :: a_o(3),a_b(3),a_vx(3),a_px(3)
    real(PMFDP)         :: a_fa(4),a_rij(4,4)
    real(PMFDP)         :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4)
    real(PMFDP)         :: l_u(3,3),a_xs(3)
    ! --------------------------------------------------------------------------

    ! inverse number of atoms
    ingra = 1.0d0 / cv_item%grps(1)
    ingrb = 1.0d0 / (cv_item%grps(2)-cv_item%grps(1))

    ! calculate geometrical centres system A (source and target) -------------------
    xs(:) = 0.0d0
    xr(:) = 0.0d0

    do  i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        ! source
        xs(:) = xs(:) + x(:,ai)

        ! reference
        xr(:) = xr(:) + cv_item%xyz_str_a%cvs(:,i)
    end do

    xs(:) = xs(:) * ingra
    xr(:) = xr(:) * ingra

    ! calculate geometrical center point B -------------------
    xb(:) = 0.0d0

    do  i = cv_item%grps(1)+1,cv_item%grps(2)
        ai = cv_item%lindexes(i)

        xb(:) = xb(:) + x(:,ai)
    end do

    xb(:) = xb(:) * ingrb

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

    do i = 1,cv_item%grps(1)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xs(1))*(cv_item%xyz_str_a%cvs(1,i) - xr(1))
        r12 = r12 + (x(1,ai) - xs(1))*(cv_item%xyz_str_a%cvs(2,i) - xr(2))
        r13 = r13 + (x(1,ai) - xs(1))*(cv_item%xyz_str_a%cvs(3,i) - xr(3))

        r21 = r21 + (x(2,ai) - xs(2))*(cv_item%xyz_str_a%cvs(1,i) - xr(1))
        r22 = r22 + (x(2,ai) - xs(2))*(cv_item%xyz_str_a%cvs(2,i) - xr(2))
        r23 = r23 + (x(2,ai) - xs(2))*(cv_item%xyz_str_a%cvs(3,i) - xr(3))

        r31 = r31 + (x(3,ai) - xs(3))*(cv_item%xyz_str_a%cvs(1,i) - xr(1))
        r32 = r32 + (x(3,ai) - xs(3))*(cv_item%xyz_str_a%cvs(2,i) - xr(2))
        r33 = r33 + (x(3,ai) - xs(3))*(cv_item%xyz_str_a%cvs(3,i) - xr(3))
    end do

    r11 = r11 * ingra
    r12 = r12 * ingra
    r13 = r13 * ingra

    r21 = r21 * ingra
    r22 = r22 * ingra
    r23 = r23 * ingra

    r31 = r31 * ingra
    r32 = r32 * ingra
    r33 = r33 * ingra

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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_orad2!')
    end if

    best = 4

    ! rotation matrix ------------------------------
    u(1,1) = f(1,best)**2 + f(2,best)**2 - f(3,best)**2 - f(4,best)**2
    u(2,1) = 2.0d0*( f(2,best)*f(3,best) - f(1,best)*f(4,best) )
    u(3,1) = 2.0d0*( f(2,best)*f(4,best) + f(1,best)*f(3,best) )

    u(1,2) = 2.0d0*( f(2,best)*f(3,best) + f(1,best)*f(4,best) )
    u(2,2) = f(1,best)**2 - f(2,best)**2 + f(3,best)**2 - f(4,best)**2
    u(3,2) = 2.0d0*( f(3,best)*f(4,best) - f(1,best)*f(2,best) )

    u(1,3) = 2.0d0*( f(2,best)*f(4,best) - f(1,best)*f(3,best) )
    u(2,3) = 2.0d0*( f(3,best)*f(4,best) + f(1,best)*f(2,best) )
    u(3,3) = f(1,best)**2 - f(2,best)**2 - f(3,best)**2 + f(4,best)**2

    ! get origin of group_a
    o(:) = 0.0d0
    ! move reference point to origin
    o(:) = o(:) - xr(:)
    ! rotate
    tmp1(1) = u(1,1)*o(1) + u(1,2)*o(2) + u(1,3)*o(3)
    tmp1(2) = u(2,1)*o(1) + u(2,2)*o(2) + u(2,3)*o(3)
    tmp1(3) = u(3,1)*o(1) + u(3,2)*o(2) + u(3,3)*o(3)
    ! move origin to new reference point (experimental structure)
    o(:) = tmp1(:) + xs(:)

    ! calculate projection of point B into new coordinate system
    px(:) = xb(:) - o(:)

    ! projection into new coordinate system
    vx(1) = px(1)*u(1,1) + px(2)*u(2,1) + px(3)*u(3,1)
    vx(2) = px(1)*u(1,2) + px(2)*u(2,2) + px(3)*u(3,2)
    vx(3) = px(1)*u(1,3) + px(2)*u(2,3) + px(3)*u(3,3)

    ! radial distance
    select case(cv_item%direction)
        case(1)
            ! x
            ctx%CVsValues(cv_item%idx) = sqrt(vx(2)**2 + vx(3)**2)
        case(2)
            ! y
            ctx%CVsValues(cv_item%idx) = sqrt(vx(1)**2 + vx(3)**2)
        case(3)
            ! z
            ctx%CVsValues(cv_item%idx) = sqrt(vx(1)**2 + vx(2)**2)
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for direction in calculate_orad2!')
    end select

    ! derivatives
    if( ctx%CVsValues(cv_item%idx) .gt. 1e-7 ) then
        sc = 1.0d0 / ctx%CVsValues(cv_item%idx)
    else
        sc = 0.0d0
    end if

    a_vx(:) = 0.0d0
    select case(cv_item%direction)
        case(1)
            ! x
            a_vx(2) = vx(2)*sc
            a_vx(3) = vx(3)*sc
        case(2)
            ! y
            a_vx(1) = vx(1)*sc
            a_vx(3) = vx(3)*sc
        case(3)
            ! z
            a_vx(1) = vx(1)*sc
            a_vx(2) = vx(2)*sc
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for direction in calculate_orad2!')
    end select

    !vx(1) = px(1)*u(1,1) + px(2)*u(2,1) + px(3)*u(3,1)
    !vx(2) = px(1)*u(1,2) + px(2)*u(2,2) + px(3)*u(3,2)
    !vx(3) = px(1)*u(1,3) + px(2)*u(2,3) + px(3)*u(3,3)
    a_px(1) = a_vx(1)*u(1,1) + a_vx(2)*u(1,2) + a_vx(3)*u(1,3)
    a_px(2) = a_vx(1)*u(2,1) + a_vx(2)*u(2,2) + a_vx(3)*u(2,3)
    a_px(3) = a_vx(1)*u(3,1) + a_vx(2)*u(3,2) + a_vx(3)*u(3,3)

    a_u(1,1) = a_vx(1)*px(1)
    a_u(2,1) = a_vx(1)*px(2)
    a_u(3,1) = a_vx(1)*px(3)
    a_u(1,2) = a_vx(2)*px(1)
    a_u(2,2) = a_vx(2)*px(2)
    a_u(3,2) = a_vx(2)*px(3)
    a_u(1,3) = a_vx(3)*px(1)
    a_u(2,3) = a_vx(3)*px(2)
    a_u(3,3) = a_vx(3)*px(3)

    !px(:) = xb(:) - o(:)
    a_o(:) = -a_px(:)
    a_b(:) = a_px(:)

    ! --------------------------------------------------------------------------

    l_u(:,:) = a_u(:,:)

    l_u(1,1) = l_u(1,1) - xr(1)*a_o(1)
    l_u(2,1) = l_u(2,1) - xr(1)*a_o(2)
    l_u(3,1) = l_u(3,1) - xr(1)*a_o(3)
    l_u(1,2) = l_u(1,2) - xr(2)*a_o(1)
    l_u(2,2) = l_u(2,2) - xr(2)*a_o(2)
    l_u(3,2) = l_u(3,2) - xr(2)*a_o(3)
    l_u(1,3) = l_u(1,3) - xr(3)*a_o(1)
    l_u(2,3) = l_u(2,3) - xr(3)*a_o(2)
    l_u(3,3) = l_u(3,3) - xr(3)*a_o(3)

    a_xs(:) = a_o(:)*ingra

    ! rotation matrix a ------------------------------
    a_fa(1) = 2.0d0*( f(1,best)*l_u(1,1) - f(4,best)*l_u(2,1) + f(3,best)*l_u(3,1))
    a_fa(2) = 2.0d0*( f(2,best)*l_u(1,1) + f(3,best)*l_u(2,1) + f(4,best)*l_u(3,1))
    a_fa(3) = 2.0d0*(-f(3,best)*l_u(1,1) + f(2,best)*l_u(2,1) + f(1,best)*l_u(3,1))
    a_fa(4) = 2.0d0*(-f(4,best)*l_u(1,1) - f(1,best)*l_u(2,1) + f(2,best)*l_u(3,1))

    a_fa(1) = a_fa(1) + 2.0d0*(f(4,best)*l_u(1,2) + f(1,best)*l_u(2,2) - f(2,best)*l_u(3,2))
    a_fa(2) = a_fa(2) + 2.0d0*(f(3,best)*l_u(1,2) - f(2,best)*l_u(2,2) - f(1,best)*l_u(3,2))
    a_fa(3) = a_fa(3) + 2.0d0*(f(2,best)*l_u(1,2) + f(3,best)*l_u(2,2) + f(4,best)*l_u(3,2))
    a_fa(4) = a_fa(4) + 2.0d0*(f(1,best)*l_u(1,2) - f(4,best)*l_u(2,2) + f(3,best)*l_u(3,2))

    a_fa(1) = a_fa(1) + 2.0d0*(-f(3,best)*l_u(1,3) + f(2,best)*l_u(2,3) + f(1,best)*l_u(3,3))
    a_fa(2) = a_fa(2) + 2.0d0*( f(4,best)*l_u(1,3) + f(1,best)*l_u(2,3) - f(2,best)*l_u(3,3))
    a_fa(3) = a_fa(3) + 2.0d0*(-f(1,best)*l_u(1,3) + f(4,best)*l_u(2,3) - f(3,best)*l_u(3,3))
    a_fa(4) = a_fa(4) + 2.0d0*( f(2,best)*l_u(1,3) + f(3,best)*l_u(2,3) + f(4,best)*l_u(3,3))

    ! derivatives of fa with respect to matrix elements
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

    ! merge xij with a_fa, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fa(1)*xij(1,mi,mj)+a_fa(2)*xij(2,mi,mj)+a_fa(3)*xij(3,mi,mj)+a_fa(4)*xij(4,mi,mj))*ingra
        end do
    end do

    ! finally gradients for group_a
    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a%cvs(1,i) - xr(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_a%cvs(2,i) - xr(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a%cvs(3,i) - xr(3)) &
                + a_xs(1)

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_a%cvs(1,i) - xr(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a%cvs(2,i) - xr(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a%cvs(3,i) - xr(3)) &
                + a_xs(2)

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a%cvs(1,i) - xr(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a%cvs(2,i) - xr(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_a%cvs(3,i) - xr(3)) &
                + a_xs(3)
    end do

    ! finally gradients for group_b
    do i = cv_item%grps(1)+1, cv_item%grps(2)
        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_b(:)*ingrb
    end do

end subroutine calculate_orad2

!===============================================================================

end module cv_orad2

