!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2016 Ivo Durnik, 
!    Copyright (C) 2016 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_nasbpp

use pmf_sizes
use pmf_constants
use pmf_dat
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNASBPP

    type(XYZFILE_TYPE)  :: xyz_str_a
    type(XYZFILE_TYPE)  :: xyz_str_b
    integer             :: lbp_par     

    contains
        procedure :: load_cv        => load_nasbpp
        procedure :: calculate_cv   => calculate_nasbpp
end type CVTypeNASBPP

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_nasbpp
!===============================================================================

subroutine load_nasbpp(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeNASBPP)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'NASBPP'
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 4
    call cv_common_init_groups(cv_item,prm_fin)

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_a',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_a) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_a)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a)

    if( cv_item%xyz_str_a%natoms .ne. cv_item%grps(1) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group A and reference structure A differs!')
    end if

    do i = 1, cv_item%grps(1)
        ar = cv_item%rindexes(i)
        if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))) .gt. 1.0 ) then
            write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))
            call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
        end if
    end do

    ! read group b ----------------------------------
    write(PMF_OUT,60)
    call cv_common_read_group(cv_item,prm_fin,2)

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_b',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_b) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_b)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_b,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_b)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_b)

    if( cv_item%xyz_str_b%natoms .ne. (cv_item%grps(2) - cv_item%grps(1)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group B and reference structure B differs!')
    end if

    do i = cv_item%grps(1)+1, cv_item%grps(2)
        ar = cv_item%rindexes(i)
        if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
            write(tmpstr,110) i-cv_item%grps(1), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b%symbols(i-cv_item%grps(1)))
            call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
        end if
    end do

    ! read group c,d ----------------------------------
    write(PMF_OUT,75)
    call cv_common_read_group(cv_item,prm_fin,3)
    call cv_common_read_group(cv_item,prm_fin,4)

    ! parameter to be calculated ------------------------
    write(PMF_OUT,80)
    if( .not. prmfile_get_string_by_key(prm_fin,'parameter',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'Type of local base pair parameter (parameter) is not specified!')
    end if
    write(PMF_OUT,90) trim(tmpstr)

    select case(trim(tmpstr))
        case('shear')
            cv_item%unit    = LengthUnit
            cv_item%lbp_par = 1
        case('stretch')
            cv_item%unit    = LengthUnit
            cv_item%lbp_par = 2
        case('stagger')
            cv_item%unit    = LengthUnit
            cv_item%lbp_par = 3
        case('buckle')
            cv_item%unit    = AngleUnit
            cv_item%lbp_par = 4
        case('propeller')
            cv_item%unit    = AngleUnit
            cv_item%lbp_par = 5
        case('opening')
            cv_item%unit    = AngleUnit
            cv_item%lbp_par = 6
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option!')
    end select

    return

 50 format('   == Base #1 ====================================')
 60 format('   == Base #2 ====================================')
 70 format('   ** reference structure: ',A)
 75 format('   == y-axis  ====================================')
 80 format('   -----------------------------------------------') 
 90 format('   ** local BP parameter : ',A)
100 format('Atom mismatch between group A and reference A atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B and reference B atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)  

end subroutine load_nasbpp

!===============================================================================
! Subroutine:  calculate_nasbpp
!===============================================================================

subroutine calculate_nasbpp(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASBPP) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,m,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33,x2sum,xr2sum
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: one_rij2,one_rij
    real(PMFDP)         :: one_rkj2,one_rkj
    real(PMFDP)         :: one_rij2rkj2,one_rijrkj
    real(PMFDP)         :: one_rij3rkj,one_rijrkj3
    real(PMFDP)         :: rij2,rij,rkj2,rkj,arg
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zsc,y0axis(3),asc,d(3)
    real(PMFDP)         :: tmp1(3),tmp2(3),ingra,ingrb, oa(3), ob(3)
    ! --------------------------------------------------------------------------

    ! inverse number of atoms
    ingra = 1.0d0 / cv_item%grps(1)
    ingrb = 1.0d0 / (cv_item%grps(2)-cv_item%grps(1))

    ! calculate geometrical centres (source and target) -------------------
    xsa(:) = 0.0d0
    xra(:) = 0.0d0

    do  i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        ! source
        xsa(:) = xsa(:) + x(:,ai)

        ! reference
        xra(:) = xra(:) + cv_item%xyz_str_a%cvs(:,i)
    end do

    xsa(:) = xsa(:) * ingra
    xra(:) = xra(:) * ingra

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

    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a%cvs(1,i) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a%cvs(2,i) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a%cvs(3,i) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a%cvs(1,i) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a%cvs(2,i) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a%cvs(3,i) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a%cvs(1,i) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a%cvs(2,i) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a%cvs(3,i) - xra(3))
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
    fa(1,1) =  r11 + r22 + r33
    fa(1,2) =  r23 - r32
    fa(1,3) =  r31 - r13
    fa(1,4) =  r12 - r21

    fa(2,1) =  r23 - r32
    fa(2,2) =  r11 - r22 - r33
    fa(2,3) =  r12 + r21
    fa(2,4) =  r13 + r31

    fa(3,1) =  r31 - r13
    fa(3,2) =  r12 + r21
    fa(3,3) = -r11 + r22 - r33
    fa(3,4) =  r23 + r32

    fa(4,1) =  r12 - r21
    fa(4,2) =  r13 + r31
    fa(4,3) =  r23 + r32
    fa(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesa(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fa, 4, eigenvaluesa, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point a!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    xsb(:) = 0.0d0
    xrb(:) = 0.0d0

    do  i = cv_item%grps(1)+1,cv_item%grps(2) 
        ai = cv_item%lindexes(i)

        ! source
        xsb(:) = xsb(:) + x(:,ai)

        ! reference
        xrb(:) = xrb(:) + cv_item%xyz_str_b%cvs(:,i-cv_item%grps(1))
    end do

    xsb(:) = xsb(:) * ingrb
    xrb(:) = xrb(:) * ingrb

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

    do i = cv_item%grps(1)+1,cv_item%grps(2)
        ai = cv_item%lindexes(i)

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b%cvs(3,i-cv_item%grps(1)) - xrb(3))
    end do

    r11 = r11 * ingrb
    r12 = r12 * ingrb
    r13 = r13 * ingrb

    r21 = r21 * ingrb
    r22 = r22 * ingrb
    r23 = r23 * ingrb

    r31 = r31 * ingrb
    r32 = r32 * ingrb
    r33 = r33 * ingrb

    ! construct matrix for quaterion fitting
    fb(1,1) =  r11 + r22 + r33
    fb(1,2) =  r23 - r32
    fb(1,3) =  r31 - r13
    fb(1,4) =  r12 - r21

    fb(2,1) =  r23 - r32
    fb(2,2) =  r11 - r22 - r33
    fb(2,3) =  r12 + r21
    fb(2,4) =  r13 + r31

    fb(3,1) =  r31 - r13
    fb(3,2) =  r12 + r21
    fb(3,3) = -r11 + r22 - r33
    fb(3,4) =  r23 + r32

    fb(4,1) =  r12 - r21
    fb(4,2) =  r13 + r31
    fb(4,3) =  r23 + r32
    fb(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvaluesb(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, fb, 4, eigenvaluesb, work, 26*4, info)

    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nasbpp for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    ua(1,1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
    ua(2,1) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
    ua(3,1) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )

    ua(1,2) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
    ua(2,2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
    ua(3,2) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )

    ua(1,3) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
    ua(2,3) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
    ua(3,3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    ! rotation matrix b  ------------------------------
    ub(1,1) = fb(1,best)**2 + fb(2,best)**2 - fb(3,best)**2 - fb(4,best)**2
    ub(2,1) = 2.0d0*( fb(2,best)*fb(3,best) - fb(1,best)*fb(4,best) )
    ub(3,1) = 2.0d0*( fb(2,best)*fb(4,best) + fb(1,best)*fb(3,best) )

    ub(1,2) = 2.0d0*( fb(2,best)*fb(3,best) + fb(1,best)*fb(4,best) )
    ub(2,2) = fb(1,best)**2 - fb(2,best)**2 + fb(3,best)**2 - fb(4,best)**2
    ub(3,2) = 2.0d0*( fb(3,best)*fb(4,best) - fb(1,best)*fb(2,best) )

    ub(1,3) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
    ub(2,3) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
    ub(3,3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2

    ! z-axis ===================================================================
    ! mutual orientation of two z-axis
    zsc = sign(1.0d0,ua(1,3)*ub(1,3)+ua(2,3)*ub(2,3)+ua(3,3)*ub(3,3))
    ! get z-axis as average of two axes
    zaxis(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    zaxis(:) = zaxis(:)/sqrt(zaxis(1)**2+zaxis(2)**2+zaxis(3)**2)

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
    ! remove projections to z-axis
    yaxis(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:) 
    ! normalize
    yaxis(:) = yaxis(:)/sqrt(yaxis(1)**2+yaxis(2)**2+yaxis(3)**2)

    ! x-axis ===================================================================
    ! is crross product of y and z axes
    xaxis(1) = yaxis(2)*zaxis(3) - yaxis(3)*zaxis(2)
    xaxis(2) = yaxis(3)*zaxis(1) - yaxis(1)*zaxis(3)
    xaxis(3) = yaxis(1)*zaxis(2) - yaxis(2)*zaxis(1)

    ! get origins of bases
    oa(:) = 0.0d0
    ! move reference point to origin
    oa(:) = oa(:) - xra(:)
    ! rotate
    tmp1(1) = ua(1,1)*oa(1) + ua(1,2)*oa(2) + ua(1,3)*oa(3)
    tmp1(2) = ua(2,1)*oa(1) + ua(2,2)*oa(2) + ua(2,3)*oa(3)
    tmp1(3) = ua(3,1)*oa(1) + ua(3,2)*oa(2) + ua(3,3)*oa(3)
    ! move origin to new reference point (experiemntal structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

    write(*,*) oa,ob-oa

    ! vector between base origins
    d(:) = oa(:) - ob(:)

    select case(cv_item%lbp_par)
        case(1)
            ! 'shear'
            ctx%CVsValues(cv_item%idx) = d(1)*xaxis(1) + d(2)*xaxis(2) + d(3)*xaxis(3)
        case(2)
            ! 'stretch'
            ctx%CVsValues(cv_item%idx) = d(1)*yaxis(1) + d(2)*yaxis(2) + d(3)*yaxis(3)
        case(3)
            ! 'stagger'
            ctx%CVsValues(cv_item%idx) = d(1)*zaxis(1) + d(2)*zaxis(2) + d(3)*zaxis(3)
        case(4)
            ! 'buckle'
            tmp1(1) = ua(2,3)*xaxis(3) - ua(3,3)*xaxis(2)
            tmp1(2) = ua(3,3)*xaxis(1) - ua(1,3)*xaxis(3)
            tmp1(3) = ua(1,3)*xaxis(2) - ua(2,3)*xaxis(1)

            tmp2(1) = zsc*ub(2,3)*xaxis(3) - zsc*ub(3,3)*xaxis(2)
            tmp2(2) = zsc*ub(3,3)*xaxis(1) - zsc*ub(1,3)*xaxis(3)
            tmp2(3) = zsc*ub(1,3)*xaxis(2) - zsc*ub(2,3)*xaxis(1)

            arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3))/ &
                (sqrt(tmp1(1)**2+tmp1(2)**2+tmp1(3)**2)*sqrt(tmp2(1)**2+tmp2(2)**2+tmp2(3)**2))

            asc = - (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*xaxis(1) - (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*xaxis(2) &
                  - (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*xaxis(3)

            if ( arg .gt.  1.0 ) then
                arg =  1.0
            else if ( arg .lt. -1.0 ) then
                arg = -1.0
            end if

            ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)
        case(5)
            ! 'propeller'
            tmp1(1) = ua(2,3)*yaxis(3) - ua(3,3)*yaxis(2)
            tmp1(2) = ua(3,3)*yaxis(1) - ua(1,3)*yaxis(3)
            tmp1(3) = ua(1,3)*yaxis(2) - ua(2,3)*yaxis(1)

            tmp2(1) = zsc*ub(2,3)*yaxis(3) - zsc*ub(3,3)*yaxis(2)
            tmp2(2) = zsc*ub(3,3)*yaxis(1) - zsc*ub(1,3)*yaxis(3)
            tmp2(3) = zsc*ub(1,3)*yaxis(2) - zsc*ub(2,3)*yaxis(1)

            arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3))/ &
                (sqrt(tmp1(1)**2+tmp1(2)**2+tmp1(3)**2)*sqrt(tmp2(1)**2+tmp2(2)**2+tmp2(3)**2))

            asc = - (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*yaxis(1) - (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*yaxis(2) &
                  - (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*yaxis(3)

            if ( arg .gt.  1.0 ) then
                arg =  1.0
            else if ( arg .lt. -1.0 ) then
                arg = -1.0
            end if

            ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)
        case(6)
            ! 'opening'
            tmp1(1) = ua(2,2)*zaxis(3) - ua(3,2)*zaxis(2)
            tmp1(2) = ua(3,2)*zaxis(1) - ua(1,2)*zaxis(3)
            tmp1(3) = ua(1,2)*zaxis(2) - ua(2,2)*zaxis(1)

            tmp2(1) = zsc*ub(2,2)*zaxis(3) - zsc*ub(3,2)*zaxis(2)
            tmp2(2) = zsc*ub(3,2)*zaxis(1) - zsc*ub(1,2)*zaxis(3)
            tmp2(3) = zsc*ub(1,2)*zaxis(2) - zsc*ub(2,2)*zaxis(1)

            arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3))/ &
                (sqrt(tmp1(1)**2+tmp1(2)**2+tmp1(3)**2)*sqrt(tmp2(1)**2+tmp2(2)**2+tmp2(3)**2))

            asc = - (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*zaxis(1) - (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*zaxis(2) &
                  - (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*zaxis(3)

            if ( arg .gt.  1.0 ) then
                arg =  1.0
            else if ( arg .lt. -1.0 ) then
                arg = -1.0
            end if

            ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option!')
    end select

 return

end subroutine calculate_nasbpp

!===============================================================================

end module cv_nasbpp

