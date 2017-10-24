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

! cos(angle) between z-axis of two reference systems superimposed to two atom groups

module cv_caxang2

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeCAXANG2

    type(XYZFILE_TYPE)  :: xyz_str_a
    integer             :: direction_a
    type(XYZFILE_TYPE)  :: xyz_str_b 
    integer             :: direction_b   

    contains
        procedure :: load_cv        => load_caxang2
        procedure :: calculate_cv   => calculate_caxang2
end type CVTypeCAXANG2

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_caxang2
!===============================================================================

subroutine load_caxang2(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeCAXANG2)                :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    character(1)                        :: cdir
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'CAXANG2'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_init_groups(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    if( cv_get_group_natoms(cv_item,1) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_a must contain at least four atoms!')
    end if

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

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a%symbols(i))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    if( .not. prmfile_get_string_by_key(prm_fin,'direction_a',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'direction_a is not specified!')
    end if

    write(PMF_OUT,90) cdir

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

    ! read group b ----------------------------------
    write(PMF_OUT,60)
    call cv_common_read_group(cv_item,prm_fin,2)

    if( cv_get_group_natoms(cv_item,2) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_b must contain at least four atoms!')
    end if

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

    if( .not. skiptest ) then
        do i = cv_item%grps(1)+1, cv_item%grps(2)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
                write(tmpstr,110) i-cv_item%grps(1), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b%symbols(i-cv_item%grps(1)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

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

    return

 50 format('   == Z-axis A ===================================')
 60 format('   == Z-axis B ===================================')
 70 format('   ** reference structure: ',A)
 90 format('   ** reference axis     : ',A)

100 format('Atom mismatch between group A and reference A atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B and reference B atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)  

end subroutine load_caxang2

!===============================================================================
! Subroutine:  calculate_caxang2
!===============================================================================

subroutine calculate_caxang2(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeCAXANG2):: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,mi,mj,ai,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3),ub(3), arg
    real(PMFDP)         :: ingra,ingrb
    real(PMFDP)         :: a_fa(4),a_fb(4),a_rij(4,4)
    real(PMFDP)         :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4)
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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_caxang2 for point a!')
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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_caxang2 for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    select case(cv_item%direction_a )
        case(1)
            ua(1) = fa(1,best)**2 + fa(2,best)**2 - fa(3,best)**2 - fa(4,best)**2
            ua(2) = 2.0d0*( fa(2,best)*fa(3,best) - fa(1,best)*fa(4,best) )
            ua(3) = 2.0d0*( fa(2,best)*fa(4,best) + fa(1,best)*fa(3,best) )
        case(2)
            ua(1) = 2.0d0*( fa(2,best)*fa(3,best) + fa(1,best)*fa(4,best) )
            ua(2) = fa(1,best)**2 - fa(2,best)**2 + fa(3,best)**2 - fa(4,best)**2
            ua(3) = 2.0d0*( fa(3,best)*fa(4,best) - fa(1,best)*fa(2,best) )
        case(3)
            ua(1) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
            ua(2) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
            ua(3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2
        case default
            call pmf_utils_exit(PMF_OUT,1,'direction_a has to be x, y, or z!')
    end select

    ! rotation matrix b  ------------------------------
    select case(cv_item%direction_b )
        case(1)
            ub(1) = fb(1,best)**2 + fb(2,best)**2 - fb(3,best)**2 - fb(4,best)**2
            ub(2) = 2.0d0*( fb(2,best)*fb(3,best) - fb(1,best)*fb(4,best) )
            ub(3) = 2.0d0*( fb(2,best)*fb(4,best) + fb(1,best)*fb(3,best) )
        case(2)
            ub(1) = 2.0d0*( fb(2,best)*fb(3,best) + fb(1,best)*fb(4,best) )
            ub(2) = fb(1,best)**2 - fb(2,best)**2 + fb(3,best)**2 - fb(4,best)**2
            ub(3) = 2.0d0*( fb(3,best)*fb(4,best) - fb(1,best)*fb(2,best) )
        case(3)
            ub(1) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
            ub(2) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
            ub(3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2
        case default
            call pmf_utils_exit(PMF_OUT,1,'direction_a has to be x, y, or z!')
    end select


    arg = ua(1)*ub(1) + ua(2)*ub(2) + ua(3)*ub(3)

    if ( arg .gt.  1.0 ) then
        arg =  1.0
    else if ( arg .lt. -1.0 ) then
        arg = -1.0
    end if

    ctx%CVsValues(cv_item%idx) = arg

    ! first derivatives --------------------------------------------------------------------------------

! with respect to ua
    select case(cv_item%direction_a )
        case(1)
            a_fa(1) = 2.0d0*( fa(1,best)*ub(1) - fa(4,best)*ub(2) + fa(3,best)*ub(3))
            a_fa(2) = 2.0d0*( fa(2,best)*ub(1) + fa(3,best)*ub(2) + fa(4,best)*ub(3))
            a_fa(3) = 2.0d0*(-fa(3,best)*ub(1) + fa(2,best)*ub(2) + fa(1,best)*ub(3))
            a_fa(4) = 2.0d0*(-fa(4,best)*ub(1) - fa(1,best)*ub(2) + fa(2,best)*ub(3))
        case(2)
            a_fa(1) = 2.0d0*( fa(4,best)*ub(1) + fa(1,best)*ub(2) - fa(2,best)*ub(3))
            a_fa(2) = 2.0d0*( fa(3,best)*ub(1) - fa(2,best)*ub(2) - fa(1,best)*ub(3))
            a_fa(3) = 2.0d0*( fa(2,best)*ub(1) + fa(3,best)*ub(2) + fa(4,best)*ub(3))
            a_fa(4) = 2.0d0*( fa(1,best)*ub(1) - fa(4,best)*ub(2) + fa(3,best)*ub(3))
        case(3)
            a_fa(1) = 2.0d0*(-fa(3,best)*ub(1) + fa(2,best)*ub(2) + fa(1,best)*ub(3))
            a_fa(2) = 2.0d0*( fa(4,best)*ub(1) + fa(1,best)*ub(2) - fa(2,best)*ub(3))
            a_fa(3) = 2.0d0*(-fa(1,best)*ub(1) + fa(4,best)*ub(2) - fa(3,best)*ub(3))
            a_fa(4) = 2.0d0*( fa(2,best)*ub(1) + fa(3,best)*ub(2) + fa(4,best)*ub(3))
        case default
            call pmf_utils_exit(PMF_OUT,1,'direction_a has to be x, y, or z!')
    end select

! with respect to ub
    select case(cv_item%direction_b )
        case(1)
            a_fb(1) = 2.0d0*( fb(1,best)*ua(1) - fb(4,best)*ua(2) + fb(3,best)*ua(3))
            a_fb(2) = 2.0d0*( fb(2,best)*ua(1) + fb(3,best)*ua(2) + fb(4,best)*ua(3))
            a_fb(3) = 2.0d0*(-fb(3,best)*ua(1) + fb(2,best)*ua(2) + fb(1,best)*ua(3))
            a_fb(4) = 2.0d0*(-fb(4,best)*ua(1) - fb(1,best)*ua(2) + fb(2,best)*ua(3))
        case(2)
            a_fb(1) = 2.0d0*( fb(4,best)*ua(1) + fb(1,best)*ua(2) - fb(2,best)*ua(3))
            a_fb(2) = 2.0d0*( fb(3,best)*ua(1) - fb(2,best)*ua(2) - fb(1,best)*ua(3))
            a_fb(3) = 2.0d0*( fb(2,best)*ua(1) + fb(3,best)*ua(2) + fb(4,best)*ua(3))
            a_fb(4) = 2.0d0*( fb(1,best)*ua(1) - fb(4,best)*ua(2) + fb(3,best)*ua(3))
        case(3)
            a_fb(1) = 2.0d0*(-fb(3,best)*ua(1) + fb(2,best)*ua(2) + fb(1,best)*ua(3))
            a_fb(2) = 2.0d0*( fb(4,best)*ua(1) + fb(1,best)*ua(2) - fb(2,best)*ua(3))
            a_fb(3) = 2.0d0*(-fb(1,best)*ua(1) + fb(4,best)*ua(2) - fb(3,best)*ua(3))
            a_fb(4) = 2.0d0*( fb(2,best)*ua(1) + fb(3,best)*ua(2) + fb(4,best)*ua(3))
        case default
            call pmf_utils_exit(PMF_OUT,1,'direction_a has to be x, y, or z!')
    end select

! derivatives of fa with respect to matrix elements
    v(:,:) = fa(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvaluesa(i) - eigenvaluesa(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + fa(mj,best)

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

! finaly gradients for group_a
    do i = 1, cv_item%grps(1)

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a%cvs(1,i) - xra(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_a%cvs(2,i) - xra(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a%cvs(3,i) - xra(3))

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_a%cvs(1,i) - xra(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_a%cvs(2,i) - xra(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a%cvs(3,i) - xra(3))

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_a%cvs(1,i) - xra(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_a%cvs(2,i) - xra(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_a%cvs(3,i) - xra(3))

    end do

! derivatives of fb with respect to matrix elements
    v(:,:) = fb(:,:)
    api(:,:) = 0.0d0
    do i=1,4
        if( i .ne. best ) api(i,i) = 1.0d0/(eigenvaluesb(i) - eigenvaluesb(best))
    end do
    call dgemm('N','N',4,4,4,1.0d0,v,4,api,4,0.0d0,bint,4)
    call dgemm('N','T',4,4,4,1.0d0,bint,4,v,4,0.0d0,api,4)

    ! and solve system of equations
    xij(:,:,:) = 0.0d0
    do mi=1,4
        do mj=1,4
            ! construct cij
            cij(:) = 0.0d0
            cij(mi) = cij(mi) + fb(mj,best)

            ! find eigenvector derivatives
            ! xi contains derivatives of eigenvector by A_ij element
            call dgemv('N',4,4,-1.0d0,api,4,cij,1,0.0d0,xij(:,mi,mj),1)
        end do
    end do

! merge xij with a_fb, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = (a_fb(1)*xij(1,mi,mj)+a_fb(2)*xij(2,mi,mj)+a_fb(3)*xij(3,mi,mj)+a_fb(4)*xij(4,mi,mj))*ingrb
        end do
    end do

! finaly gradients for group_b
    do i = cv_item%grps(1) + 1, cv_item%grps(2)

        ai = cv_item%lindexes(i)

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_b%cvs(1,i-cv_item%grps(1)) - xrb(1)) &
                + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_str_b%cvs(2,i-cv_item%grps(1)) - xrb(2)) &
                + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_b%cvs(3,i-cv_item%grps(1)) - xrb(3))

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_str_b%cvs(1,i-cv_item%grps(1)) - xrb(1)) &
                + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_str_b%cvs(2,i-cv_item%grps(1)) - xrb(2)) &
                + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_b%cvs(3,i-cv_item%grps(1)) - xrb(3))

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_str_b%cvs(1,i-cv_item%grps(1)) - xrb(1)) &
                + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_str_b%cvs(2,i-cv_item%grps(1)) - xrb(2)) &
                + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_str_b%cvs(3,i-cv_item%grps(1)) - xrb(3))
    end do

 return

end subroutine calculate_caxang2

!===============================================================================

end module cv_caxang2

