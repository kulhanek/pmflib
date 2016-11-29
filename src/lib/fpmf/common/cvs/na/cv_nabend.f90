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

module cv_nabend

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNABEND

    type(XYZFILE_TYPE)  :: xyz_str_a
    type(XYZFILE_TYPE)  :: xyz_str_b    

    contains
        procedure :: load_cv        => load_nabend
        procedure :: calculate_cv   => calculate_nabend
end type CVTypeNABEND

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_nabend
!===============================================================================

subroutine load_nabend(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeNABEND)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'NABEND'
    cv_item%unit          = AngleUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 4
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

    ! read group c,d ----------------------------------
    write(PMF_OUT,75)
    call cv_common_read_group(cv_item,prm_fin,3)
    if( cv_get_group_natoms(cv_item,3) .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_c must contain at least one atom!')
    end if
    call cv_common_read_group(cv_item,prm_fin,4)
    if( cv_get_group_natoms(cv_item,4) .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_d must contain at least one atom!')
    end if

    return

 50 format('   == Helical axis A =============================')
 60 format('   == Helical axis B =============================')
 75 format('   == Hinge axis  ================================')
 70 format('   ** reference structure: ',A)

100 format('Atom mismatch between group A and reference A atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B and reference B atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)  

end subroutine load_nabend

!===============================================================================
! Subroutine:  calculate_nabend
!===============================================================================

subroutine calculate_nabend(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNABEND) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,mi,mj,ai,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: uza(3),uzb(3),d1(3),d2(3)
    real(PMFDP)         :: arg,tmp12,tmp22,o_tmp12,o_tmp22,o_tmp1vtmp2v,argo_tmp12,argo_tmp22 
    real(PMFDP)         :: haxis(3),asc,f1
    real(PMFDP)         :: tmp1(3),tmp2(3),ingra,ingrb,ingrc,ingrd
    real(PMFDP)         :: a_tmp1(3),a_tmp2(3),a_uza(3),a_uzb(3)
    real(PMFDP)         :: a_fa(4),a_fb(4),a_haxis(3),a_rij(4,4)
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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nabend for point a!')
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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_nabend for point b!')
    end if

    best = 4

    ! rotation matrix a ------------------------------
    uza(1) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
    uza(2) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
    uza(3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    ! rotation matrix b  ------------------------------
    uzb(1) = 2.0d0*( fb(2,best)*fb(4,best) - fb(1,best)*fb(3,best) )
    uzb(2) = 2.0d0*( fb(3,best)*fb(4,best) + fb(1,best)*fb(2,best) )
    uzb(3) = fb(1,best)**2 - fb(2,best)**2 - fb(3,best)**2 + fb(4,best)**2


    ! inverse number of atoms
    ingrc = 1.0d0 / (cv_item%grps(3)-cv_item%grps(2))
    ingrd = 1.0d0 / (cv_item%grps(4)-cv_item%grps(3))

    ! calculate hinge axis
    d1(:) = 0.0
    do  i = cv_item%grps(2)+1, cv_item%grps(3)
        ai = cv_item%lindexes(i)
        d1(:) = d1(:) + x(:,ai)
    end do
    d1(:) = d1(:) * ingrc 
 
    d2(:) = 0.0d0
    do  i = cv_item%grps(3) + 1 , cv_item%grps(4)
        ai = cv_item%lindexes(i)
        d2(:) = d2(:) + x(:,ai)
    end do
    d2(:) = d2(:) * ingrd

    haxis(:) = d2(:) - d1(:)

    ! 'bending'
    tmp1(1) = uza(2)*haxis(3) - uza(3)*haxis(2)
    tmp1(2) = uza(3)*haxis(1) - uza(1)*haxis(3)
    tmp1(3) = uza(1)*haxis(2) - uza(2)*haxis(1)

    tmp2(1) = uzb(2)*haxis(3) - uzb(3)*haxis(2)
    tmp2(2) = uzb(3)*haxis(1) - uzb(1)*haxis(3)
    tmp2(3) = uzb(1)*haxis(2) - uzb(2)*haxis(1)

    tmp12 = tmp1(1)**2 + tmp1(2)**2 + tmp1(3)**2
    tmp22 = tmp2(1)**2 + tmp2(2)**2 + tmp2(3)**2

    o_tmp12 = 1.0d0 / tmp12
    o_tmp22 = 1.0d0 / tmp22

    o_tmp1vtmp2v = sqrt( o_tmp12 * o_tmp22 )

    arg = (tmp1(1)*tmp2(1) + tmp1(2)*tmp2(2) + tmp1(3)*tmp2(3)) * o_tmp1vtmp2v

    asc = - (tmp1(2)*tmp2(3) - tmp1(3)*tmp2(2))*haxis(1) - (tmp1(3)*tmp2(1) - tmp1(1)*tmp2(3))*haxis(2) &
            - (tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1))*haxis(3)

    if ( arg .gt.  1.0 ) then
        arg =  1.0
    else if ( arg .lt. -1.0 ) then
        arg = -1.0
    end if

    ctx%CVsValues(cv_item%idx) = sign(1.0d0,asc)*acos(arg)

    ! first derivatives --------------------------------------------------------------------------------

    argo_tmp12 = arg*o_tmp12;
    argo_tmp22 = arg*o_tmp22;

    f1 = sin(ctx%CVsValues(cv_item%idx))
    if( abs(f1) .lt. 1.e-12 ) then
        ! avoid division by zero
        f1 = -1.e12
    else
        f1 = -1.0d0 / f1
    end if

! with respect to tmp1 and tmp2
    a_tmp1(:) = f1*(tmp2(:)*o_tmp1vtmp2v - tmp1(:)*argo_tmp12)
    a_tmp2(:) = f1*(tmp1(:)*o_tmp1vtmp2v - tmp2(:)*argo_tmp22)

! with respect to uza, uzb, and haxis
    a_uza(1) = haxis(2)*a_tmp1(3) - a_tmp1(2)*haxis(3) 
    a_uza(2) = haxis(3)*a_tmp1(1) - a_tmp1(3)*haxis(1)
    a_uza(3) = haxis(1)*a_tmp1(2) - a_tmp1(1)*haxis(2)

    a_uzb(1) = haxis(2)*a_tmp2(3) - a_tmp2(2)*haxis(3) 
    a_uzb(2) = haxis(3)*a_tmp2(1) - a_tmp2(3)*haxis(1)
    a_uzb(3) = haxis(1)*a_tmp2(2) - a_tmp2(1)*haxis(2) 

    a_haxis(1) = a_tmp1(2)*uza(3) - uza(2)*a_tmp1(3) + a_tmp2(2)*uzb(3) - uzb(2)*a_tmp2(3) 
    a_haxis(2) = a_tmp1(3)*uza(1) - uza(3)*a_tmp1(1) + a_tmp2(3)*uzb(1) - uzb(3)*a_tmp2(1)
    a_haxis(3) = a_tmp1(1)*uza(2) - uza(1)*a_tmp1(2) + a_tmp2(1)*uzb(2) - uzb(1)*a_tmp2(2)
 
! with respect to fa and fb eigenvectors
    
!     uza(1) = 2.0d0*( fa(2,best)*fa(4,best) - fa(1,best)*fa(3,best) )
!     uza(2) = 2.0d0*( fa(3,best)*fa(4,best) + fa(1,best)*fa(2,best) )
!     uza(3) = fa(1,best)**2 - fa(2,best)**2 - fa(3,best)**2 + fa(4,best)**2

    a_fa(1) = 2.0d0*(-fa(3,best)*a_uza(1) + fa(2,best)*a_uza(2) + fa(1,best)*a_uza(3)) 
    a_fa(2) = 2.0d0*( fa(4,best)*a_uza(1) + fa(1,best)*a_uza(2) - fa(2,best)*a_uza(3))
    a_fa(3) = 2.0d0*(-fa(1,best)*a_uza(1) + fa(4,best)*a_uza(2) - fa(3,best)*a_uza(3))
    a_fa(4) = 2.0d0*( fa(2,best)*a_uza(1) + fa(3,best)*a_uza(2) + fa(4,best)*a_uza(3))

    a_fb(1) = 2.0d0*(-fb(3,best)*a_uzb(1) + fb(2,best)*a_uzb(2) + fb(1,best)*a_uzb(3)) 
    a_fb(2) = 2.0d0*( fb(4,best)*a_uzb(1) + fb(1,best)*a_uzb(2) - fb(2,best)*a_uzb(3))
    a_fb(3) = 2.0d0*(-fb(1,best)*a_uzb(1) + fb(4,best)*a_uzb(2) - fb(3,best)*a_uzb(3))
    a_fb(4) = 2.0d0*( fb(2,best)*a_uzb(1) + fb(3,best)*a_uzb(2) + fb(4,best)*a_uzb(3))

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

! finaly gradients for group_c, group_d
    do i = cv_item%grps(2) + 1, cv_item%grps(3)
        ai = cv_item%lindexes(i)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - a_haxis(:)*ingrc
    end do

    do i = cv_item%grps(3) + 1, cv_item%grps(4)
        ai = cv_item%lindexes(i)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + a_haxis(:)*ingrd
    end do


 return

end subroutine calculate_nabend

!===============================================================================

end module cv_nabend

