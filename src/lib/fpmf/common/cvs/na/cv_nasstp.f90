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

module cv_nasstp

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeNASSTP

    type(XYZFILE_TYPE)  :: xyz_str_a1
    type(XYZFILE_TYPE)  :: xyz_str_b1
    type(XYZFILE_TYPE)  :: xyz_str_a2
    type(XYZFILE_TYPE)  :: xyz_str_b2
    integer             :: lst_par

    contains
        procedure :: load_cv        => load_nasstp
        procedure :: calculate_cv   => calculate_nasstp
end type CVTypeNASSTP

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_nasstp
!===============================================================================

subroutine load_nasstp(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_periodic_table

    implicit none
    class(CVTypeNASSTP)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'NASSTP'
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 8
    call cv_common_init_groups(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

!===============================================================================

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    if( cv_get_group_natoms(cv_item,1) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_a must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_a1',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_a1) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_a1)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a1,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a1)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a1)

    if( cv_item%xyz_str_a1%natoms .ne. cv_item%grps(1) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group a and reference structure A1 differs!')
    end if

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a1%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a1%symbols(i))
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
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_b1',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_b1) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_b1)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_b1,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_b1)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_b1)

    if( cv_item%xyz_str_b1%natoms .ne. (cv_item%grps(2) - cv_item%grps(1)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group b and reference structure B1 differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(1)+1, cv_item%grps(2)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b1%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
                write(tmpstr,110) i-cv_item%grps(1), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b1%symbols(i-cv_item%grps(1)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group c,d ----------------------------------
    write(PMF_OUT,75)
    call cv_common_read_group(cv_item,prm_fin,3)
    if( cv_get_group_natoms(cv_item,3) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_c can contain only one atom!')
    end if

    call cv_common_read_group(cv_item,prm_fin,4)
    if( cv_get_group_natoms(cv_item,4) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_d can contain only one atom!')
    end if

!===============================================================================

    ! read group a ----------------------------------
    write(PMF_OUT,55)
    call cv_common_read_group(cv_item,prm_fin,5)

    if( cv_get_group_natoms(cv_item,5) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_e must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_a2',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_a2) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_a2)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_a2,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_a2)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_a2)

    if( cv_item%xyz_str_a2%natoms .ne. (cv_item%grps(5) - cv_item%grps(4) ) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group e and reference structure A2 differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(4)+1, cv_item%grps(5)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_a2%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,120) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str_a2%symbols(i))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,65)
    call cv_common_read_group(cv_item,prm_fin,6)

    if( cv_get_group_natoms(cv_item,6) .le. 3 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_f must contain at least four atoms!')
    end if

    ! read reference structure ------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'reference_b2',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of reference structure (reference_b2) is not specified!')
    end if
    write(PMF_OUT,70) trim(tmpstr)

    call init_xyz(cv_item%xyz_str_b2)
    call open_xyz(PMF_XYZ,tmpstr,cv_item%xyz_str_b2,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str_b2)
    call close_xyz(PMF_XYZ,cv_item%xyz_str_b2)

    if( cv_item%xyz_str_b2%natoms .ne. (cv_item%grps(6) - cv_item%grps(5)) ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in group f and reference structure B2 differs!')
    end if

    if( .not. skiptest ) then
        do i = cv_item%grps(5)+1, cv_item%grps(6)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str_b2%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
                write(tmpstr,130) i-cv_item%grps(1), frmass(ar), SearchMassBySymbol(cv_item%xyz_str_b2%symbols(i-cv_item%grps(1)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    ! read group g,h ----------------------------------
    write(PMF_OUT,76)
    call cv_common_read_group(cv_item,prm_fin,7)
    if( cv_get_group_natoms(cv_item,7) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_g can contain only one atom!')
    end if

    call cv_common_read_group(cv_item,prm_fin,8)
    if( cv_get_group_natoms(cv_item,8) .ne.1 ) then
        call pmf_utils_exit(PMF_OUT,1,'group_h can contain only one atom!')
    end if

!===============================================================================

    ! parameter to be calculated ------------------------
    write(PMF_OUT,80)
    if( .not. prmfile_get_string_by_key(prm_fin,'parameter',tmpstr) ) then
        call pmf_utils_exit(PMF_OUT,1,'Type of local base pair step parameter (parameter) is not specified!')
    end if
    write(PMF_OUT,90) trim(tmpstr)

    select case(trim(tmpstr))
        case('shift')
            cv_item%unit    = LengthUnit
            cv_item%lst_par= 1
        case('slide')
            cv_item%unit    = LengthUnit
            cv_item%lst_par= 2
        case('rise')
            cv_item%unit    = LengthUnit
            cv_item%lst_par= 3
        case('tilt')
            cv_item%unit    = AngleUnit
            cv_item%lst_par= 4
        case('roll')
            cv_item%unit    = AngleUnit
            cv_item%lst_par= 5
        case('twist')
            cv_item%unit    = AngleUnit
            cv_item%lst_par= 6
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unrecognized value for parameter option!')
    end select

    return

 50 format('   == Base #A1 ===================================')
 60 format('   == Base #B1 ===================================')
 55 format('   == Base #A2 ===================================')
 65 format('   == Base #B2 ===================================')
 70 format('   ** reference structure: ',A)
 75 format('   == y1-axis  ===================================')
 76 format('   == y2-axis  ===================================')
 80 format('   -----------------------------------------------')
 90 format('   ** local step parameter : ',A)

100 format('Atom mismatch between group A1 and reference A1 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B1 and reference B1 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

120 format('Atom mismatch between group A2 and reference A2 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
130 format('Atom mismatch between group B2 and reference B2 atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

end subroutine load_nasstp

!===============================================================================
! Subroutine:  calculate_nasstp
!===============================================================================

subroutine calculate_nasstp(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------

    ! --------------------------------------------------------------------------

    !call calculate_nasstp_getbp1(cv_item,x,mst1,morg1)
    !call calculate_nasstp_getbp2(cv_item,x,mst2,morg2)


end subroutine calculate_nasstp

!===============================================================================
! Subroutine:  calculate_nasstp_getbp1
!===============================================================================

subroutine calculate_nasstp_getbp1(cv_item,x,mst,morg)
    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    real(PMFDP)         :: mst(3,3)
    real(PMFDP)         :: morg(3)
    ! -----------------------------------------------
    integer             :: i,ai,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: o_zaxis2,o_zaxis
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zaxisr(3),zsc,y0axis(3)
    real(PMFDP)         :: tmp1(3),ingra,ingrb,oa(3),ob(3)
    real(PMFDP)         :: yaxisr(3),o_yaxis,o_yaxis2
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
        xra(:) = xra(:) + cv_item%xyz_str_a1%cvs(:,i)
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

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))
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
        xrb(:) = xrb(:) + cv_item%xyz_str_b1%cvs(:,i-cv_item%grps(1))
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

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))
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
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_yaxis2 = 1.0d0 / (yaxisr(1)**2 + yaxisr(2)**2 + yaxisr(3)**2)
    o_yaxis  = sqrt(o_yaxis2)
    yaxis(:) = yaxisr(:) * o_yaxis

    ! x-axis ===================================================================
    ! is cross product of y and z axes
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
    ! move origin to new reference point (experimental structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

    ! position of bp origin
    morg(:) = 0.5d0*(oa(:) + ob(:))

end subroutine calculate_nasstp_getbp1

!===============================================================================
! Subroutine:  calculate_nasstp_getbp2
!===============================================================================

subroutine calculate_nasstp_getbp2(cv_item,x,mst,morg)
    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypeNASSTP) :: cv_item
    real(PMFDP)         :: x(:,:)
    real(PMFDP)         :: mst(3,3)
    real(PMFDP)         :: morg(3)
    ! -----------------------------------------------
    integer             :: i,ai,info, best
    real(PMFDP)         :: xsa(3),xra(3),xsb(3),xrb(3)
    real(PMFDP)         :: fa(4,4),fb(4,4),eigenvaluesa(4),eigenvaluesb(4),work(26*4)
    real(PMFDP)         :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)         :: ua(3,3),ub(3,3)
    real(PMFDP)         :: o_zaxis2,o_zaxis
    real(PMFDP)         :: xaxis(3),yaxis(3),zaxis(3),zaxisr(3),zsc,y0axis(3)
    real(PMFDP)         :: tmp1(3),ingra,ingrb,oa(3),ob(3)
    real(PMFDP)         :: yaxisr(3),o_yaxis,o_yaxis2
    ! --------------------------------------------------------------------------

    ! FIXME

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
        xra(:) = xra(:) + cv_item%xyz_str_a1%cvs(:,i)
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

        r11 = r11 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r12 = r12 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r13 = r13 + (x(1,ai) - xsa(1))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r21 = r21 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r22 = r22 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r23 = r23 + (x(2,ai) - xsa(2))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))

        r31 = r31 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(1,i) - xra(1))
        r32 = r32 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(2,i) - xra(2))
        r33 = r33 + (x(3,ai) - xsa(3))*(cv_item%xyz_str_a1%cvs(3,i) - xra(3))
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
        xrb(:) = xrb(:) + cv_item%xyz_str_b1%cvs(:,i-cv_item%grps(1))
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

        r11 = r11 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r12 = r12 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r13 = r13 + (x(1,ai) - xsb(1))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r21 = r21 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r22 = r22 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r23 = r23 + (x(2,ai) - xsb(2))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))

        r31 = r31 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(1,i-cv_item%grps(1)) - xrb(1))
        r32 = r32 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(2,i-cv_item%grps(1)) - xrb(2))
        r33 = r33 + (x(3,ai) - xsb(3))*(cv_item%xyz_str_b1%cvs(3,i-cv_item%grps(1)) - xrb(3))
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
    zaxisr(:) = 0.5d0*ua(:,3) + 0.5d0*zsc*ub(:,3)
    ! normalize
    o_zaxis2 = 1.0d0 / (zaxisr(1)**2 + zaxisr(2)**2 + zaxisr(3)**2)
    o_zaxis  = sqrt(o_zaxis2)
    zaxis(:) = zaxisr(:) * o_zaxis

    ! y-axis ===================================================================
    y0axis(:) = x(:,cv_item%lindexes(cv_item%grps(3))) - x(:,cv_item%lindexes(cv_item%grps(4)))
    ! remove projections to z-axis
    yaxisr(:) = y0axis(:) - (y0axis(1)*zaxis(1)+y0axis(2)*zaxis(2)+y0axis(3)*zaxis(3))*zaxis(:)
    ! normalize
    o_yaxis2 = 1.0d0 / (yaxisr(1)**2 + yaxisr(2)**2 + yaxisr(3)**2)
    o_yaxis  = sqrt(o_yaxis2)
    yaxis(:) = yaxisr(:) * o_yaxis

    ! x-axis ===================================================================
    ! is cross product of y and z axes
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
    ! move origin to new reference point (experimental structure)
    oa(:) = tmp1(:) + xsa(:)

    ob(:) = 0.0d0
    ob(:) = ob(:) - xrb(:)
    tmp1(1) = ub(1,1)*ob(1) + ub(1,2)*ob(2) + ub(1,3)*ob(3)
    tmp1(2) = ub(2,1)*ob(1) + ub(2,2)*ob(2) + ub(2,3)*ob(3)
    tmp1(3) = ub(3,1)*ob(1) + ub(3,2)*ob(2) + ub(3,3)*ob(3)
    ob(:) = tmp1(:) + xsb(:)

    ! position of bp origin
    morg(:) = 0.5d0*(oa(:) + ob(:))

end subroutine calculate_nasstp_getbp2

!===============================================================================

end module cv_nasstp

