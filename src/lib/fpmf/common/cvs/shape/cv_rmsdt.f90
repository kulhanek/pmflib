!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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

module cv_rmsdt

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeRMSDT

    type(XYZFILE_TYPE)  :: xyz_str

    contains
        procedure :: load_cv        => load_rmsdt
        procedure :: calculate_cv   => calculate_rmsdt
end type CVTypeRMSDT

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_rmsdt
!===============================================================================

subroutine load_rmsdt(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common
    use smf_xyzfile_type
    use smf_xyzfile
    use smf_periodic_table

    implicit none
    class(CVTypeRMSDT)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_VALUE)    :: file_name
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'RMSDT'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_read_groups(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

    ! read target structure -------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'target',file_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of target structure (target) is not specified!')
    end if
    write(PMF_OUT,50) trim(file_name)

    call init_xyz(cv_item%xyz_str)
    call open_xyz(PMF_XYZ,file_name,cv_item%xyz_str,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_str)
    call close_xyz(PMF_XYZ,cv_item%xyz_str)

    if( cv_item%xyz_str%natoms .ne. cv_item%natoms ) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group and target structure differs!')
    end if

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_str%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_str%symbols(i))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    return

 50 format('   ** target structure   : ',A)
100 format('Atom mismatch between group A and reference A atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

end subroutine load_rmsdt

!===============================================================================
! Subroutine:  calculate_rmsdt
!===============================================================================

subroutine calculate_rmsdt(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeRMSDT)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,ai,info,best
    real(PMFDP)    :: x1,x2,x3,xr1,xr2,xr3,amass,totmass,itotmass
    real(PMFDP)    :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)    :: f(4,4),u(3,3),ut(3,3)
    real(PMFDP)    :: eigenvalues(4),work(26*4)
    real(PMFDP)    :: x2sum,xr2sum,sc
    ! -----------------------------------------------------------------------------

    ! calculate geometrical centres (source and target) -------------------
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

    if( totmass .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass is zero in calculate_rmsdt!')
    end if

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

    x2sum = 0.0d0
    xr2sum = 0.0d0

    do i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)

        x2sum = x2sum + amass*((x(1,ai) - x1)**2 &
                             + (x(2,ai) - x2)**2 &
                             + (x(3,ai) - x3)**2)
        xr2sum = xr2sum + amass*((cv_item%xyz_str%cvs(1,i) - xr1)**2 &
                               + (cv_item%xyz_str%cvs(2,i) - xr2)**2 &
                               + (cv_item%xyz_str%cvs(3,i) - xr3)**2)

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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_rmsdt!')
    end if

    best = 4

    ! calculate rmsdt
    ctx%CVsValues(cv_item%idx) = sqrt((x2sum + xr2sum - 2.0d0*eigenvalues(best))*itotmass)

    ! calculate derivatives -------------------------------------------------------

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

    ! finaly derivatives -------
    if( ctx%CVsValues(cv_item%idx) .gt. 1e-7 ) then
        sc = itotmass / ctx%CVsValues(cv_item%idx)
    else
        sc = 0.0d0
    end if

    do  i = 1, cv_item%natoms
        ai = cv_item%lindexes(i)
        amass = mass(ai)

        r11 = x(1,ai) - x1
        r12 = x(2,ai) - x2
        r13 = x(3,ai) - x3

        r21 = cv_item%xyz_str%cvs(1,i) - xr1
        r22 = cv_item%xyz_str%cvs(2,i) - xr2
        r23 = cv_item%xyz_str%cvs(3,i) - xr3

        ! transform r2 point
        r31 = ut(1,1)*r21 + ut(1,2)*r22 + ut(1,3)*r23
        r32 = ut(2,1)*r21 + ut(2,2)*r22 + ut(2,3)*r23
        r33 = ut(3,1)*r21 + ut(3,2)*r22 + ut(3,3)*r23

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + amass*sc*(r11 - r31)
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + amass*sc*(r12 - r32)
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + amass*sc*(r13 - r33)
    end do

    return

end subroutine calculate_rmsdt

!===============================================================================

end module cv_rmsdt

