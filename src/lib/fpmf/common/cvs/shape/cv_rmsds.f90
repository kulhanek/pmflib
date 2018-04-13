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

module cv_rmsds

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common
use smf_xyzfile
use smf_xyzfile_type

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeRMSDS

    type(XYZFILE_TYPE)  :: xyz_strA
    type(XYZFILE_TYPE)  :: xyz_strB

    contains
        procedure :: load_cv        => load_rmsds
        procedure :: calculate_cv   => calculate_rmsds
end type CVTypeRMSDS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_rmsds
!===============================================================================

subroutine load_rmsds(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use smf_xyzfile_type
    use smf_xyzfile
    use smf_periodic_table

    implicit none
    class(CVTypeRMSDS)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(len=PRMFILE_MAX_VALUE)    :: file_name
    character(len=PRMFILE_MAX_VALUE)    :: tmpstr
    integer                             :: i,ar
    logical                             :: lresult, skiptest
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'RMSDS'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! this is important for testing
    skiptest = .false.
    lresult = prmfile_get_logical_by_key(prm_fin,'skip_mass_test',skiptest)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_init_groups(cv_item,prm_fin)

    ! read group a ----------------------------------
    write(PMF_OUT,50)
    call cv_common_read_group(cv_item,prm_fin,1)

    ! read target structure -------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'fit',file_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of fit structure (fit) is not specified!')
    end if
    write(PMF_OUT,55) trim(file_name)

    call init_xyz(cv_item%xyz_strA)
    call open_xyz(PMF_XYZ,file_name,cv_item%xyz_strA,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_strA)
    call close_xyz(PMF_XYZ,cv_item%xyz_strA)

    if( cv_item%xyz_strA%natoms .ne. cv_item%grps(1)) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group A and target structure differs!')
    end if

    ! read group b ----------------------------------
    write(PMF_OUT,60)
    call cv_common_read_group(cv_item,prm_fin,2)

    ! read reference structure -------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'target',file_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'File name of target structure (target) is not specified!')
    end if
    write(PMF_OUT,65) trim(file_name)

    call init_xyz(cv_item%xyz_strB)
    call open_xyz(PMF_XYZ,file_name,cv_item%xyz_strB,'OLD')
    call read_xyz(PMF_XYZ,cv_item%xyz_strB)
    call close_xyz(PMF_XYZ,cv_item%xyz_strB)

    if( cv_item%xyz_strB%natoms .ne. (cv_item%grps(2)-cv_item%grps(1))) then
        call pmf_utils_exit(PMF_OUT,1,'Number of atoms in the group B and fit structure differs!')
    end if

    if( .not. skiptest ) then
        do i = 1, cv_item%grps(1)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_strA%symbols(i))) .gt. 1.0 ) then
                write(tmpstr,100) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_strA%symbols(i))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
        do i = cv_item%grps(1)+1, cv_item%grps(2)
            ar = cv_item%rindexes(i)
            if( dabs(frmass(ar) - SearchMassBySymbol(cv_item%xyz_strB%symbols(i-cv_item%grps(1)))) .gt. 1.0 ) then
                write(tmpstr,110) i, frmass(ar), SearchMassBySymbol(cv_item%xyz_strB%symbols(i-cv_item%grps(1)))
                call pmf_utils_exit(PMF_OUT,1,trim(tmpstr))
            end if
        end do
    end if

    return

 50 format('   == Fit structure A ============================')
 55 format('   ** fit structure      : ',A)
 60 format('   == Target structure B =========================')
 65 format('   ** target structure   : ',A)

100 format('Atom mismatch between group A and fit structure atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)
110 format('Atom mismatch between group B and target structure atoms! atom: ',I6,', group mass: ',F10.3, ', ref mass: ',F10.3)

end subroutine load_rmsds

!===============================================================================
! Subroutine:  calculate_rmsds
!===============================================================================

subroutine calculate_rmsds(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeRMSDS)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: i,ai,si,info,best,mi,mj
    real(PMFDP)    :: xA1,xA2,xA3,xrA1,xrA2,xrA3,amass,totmassA,itotmassA
!    real(PMFDP)    :: xB1,xB2,xB3,xrB1,xrB2,xrB3,
    real(PMFDP)    :: totmassB,itotmassB
    real(PMFDP)    :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(PMFDP)    :: f(4,4),u(3,3),a_f(4),a_u(3,3),a_rij(4,4)
    real(PMFDP)    :: eigenvalues(4),work(26*4)
    real(PMFDP)    :: x2sum,xr2sum,sc,fit_rmsdt,value
    real(PMFDP)    :: v(4,4),api(4,4),cij(4),xij(4,4,4),bint(4,4)
    ! -----------------------------------------------------------------------------

    ! calculate geometrical centres (source and target) -------------------
    xA1 = 0.0d0
    xA2 = 0.0d0
    xA3 = 0.0d0
    xrA1 = 0.0d0
    xrA2 = 0.0d0
    xrA3 = 0.0d0
    totmassA = 0.0d0

    do  i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        ! source
        xA1 = xA1 + x(1,ai)*amass
        xA2 = xA2 + x(2,ai)*amass
        xA3 = xA3 + x(3,ai)*amass

        ! reference
        xrA1 = xrA1 + cv_item%xyz_strA%cvs(1,i)*amass
        xrA2 = xrA2 + cv_item%xyz_strA%cvs(2,i)*amass
        xrA3 = xrA3 + cv_item%xyz_strA%cvs(3,i)*amass

        totmassA = totmassA + amass
    end do

    if( totmassA .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmassA is zero in calculate_rmsds!')
    end if

    itotmassA = 1.0d0 / totmassA
    xA1 = xA1 * itotmassA
    xA2 = xA2 * itotmassA
    xA3 = xA3 * itotmassA
    xrA1 = xrA1 * itotmassA
    xrA2 = xrA2 * itotmassA
    xrA3 = xrA3 * itotmassA

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

    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        amass = mass(ai)

        x2sum = x2sum + amass*((x(1,ai) - xA1)**2 &
                             + (x(2,ai) - xA2)**2 &
                             + (x(3,ai) - xA3)**2)
        xr2sum = xr2sum + amass*((cv_item%xyz_strA%cvs(1,i) - xrA1)**2 &
                               + (cv_item%xyz_strA%cvs(2,i) - xrA2)**2 &
                               + (cv_item%xyz_strA%cvs(3,i) - xrA3)**2)

        r11 = r11 + amass*(x(1,ai) - xA1)*(cv_item%xyz_strA%cvs(1,i) - xrA1)
        r12 = r12 + amass*(x(1,ai) - xA1)*(cv_item%xyz_strA%cvs(2,i) - xrA2)
        r13 = r13 + amass*(x(1,ai) - xA1)*(cv_item%xyz_strA%cvs(3,i) - xrA3)

        r21 = r21 + amass*(x(2,ai) - xA2)*(cv_item%xyz_strA%cvs(1,i) - xrA1)
        r22 = r22 + amass*(x(2,ai) - xA2)*(cv_item%xyz_strA%cvs(2,i) - xrA2)
        r23 = r23 + amass*(x(2,ai) - xA2)*(cv_item%xyz_strA%cvs(3,i) - xrA3)

        r31 = r31 + amass*(x(3,ai) - xA3)*(cv_item%xyz_strA%cvs(1,i) - xrA1)
        r32 = r32 + amass*(x(3,ai) - xA3)*(cv_item%xyz_strA%cvs(2,i) - xrA2)
        r33 = r33 + amass*(x(3,ai) - xA3)*(cv_item%xyz_strA%cvs(3,i) - xrA3)
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
        call pmf_utils_exit(PMF_OUT,1,'Unable to diagonalize matrix in calculate_rmsds!')
    end if

    best = 4

    ! rotation matrix ------------------------------
    ! already transposed
    u(1,1) = f(1,best)**2 + f(2,best)**2 - f(3,best)**2 - f(4,best)**2
    u(2,1) = 2.0d0*( f(2,best)*f(3,best) - f(1,best)*f(4,best) )
    u(3,1) = 2.0d0*( f(2,best)*f(4,best) + f(1,best)*f(3,best) )

    u(1,2) = 2.0d0*( f(2,best)*f(3,best) + f(1,best)*f(4,best) )
    u(2,2) = f(1,best)**2 - f(2,best)**2 + f(3,best)**2 - f(4,best)**2
    u(3,2) = 2.0d0*( f(3,best)*f(4,best) - f(1,best)*f(2,best) )

    u(1,3) = 2.0d0*( f(2,best)*f(4,best) - f(1,best)*f(3,best) )
    u(2,3) = 2.0d0*( f(3,best)*f(4,best) + f(1,best)*f(2,best) )
    u(3,3) = f(1,best)**2 - f(2,best)**2 - f(3,best)**2 + f(4,best)**2

    ! value from fitting
    fit_rmsdt = sqrt((x2sum + xr2sum - 2.0d0*eigenvalues(best))*itotmassA)

!     ! calculate geometrical centres (source and target) -------------------
!     xB1 = 0.0d0
!     xB2 = 0.0d0
!     xB3 = 0.0d0
!     xrB1 = 0.0d0
!     xrB2 = 0.0d0
!     xrB3 = 0.0d0
!     totmassB = 0.0d0
! 
!     do  i = cv_item%grps(1)+1, cv_item%grps(2)
!         ai = cv_item%lindexes(i)
!         si = i - cv_item%grps(1)
!         amass = mass(ai)
!         ! source
!         xB1 = xB1 + x(1,ai)*amass
!         xB2 = xB2 + x(2,ai)*amass
!         xB3 = xB3 + x(3,ai)*amass
! 
!         ! reference
!         xrB1 = xrB1 + cv_item%xyz_strB%cvs(1,si)*amass
!         xrB2 = xrB2 + cv_item%xyz_strB%cvs(2,si)*amass
!         xrB3 = xrB3 + cv_item%xyz_strB%cvs(3,si)*amass
! 
!         totmassB = totmassB + amass
!     end do
! 
!     if( totmassB .le. 0 ) then
!         call pmf_utils_exit(PMF_OUT,1,'totmassB is zero in calculate_rmsds!')
!     end if
! 
!     itotmassB = 1.0d0 / totmassB
!     xB1 = xB1 * itotmassB
!     xB2 = xB2 * itotmassB
!     xB3 = xB3 * itotmassB
!     xrB1 = xrB1 * itotmassB
!     xrB2 = xrB2 * itotmassB
!     xrB3 = xrB3 * itotmassB

    a_u(:,:) = 0.0d0
    
    ! calculate rmsds
    value = 0.0d0
    totmassB = 0.0d0
    do  i = cv_item%grps(1)+1, cv_item%grps(2)
        ai = cv_item%lindexes(i)
        si = i - cv_item%grps(1)        
        amass = mass(ai)

        r11 = x(1,ai) - xA1
        r12 = x(2,ai) - xA2
        r13 = x(3,ai) - xA3

        r21 = cv_item%xyz_strB%cvs(1,si) - xrA1
        r22 = cv_item%xyz_strB%cvs(2,si) - xrA2
        r23 = cv_item%xyz_strB%cvs(3,si) - xrA3

        ! transform r2 point
        r31 = u(1,1)*r21 + u(1,2)*r22 + u(1,3)*r23
        r32 = u(2,1)*r21 + u(2,2)*r22 + u(2,3)*r23
        r33 = u(3,1)*r21 + u(3,2)*r22 + u(3,3)*r23

        value = value + amass*( (r11-r31)**2 + (r12-r32)**2 + (r13-r33)**2 )
        
        a_u(1,1) = a_u(1,1) - amass*(r11-r31)*r21
        a_u(1,2) = a_u(1,2) - amass*(r11-r31)*r22
        a_u(1,3) = a_u(1,3) - amass*(r11-r31)*r23

        a_u(2,1) = a_u(2,1) - amass*(r12-r32)*r21
        a_u(2,2) = a_u(2,2) - amass*(r12-r32)*r22
        a_u(2,3) = a_u(2,3) - amass*(r12-r32)*r23
        
        a_u(3,1) = a_u(3,1) - amass*(r13-r33)*r21
        a_u(3,2) = a_u(3,2) - amass*(r13-r33)*r22
        a_u(3,3) = a_u(3,3) - amass*(r13-r33)*r23  
        
        totmassB = totmassB + amass
    end do
    
    itotmassB = 1.0d0 / totmassB
    value = value * itotmassB
    
    ctx%CVsValues(cv_item%idx) = sqrt(value)

    ! calculate derivatives -------------------------------------------------------
    
    ! finaly derivatives group_b
    if( ctx%CVsValues(cv_item%idx) .gt. 1e-7 ) then
        sc = itotmassB / ctx%CVsValues(cv_item%idx)
    else
        sc = 0.0d0
    end if    
    
    ! rotation matrix a ------------------------------
!     ua(1,1) = f(1,best)**2 + f(2,best)**2 - f(3,best)**2 - f(4,best)**2
!     ua(2,1) = 2.0d0*( f(2,best)*f(3,best) - f(1,best)*f(4,best) )
!     ua(3,1) = 2.0d0*( f(2,best)*f(4,best) + f(1,best)*f(3,best) )

    a_f(1) = 2.0d0*( f(1,best)*a_u(1,1) - f(4,best)*a_u(2,1) + f(3,best)*a_u(3,1))
    a_f(2) = 2.0d0*( f(2,best)*a_u(1,1) + f(3,best)*a_u(2,1) + f(4,best)*a_u(3,1))
    a_f(3) = 2.0d0*(-f(3,best)*a_u(1,1) + f(2,best)*a_u(2,1) + f(1,best)*a_u(3,1))
    a_f(4) = 2.0d0*(-f(4,best)*a_u(1,1) - f(1,best)*a_u(2,1) + f(2,best)*a_u(3,1))

!     ua(1,2) = 2.0d0*( f(2,best)*f(3,best) + f(1,best)*f(4,best) )
!     ua(2,2) = f(1,best)**2 - f(2,best)**2 + f(3,best)**2 - f(4,best)**2
!     ua(3,2) = 2.0d0*( f(3,best)*f(4,best) - f(1,best)*f(2,best) )

    a_f(1) = a_f(1) + 2.0d0*(f(4,best)*a_u(1,2) + f(1,best)*a_u(2,2) - f(2,best)*a_u(3,2)) 
    a_f(2) = a_f(2) + 2.0d0*(f(3,best)*a_u(1,2) - f(2,best)*a_u(2,2) - f(1,best)*a_u(3,2))
    a_f(3) = a_f(3) + 2.0d0*(f(2,best)*a_u(1,2) + f(3,best)*a_u(2,2) + f(4,best)*a_u(3,2))
    a_f(4) = a_f(4) + 2.0d0*(f(1,best)*a_u(1,2) - f(4,best)*a_u(2,2) + f(3,best)*a_u(3,2))

!     ua(1,3) = 2.0d0*( f(2,best)*f(4,best) - f(1,best)*f(3,best) )
!     ua(2,3) = 2.0d0*( f(3,best)*f(4,best) + f(1,best)*f(2,best) )
!     ua(3,3) = f(1,best)**2 - f(2,best)**2 - f(3,best)**2 + f(4,best)**2

    a_f(1) = a_f(1) + 2.0d0*(-f(3,best)*a_u(1,3) + f(2,best)*a_u(2,3) + f(1,best)*a_u(3,3)) 
    a_f(2) = a_f(2) + 2.0d0*( f(4,best)*a_u(1,3) + f(1,best)*a_u(2,3) - f(2,best)*a_u(3,3))
    a_f(3) = a_f(3) + 2.0d0*(-f(1,best)*a_u(1,3) + f(4,best)*a_u(2,3) - f(3,best)*a_u(3,3))
    a_f(4) = a_f(4) + 2.0d0*( f(2,best)*a_u(1,3) + f(3,best)*a_u(2,3) + f(4,best)*a_u(3,3))
    
    
! derivatives of f with respect to matrix elements
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
    
! merge xij with a_f, and update by prefactor
    do mi=1,4
        do mj=1,4
            a_rij(mi,mj) = ( a_f(1)*xij(1,mi,mj) + a_f(2)*xij(2,mi,mj)+ &
                             a_f(3)*xij(3,mi,mj) + a_f(4)*xij(4,mi,mj) )
        end do
    end do

    ! finaly gradients for group_a
    do i = 1, cv_item%grps(1)

        ai = cv_item%lindexes(i)
        amass = mass(ai)        

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) &
                + sc*amass*( ( a_rij(1,1)+a_rij(2,2)-a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_strA%cvs(1,i) - xrA1) &
                           + ( a_rij(1,4)+a_rij(2,3)+a_rij(3,2)+a_rij(4,1))*(cv_item%xyz_strA%cvs(2,i) - xrA2) &
                           + (-a_rij(1,3)+a_rij(2,4)-a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_strA%cvs(3,i) - xrA3) ) 

        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) &
                + sc*amass*( (-a_rij(1,4)+a_rij(2,3)+a_rij(3,2)-a_rij(4,1))*(cv_item%xyz_strA%cvs(1,i) - xrA1) &
                           + ( a_rij(1,1)-a_rij(2,2)+a_rij(3,3)-a_rij(4,4))*(cv_item%xyz_strA%cvs(2,i) - xrA2) &
                           + ( a_rij(1,2)+a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_strA%cvs(3,i) - xrA3) )

        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) &
                + sc*amass*( ( a_rij(1,3)+a_rij(2,4)+a_rij(3,1)+a_rij(4,2))*(cv_item%xyz_strA%cvs(1,i) - xrA1) &
                           + (-a_rij(1,2)-a_rij(2,1)+a_rij(3,4)+a_rij(4,3))*(cv_item%xyz_strA%cvs(2,i) - xrA2) &
                           + ( a_rij(1,1)-a_rij(2,2)-a_rij(3,3)+a_rij(4,4))*(cv_item%xyz_strA%cvs(3,i) - xrA3) )

    end do    
    
    
    ! finaly derivatives group_b
    if( ctx%CVsValues(cv_item%idx) .gt. 1e-7 ) then
        sc = itotmassB / ctx%CVsValues(cv_item%idx)
    else
        sc = 0.0d0
    end if

    do  i = cv_item%grps(1)+1, cv_item%grps(2)
        ai = cv_item%lindexes(i)
        si = i - cv_item%grps(1)          
        amass = mass(ai)

        r11 = x(1,ai) - xA1
        r12 = x(2,ai) - xA2
        r13 = x(3,ai) - xA3

        r21 = cv_item%xyz_strB%cvs(1,si) - xrA1
        r22 = cv_item%xyz_strB%cvs(2,si) - xrA2
        r23 = cv_item%xyz_strB%cvs(3,si) - xrA3

        ! transform r2 point
        r31 = u(1,1)*r21 + u(1,2)*r22 + u(1,3)*r23
        r32 = u(2,1)*r21 + u(2,2)*r22 + u(2,3)*r23
        r33 = u(3,1)*r21 + u(3,2)*r22 + u(3,3)*r23

        ctx%CVsDrvs(1,ai,cv_item%idx) = ctx%CVsDrvs(1,ai,cv_item%idx) + amass*sc*(r11 - r31)
        ctx%CVsDrvs(2,ai,cv_item%idx) = ctx%CVsDrvs(2,ai,cv_item%idx) + amass*sc*(r12 - r32)
        ctx%CVsDrvs(3,ai,cv_item%idx) = ctx%CVsDrvs(3,ai,cv_item%idx) + amass*sc*(r13 - r33)
    end do

    return

end subroutine calculate_rmsds

!===============================================================================

end module cv_rmsds

