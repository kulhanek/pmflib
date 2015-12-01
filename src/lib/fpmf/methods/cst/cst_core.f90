!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
!    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_core

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  cst_core_main_lf
!===============================================================================

subroutine cst_core_main_lf

    use cst_constraints
    use cst_lambdas
    use cst_output

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_increment
    call cst_lambdas_calculate
    call cst_core_analyze
    call cst_output_write

end subroutine cst_core_main_lf

!===============================================================================
! Subroutine:  cst_core_main_vv_shake
!===============================================================================

subroutine cst_core_main_vv_shake

    use cst_lambdas
    use cst_output
    use cst_velocities

    implicit none
    ! --------------------------------------------------------------------------

    call cst_lambdas_calculate
    call cst_velocities_correct_a
    call cst_core_analyze
    call cst_output_write

end subroutine cst_core_main_vv_shake

!===============================================================================
! Subroutine:  cst_core_main_vv_rattle
!===============================================================================

subroutine cst_core_main_vv_rattle()

    use cst_constraints
    use cst_velocities

    implicit none
    ! --------------------------------------------------------------------------

    call cst_constraints_increment
    call cst_velocities_correct_b

end subroutine cst_core_main_vv_rattle

!===============================================================================
! Subroutine:  cst_core_analyze
!===============================================================================

subroutine cst_core_analyze

 use pmf_utils
 use pmf_dat
 use cst_dat

 implicit none
 integer                :: i,ci,j,cj,k,info
 real(PMFDP)            :: fzv,isrz,lam,mu
 ! -----------------------------------------------------------------------------

! reset accumulators ---------------------------------------------------
 if ( faccurst .eq. 0 ) then
    faccumulation = 0
    isrztotal = 0
    isrztotals = 0
    faccurst = -1

    lambdatotal(:) = 0.0d0
    lambdatotals(:) = 0.0d0

    if( has_lambdav ) then
        lambdavtotal(:) = 0.0d0
        lambdavtotals(:) = 0.0d0
    end if

    CONList(:)%sdevtot = 0.0d0

    write(CON_OUT,'(A)') '#-------------------------------------------------------------------------------'
    write(CON_OUT,'(A)') '# INFO: ALL ACCUMULATORS WERE RESETED                                           '
    write(CON_OUT,'(A)') '#       PRODUCTION STAGE OF ACCUMULATION IS STARTED                             '
    write(CON_OUT,'(A)') '#-------------------------------------------------------------------------------'
 end if

 ! cummulate results -----------------------------------------------------
 if( faccurst .gt. 0 ) then
    faccurst = faccurst - 1
 end if

 ! it is not production part
 if( fstep .le. 0 ) return

 faccumulation = faccumulation + 1

 ! this will occur for velocity verlet algorithm
 if( faccumulation .le. 0 ) return

 ! calaculate final constraint deviations
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx
    CONList(i)%deviation = get_deviation(CONList(i)%cv,CVContextP%CVsValues(ci),CONList(i)%value)
    CONList(i)%sdevtot = CONList(i)%sdevtot + CONList(i)%deviation**2
 end do

 ! calculate Z matrix --------------------------------------------------------
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx
    do j=1,NumOfCONs
        cj = CONList(j)%cvindx
        fzv = 0.0
        do k=1,NumOfLAtoms
            fzv = fzv + MassInv(k)*dot_product(CVContextP%CVsDrvs(:,k,ci),CVContextP%CVsDrvs(:,k,cj))
        end do
        fz(i,j) = fzv
    end do
 end do

 ! calculate Z determinant ------------------------------------
 if( NumOfCONs .gt. 1 ) then
    ! LU decomposition
    call dgetrf(NumOfCONs,NumOfCONs,fz,NumOfCONs,indx,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[CST] LU decomposition failed in cst_main!')
    end if
    fzdet = 1.0d0
    ! and finaly determinant
    do i=1,NumOfCONs
        if( indx(i) .ne. i ) then
            fzdet = - fzdet * fz(i,i)
        else
            fzdet = fzdet * fz(i,i)
        end if
    end do
 else
    fzdet = fz(1,1)
 end if

 ! calculate metric tensor correction --------------------------------------------------------
 isrz = 1.0/sqrt(fzdet)
 isrztotal  =  isrztotal  + isrz
 isrztotals =  isrztotals + 1.0d0/fzdet     ! e.g. isrz**2

 do i=1,NumOfCONs
    ! lambda ----------------------------------
    lam             = lambda(i)
    lambdatotal(i)  = lambdatotal(i)  + lam
    lambdatotals(i) = lambdatotals(i) + lam**2

    if( has_lambdav ) then
        mu                = lambdav(i)
        lambdavtotal(i)   = lambdavtotal(i)  + mu
        lambdavtotals(i)  = lambdavtotals(i) + mu**2
    end if
 end do

end subroutine cst_core_analyze

!===============================================================================

end module cst_core

