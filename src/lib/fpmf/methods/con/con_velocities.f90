!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

module con_velocities

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! Subroutine:  con_velocities_correct_a
!===============================================================================

subroutine con_velocities_correct_a

 use pmf_dat
 use con_dat
 use con_constraints

 implicit none
 integer        :: i,ci,k
 ! -----------------------------------------------

 ! correct velocities
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx
    do k=1,NumOfLAtoms
        VelP(:,k) = VelP(:,k) - 0.5d0*lambda(i)*PMF_CL2L*fdt*MassInv(k)*CVContext%CVsDrvs(:,k,ci)*PMF_VDT2DT
    end do
 end do

end subroutine con_velocities_correct_a

!===============================================================================
! Subroutine:  con_velocities_correct_b
!===============================================================================

subroutine con_velocities_correct_b

 use pmf_utils
 use pmf_dat
 use con_dat
 use con_constraints

 implicit none
 integer        :: i,ci,j,cj,k,m,info
 real(PMFDP)    :: tmp,tmp1
 character(80)  :: msg
 ! -----------------------------------------------

 ! CVsDrvs has to be ready - in cpmd it is calculated in force part

 ! construct right hand side
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx
    tmp = 0.0d0
    do k=1,NumOfLAtoms
        do m=1,3
            tmp = tmp + CVContext%CVsDrvs(m,k,ci)*VelP(m,k)*PMF_DT2VDT
        end do
    end do
    lambdav(i) = tmp
 end do

 ! construct matrix part
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx
    do j=1,NumOfCONs
        cj = CONList(j)%cvindx
        tmp1 = 0.0d0
        do k=1,NumOfLAtoms
            tmp = 0.0d0
            do m=1,3
                tmp = tmp + CVContext%CVsDrvs(m,k,ci)*CVContext%CVsDrvs(m,k,cj)
            end do
            tmp1 = tmp1 + 0.5d0*fdt*MassInv(k)*tmp
        end do
        matv(i,j) = tmp1
    end do
 end do

 ! solve LE
 if( NumOfCONs .gt. 1 ) then
    call dgetrf(NumOfCONs,NumOfCONs,matv,NumOfCONs,indx,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[CON] LU decomposition failed in con_correct_velocities_b!')
    end if
    call dgetrs('N',NumOfCONs,1,matv,NumOfCONs,indx,lambdav,NumOfCONs,info)
    if( info .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'[CON] Solution of LE failed in con_correct_velocities_b!')
    end if
 else
    lambdav(1) = lambdav(1) / matv(1,1)
 end if

 ! correct velocities
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx 
    do k=1,NumOfLAtoms
        do m=1,3
            VelP(m,k) = VelP(m,k) - 0.5d0*lambdav(i)*fdt*MassInv(k)*CVContext%CVsDrvs(m,k,ci)*PMF_VDT2DT
        end do
    end do
 end do

 ! transform to kcal/mol unit
 lambdav(:) = lambdav(:) * PMF_L2CL

 ! final check of convergence
 do i=1,NumOfCONs
    ci = CONList(i)%cvindx 
    tmp = 0.0d0
    do k=1,NumOfLAtoms
        do m=1,3
            tmp = tmp + CVContext%CVsDrvs(m,k,ci)*VelP(m,k)*PMF_DT2VDT
        end do
    end do
    if( abs(tmp) .gt. 1e-9 ) then
        write(msg,'(A,E16.6,A)') '(',tmp,'>1e-9)'
        call pmf_utils_exit(PMF_OUT,1,'[CON] RATTLE convergence was not achived! '//trim(msg))
    end if
 end do

end subroutine con_velocities_correct_b

!===============================================================================

end module con_velocities

