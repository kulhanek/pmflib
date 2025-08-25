!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf_constraints

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! logical function abf_constraints_checkatom(atomid)
!===============================================================================

logical function abf_constraints_checkatom(atomid,stage)

    use pmf_dat
    use abf_dat

    implicit none
    integer    :: atomid
    integer    :: stage     ! 0 - check atom for SHAKE exclusion from MD engine
                            ! 1 - check atom for SHAKE inclusion into PMFLib CV list
    ! -----------------------------------------------
    integer    :: i
    ! --------------------------------------------------------------------------

    abf_constraints_checkatom = .false.

    select case(fmdconmode)
        case(0)
            ! do nothing
            return
        case(1)
            ! remove SHAKE in collision from MD engine
            if( stage .eq. 1 ) return
        case(2)
            ! add SHAKE in collision to PMFLib CV list and remove it from MD engine
    end select

    ! test if the atom is in collision with ABF CV
    do i=1,NumOfABFAtoms
        if( ABFAtoms(i) .eq. atomid ) then
            abf_constraints_checkatom = .true.
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,*) 'abf_constraints_checkatom-> conflict ',atomid
            end if
            return
        end if
    end do

    return

end function abf_constraints_checkatom

!===============================================================================
! Function:  abf_constraints_allocate
!===============================================================================

subroutine abf_constraints_allocate(num)

    use pmf_utils
    use pmf_dat
    use abf_dat

    implicit none
    integer    :: num ! number of shake constraints
    ! -----------------------------------------------
    integer    :: i,alloc_failed
    ! -----------------------------------------------------------------------------

    NumOfABFSHAKECONs = num
    if( NumOfABFSHAKECONs .eq. 0 ) return

    allocate(ABFSHAKECONList(NumOfABFSHAKECONs),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        write(PMF_OUT,*) 'Unable to allocate memory for SHAKE constraints!'
        call pmf_utils_exit(PMF_OUT, 1)
    end if

    do i=1,NumOfABFSHAKECONs
        ABFSHAKECONList(i)%at1   = 0
        ABFSHAKECONList(i)%at2   = 0
        ABFSHAKECONList(i)%value = 0.0d0
    end do

return

end subroutine abf_constraints_allocate

!===============================================================================
! Function:  abf_constraints_set
!===============================================================================

subroutine abf_constraints_set(id,at1,at2,value)

    use pmf_dat
    use abf_dat

    implicit none
    integer        :: id       ! id of constraint
    integer        :: at1      ! id of first atom
    integer        :: at2      ! id of second atom
    real(PMFDP)    :: value    ! value of DS constraint
    ! -----------------------------------------------------------------------------

    ABFSHAKECONList(id)%at1   = at1
    ABFSHAKECONList(id)%at2   = at2
    ABFSHAKECONList(id)%value = value

return

end subroutine abf_constraints_set


!===============================================================================
! Subroutine:  abf_constraints_cv_info
!===============================================================================

subroutine abf_constraints_cv_info(shake_item)

    use abf_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit
    use prmfile

    implicit none
    type(ABFTypeSHAKE)  :: shake_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,100) trim(shake_item%cv%name)
    write(PMF_OUT,110) trim(shake_item%cv%ctype)
    write(PMF_OUT,120) shake_item%at1
    write(PMF_OUT,130) shake_item%at2
    write(PMF_OUT,140) shake_item%cv%get_rvalue(shake_item%value), &
                    trim(shake_item%cv%get_ulabel())
    write(PMF_OUT,150) shake_item%cv%get_rvalue(CVContext%CVsValues(shake_item%cvindx)), &
                    trim(shake_item%cv%get_ulabel())
    return

100 format('    ** Name              : ',a)
110 format('    ** Type              : ',a)
120 format('    ** Atom A            : ',I6)
130 format('    ** Atom B            : ',I6)
140 format('    ** Target value      : ',E16.7,' [',A,']')
150 format('    ** Current value     : ',E16.7,' [',A,']')

end subroutine abf_constraints_cv_info

!===============================================================================
! subroutine:  abf_constraints_calc_ZmatInv
!===============================================================================

subroutine abf_constraints_calc_ZmatInv(cvsdrv)

    use pmf_utils
    use abf_dat

    implicit none
    real(PMFDP)         :: cvsdrv(:,:,:)
    ! --------------------------------------------
    integer             :: i,ci,j,cj,k,info
    ! -----------------------------------------------------------------------------

    ! calculate Z matrix
    do i=1,NumOfABFSHAKECONs
        ci = ABFSHAKECONList(i)%cvindx
        do j=1,NumOfABFSHAKECONs
            cj = ABFSHAKECONList(j)%cvindx
            zinvcst(i,j) = 0.0d0
            do k=1,NumOfLAtoms
                zinvcst(i,j) = zinvcst(i,j) + MassInv(k)*dot_product(cvsdrv(:,k,ci),cvsdrv(:,k,cj))
            end do
        end do
    end do

    ! and now its inversion - we will use LAPAC and LU decomposition
    if (NumOfABFSHAKECONs .gt. 1) then
        call dgetrf(NumOfABFSHAKECONs,NumOfABFSHAKECONs,zinvcst,NumOfABFSHAKECONs,indxcst,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] LU decomposition failed in abf_core_calc_Zmat!')
        end if


        call dgetri(NumOfABFSHAKECONs,zinvcst,NumOfABFSHAKECONs,indxcst,vvcst,NumOfABFSHAKECONs,info)
        if( info .ne. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'[ABF] Matrix inversion failed in abf_core_calc_Zmat!')
        end if
    else
        zinvcst(1,1)  = 1.0d0/zinvcst(1,1)
    end if

    return

end subroutine abf_constraints_calc_ZmatInv

!===============================================================================
! subroutine:  abf_constraints_calc_Pmat
!===============================================================================

subroutine abf_constraints_calc_Pmat(cvsdrv)

    use pmf_utils
    use abf_dat

    implicit none
    real(PMFDP)         :: cvsdrv(:,:,:)
    ! --------------------------------------------
    integer             :: i,ci,j,k,l,m,n
    real(PMFDP)         :: mc(NumOfLAtoms*3,NumOfABFSHAKECONs)
    real(PMFDP)         :: c(NumOfABFSHAKECONs,NumOfLAtoms*3)
    real(PMFDP)         :: mh(3*NumOfLAtoms,NumOfABFSHAKECONs)
    ! -----------------------------------------------------------------------------

    call abf_constraints_calc_ZmatInv(cvsdrv)

!    write(*,*) 'Zinv='
!    do i=1,NumOfABFSHAKECONs
!        do j=1,NumOfABFSHAKECONs
!            write(*,'(F10.6,1X)',ADVANCE='NO') zinvcst(i,j)
!        end do
!        write(*,*)
!    end do
!    write(*,*)

    do i=1,NumOfABFSHAKECONs
        ci = ABFSHAKECONList(i)%cvindx
        l = 1
        do j=1,NumOfLAtoms
          !  write(*,*) 'm=',1.0/MassInv(j),MassInv(j)
            do k=1,3
                mc(l,i) = cvsdrv(k,j,ci)
                c(i,l)  = MassInv(j)*cvsdrv(k,j,ci)
                l = l + 1
            end do
        end do
    end do

!    write(*,*) 'MC='
!        l = 1
!        do j=1,NumOfLAtoms
!            do k=1,3
!                do i=1,NumOfABFSHAKECONs
!                    write(*,'(F10.6,1X)',ADVANCE='NO') mc(l,i)
!
!                end do
!                write(*,*)
!                l = l + 1
!        end do
!
!    end do
!    write(*,*)

    ! M^-1 x C^T x Z^-1
    call dgemm('N','N',NumOfLAtoms*3,NumOfABFSHAKECONs,NumOfABFSHAKECONs,1.0d0,mc,NumOfLAtoms*3,&
               zinvcst,NumOfABFSHAKECONs,0.0d0,mh,NumOfLAtoms*3)

!    write(*,*) 'MH='
!        l = 1
!        do j=1,NumOfLAtoms
!            do k=1,3
!                do i=1,NumOfABFSHAKECONs
!                write(*,'(F10.6,1X)',ADVANCE='NO') mh(l,i)
!
!                end do
!                write(*,*)
!                l = l + 1
!        end do
!    end do
!    write(*,*)


    m = 1
    do i=1,NumOfLAtoms
        do j=1,3
            n = 1
            do k=1,NumOfLAtoms
                do l=1,3
                    if( m .eq. n ) then
                        pcst(j,i,l,k) = 1.0
                    else
                        pcst(j,i,l,k) = 0.0
                    end if
                    n = n + 1
                end do
            end do
            m = m + 1
        end do
    end do

!    do i=1,NumOfLAtoms
!        do j=1,3
!            do k=1,NumOfLAtoms
!                do l=1,3
!                    write(*,'(F10.6,1X)',ADVANCE='NO') pcst(j,i,l,k)
!                end do
!            end do
!            write(*,*)
!        end do
!    end do
!    write(*,*)


!    write(*,*) 'C='
!    do i=1,NumOfABFSHAKECONs
!        l = 1
!        do j=1,NumOfLAtoms
!            do k=1,3
!
!                    write(*,'(F10.6,1X)',ADVANCE='NO') c(i,l)
!                    l = l + 1
!                end do
!
!
!        end do
!        write(*,*)
!    end do
!    write(*,*)

    !
    call dgemm('N','N',NumOfLAtoms*3,NumOfLAtoms*3,NumOfABFSHAKECONs,-1.0d0,mh,NumOfLAtoms*3,c,&
                NumOfABFSHAKECONs,1.0d0,pcst,NumOfLAtoms*3)





!    do i=1,NumOfLAtoms
!        do j=1,3
!            do k=1,NumOfLAtoms
!                do l=1,3
!                    write(*,'(F10.6,1X)',ADVANCE='NO') pcst(j,i,l,k)
!                end do
!            end do
!            write(*,*)
!        end do
!    end do

end subroutine abf_constraints_calc_Pmat

!===============================================================================
! subroutine:  abf_constraints_calc_Pmat
!===============================================================================

subroutine abf_constraints_calc_Pfix(frcoldp,frcnewp,cvsdrv)

    use pmf_utils
    use abf_dat

    implicit none
    real(PMFDP)         :: frcoldp(:,:)
    real(PMFDP)         :: frcnewp(:,:)
    real(PMFDP)         :: cvsdrv(:,:,:)
    ! --------------------------------------------
    real(PMFDP)         :: d1, d2, d3
    integer             :: i, ci, j, k, l
    ! -----------------------------------------------------------------------------

   d1 = (frcold(1,1)-frcold(1,2))**2 + (frcold(2,1)-frcold(2,2))**2 + (frcold(3,1)-frcold(3,2))**2


    call dgemv('N',NumOfLAtoms*3,NumOfLAtoms*3,1.0d0,pcst,NumOfLAtoms*3,frcoldp,&
                1,0.0d0,frcnewp,1)

   d2 = (frcnewp(1,1)-frcnewp(1,2))**2 + (frcnewp(2,1)-frcnewp(2,2))**2 + (frcnewp(3,1)-frcnewp(3,2))**2

    do i=1,NumOfABFSHAKECONs
        ci = ABFSHAKECONList(i)%cvindx
        d3 = 0
        do j=1,NumOfLAtoms
            do k=1,3
                d3 = d3 + MassInv(j)*cvsdrv(k,j,ci)*frcnewp(k,j)
            end do
        end do
        write(*,*) 'd3=', d3
    end do


    write(*,*) 'HERE=', d1, d2

    stop

end subroutine abf_constraints_calc_Pfix




!===============================================================================

end module abf_constraints

