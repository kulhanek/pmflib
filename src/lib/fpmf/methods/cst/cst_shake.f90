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

module cst_shake

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! logical function cst_shake_checkatom(atomid)
!===============================================================================

logical function cst_shake_checkatom(atomid)

    use pmf_dat
    use cst_dat

    implicit none
    integer    :: atomid
    ! -----------------------------------------------
    integer    :: i
    ! --------------------------------------------------------------------------

    do i=1,NumOfConAtoms
        if( ConAtoms(i) .eq. atomid ) then
            cst_shake_checkatom = .true.
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,*) 'cst_shake_checkatom-> conflict ',atomid
            end if
            return
        end if
    end do

    cst_shake_checkatom = .false.

    return

end function cst_shake_checkatom

!===============================================================================
! Function:  cst_shake_allocate
!===============================================================================

subroutine cst_shake_allocate(num)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer    :: num ! number of shake constraints
    ! -----------------------------------------------
    integer    :: i,alloc_failed
    ! -----------------------------------------------------------------------------

    NumOfSHAKECONs = num
    if( NumOfSHAKECONs .eq. 0 ) return

    allocate(SHAKECONList(NumOfSHAKECONs),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        write(PMF_OUT,*) 'Unable to allocate memory for SHAKE constraints!'
        call pmf_utils_exit(PMF_OUT, 1)
    end if

    do i=1,NumOfSHAKECONs
        SHAKECONList(i)%at1   = 0
        SHAKECONList(i)%at2   = 0
        SHAKECONList(i)%value = 0.0d0
    end do

return

end subroutine cst_shake_allocate

!===============================================================================
! Function:  cst_shake_set
!===============================================================================

subroutine cst_shake_set(id,at1,at2,value)

    use pmf_dat
    use cst_dat

    implicit none
    integer        :: id       ! id of constraint
    integer        :: at1      ! id of first atom
    integer        :: at2      ! id of second atom
    real(PMFDP)    :: value    ! value of DS constraint
    ! -----------------------------------------------------------------------------

    SHAKECONList(id)%at1   = at1
    SHAKECONList(id)%at2   = at2
    SHAKECONList(id)%value = value

return

end subroutine cst_shake_set

!===============================================================================

end module cst_shake

