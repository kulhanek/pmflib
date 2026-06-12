!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2025-2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cst_mdcon_cvs

use pmf_sizes
use pmf_constants

implicit none
contains

!===============================================================================
! logical function cst_mdcon_cvs_checkatom(atomid)
!===============================================================================

logical function cst_mdcon_cvs_checkatom(atomid,stage)

    use pmf_dat
    use cst_dat

    implicit none
    integer     :: atomid
    integer     :: stage
    ! -----------------------------------------------
    integer     :: i
    ! --------------------------------------------------------------------------

    cst_mdcon_cvs_checkatom = .false.

    select case(fmdconmode)
!        case(0)                ! not applicable for CST
!            ! do nothing
!            return
        case(1)
            ! exclude MD constraint in collision from MD engine
            if( stage .eq. 1 ) return
        case(2)
            ! add MD constraint in collision to PMFLib CV list and remove it from MD engine
    end select

    do i=1,NumOfCONAtoms
        if( CONAtoms(i) .eq. atomid ) then
            cst_mdcon_cvs_checkatom = .true.
            if( fdebug ) then
                write(PMF_DEBUG+fmytaskid,*) 'cst_mdcon_cvs_checkatom-> conflict ',atomid
            end if
            return
        end if
    end do

    return

end function cst_mdcon_cvs_checkatom
!===============================================================================
! Function:  cst_mdcon_cvs_allocate
!===============================================================================

subroutine cst_mdcon_cvs_setnumofexcluded(num)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer    :: num ! number of excluded constraints
    ! -----------------------------------------------------------------------------

    NumOfExcMDCONs = num

end subroutine cst_mdcon_cvs_setnumofexcluded

!===============================================================================
! Function:  cst_mdcon_cvs_allocate
!===============================================================================

subroutine cst_mdcon_cvs_allocate(num)

    use pmf_utils
    use pmf_dat
    use cst_dat

    implicit none
    integer    :: num ! number of MD constraints
    ! -----------------------------------------------
    integer    :: i,alloc_failed
    ! -----------------------------------------------------------------------------

    NumOfMDCONs = num
    if( NumOfMDCONs .eq. 0 ) return

    allocate(MDCONList(NumOfMDCONs),stat=alloc_failed)

    if( alloc_failed .ne. 0 ) then
        write(PMF_OUT,*) 'Unable to allocate memory for SHAKE constraints!'
        call pmf_utils_exit(PMF_OUT, 1)
    end if

    do i=1,NumOfMDCONs
        MDCONList(i)%at1   = 0
        MDCONList(i)%at2   = 0
        MDCONList(i)%value = 0.0d0
    end do

return

end subroutine cst_mdcon_cvs_allocate

!===============================================================================
! Function:  cst_mdcon_cvs_set
!===============================================================================

subroutine cst_mdcon_cvs_set(id,at1,at2,value)

    use pmf_dat
    use cst_dat

    implicit none
    integer        :: id       ! id of constraint
    integer        :: at1      ! id of first atom
    integer        :: at2      ! id of second atom
    real(PMFDP)    :: value    ! value of DS constraint
    ! -----------------------------------------------------------------------------

    MDCONList(id)%at1   = at1
    MDCONList(id)%at2   = at2
    MDCONList(id)%value = value

return

end subroutine cst_mdcon_cvs_set

!===============================================================================

end module cst_mdcon_cvs

