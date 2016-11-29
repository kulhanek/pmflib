!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2008 Silvia Cereda, sc578@cam.ac.uk &
!                       Petr Kulhanek, kulhanek@enzim.hu
!
!    This library is free software; you can repostribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is postributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, 
!    Boston, MA  02110-1301  USA
!===============================================================================

module cv_epot

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeEPOT
    contains
        procedure :: load_cv        => load_epot
        procedure :: calculate_cv   => calculate_epot
end type CVTypeEPOT

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_epot
!===============================================================================

subroutine load_epot(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypeEPOT)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: i, alloc_failed
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'EPOT'
    cv_item%unit          = EnergyUnit
    cv_item%gradforanycrd = .false.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1

    ! allocate grps array
    allocate(cv_item%grps(cv_item%ngrps),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    cv_item%natoms  = fnatoms
    cv_item%grps(1) = cv_item%natoms

    ! allocate arrays
    allocate(cv_item%rindexes(cv_item%natoms),cv_item%lindexes(cv_item%natoms),stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory!')
    end if

    do i=1,cv_item%natoms
        cv_item%rindexes(i) = i
    end do

end subroutine load_epot

!===============================================================================
! Subroutine:  calculate_epot
!===============================================================================

subroutine calculate_epot(cv_item,x,ctx)

    use pmf_dat

    implicit none
    class(CVTypeEPOT)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------------------------------------

    ctx%CVsValues(cv_item%idx) = PotEne
    ctx%CVsDrvs(:,:,cv_item%idx) = ctx%CVsDrvs(:,:,cv_item%idx) - Frc(:,:)

    ! disable unused variable warning
    ignored_arg__ = size(x) .ne. 0

end subroutine calculate_epot

!===============================================================================

end module cv_epot

