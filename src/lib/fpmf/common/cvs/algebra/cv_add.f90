!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_add

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeADD

    integer :: left_cv
    integer :: right_cv

    contains
        procedure :: load_cv        => load_add
        procedure :: calculate_cv   => calculate_add
end type CVTypeADD

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_add
!===============================================================================

subroutine load_add(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils
    use pmf_unit

    implicit none
    class(CVTypeADD)                    :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: alloc_failed, i
    character(PMF_MAX_CV_NAME)          :: cv_name
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'ADD'
    cv_item%gradforanycrd = .true.
    cv_item%isalgebraic   = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load sub specific data  -----------------------
    allocate(cv_item%algebraicidxs(2), stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for algebraicidxs!')
    endif

    if( .not. prmfile_get_string_by_key(prm_fin,'left',cv_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'Left CV (left) is not specified!')
    else
        write(PMF_OUT,10) cv_name
        cv_item%left_cv = cv_common_find_cv(cv_name)
        cv_item%algebraicidxs(1) = cv_item%left_cv
    end if

    if( .not. prmfile_get_string_by_key(prm_fin,'right',cv_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'Right CV (right) is not specified!')
    else
        write(PMF_OUT,20) cv_name
        cv_item%right_cv = cv_common_find_cv(cv_name)
        cv_item%algebraicidxs(2) = cv_item%right_cv
    end if

    ! this CV does not have any groups
    cv_item%ngrps = 0

    ! set atoms as a copy from left and right CVs, this information is required fro proper detection of
    ! SHAKE constraints in collisions
    cv_item%natoms = CVList(cv_item%left_cv)%cv%natoms + CVList(cv_item%right_cv)%cv%natoms

    allocate(cv_item%rindexes(cv_item%natoms), &
             cv_item%lindexes(cv_item%natoms), &
             stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for atoms!')
    endif

    do i=1,CVList(cv_item%left_cv)%cv%natoms
        cv_item%rindexes(i) = CVList(cv_item%left_cv)%cv%rindexes(i)
    end do

    do i=1,CVList(cv_item%right_cv)%cv%natoms
        cv_item%rindexes(CVList(cv_item%left_cv)%cv%natoms+i) = CVList(cv_item%right_cv)%cv%rindexes(i)
    end do

    ! check units
    if( pmf_unit_compare_units( CVList(cv_item%left_cv)%cv%unit,  CVList(cv_item%right_cv)%cv%unit) .eqv. .false. ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> Specified CVs are not compatible (units differ)!')
    end if

    if( CVList(cv_item%left_cv)%cv%gradforanycrd .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> Left CV does not have gradforanycrd == true!')
    end if

    if( CVList(cv_item%right_cv)%cv%gradforanycrd .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> Right CV does not have gradforanycrd == true!')
    end if

    ! set CV unit
    cv_item%unit = CVList(cv_item%left_cv)%cv%unit

    return

    10 format('   ** left  CV           : ',A)
    20 format('   ** right CV           : ',A)

end subroutine load_add

!===============================================================================
! Subroutine:  calculate_add
!===============================================================================

subroutine calculate_add(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeADD)    :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: cv_idx
    ! --------------------------------------------------------------------------

    cv_idx = cv_item%idx
    ctx%CVsValues(cv_idx) = ctx%CVsValues(cv_item%left_cv) + ctx%CVsValues(cv_item%right_cv)
    ctx%CVsDrvs(:,:,cv_idx) = ctx%CVsDrvs(:,:,cv_item%left_cv) + ctx%CVsDrvs(:,:,cv_item%right_cv)

    ! disable unused variable warning
    ignored_arg__ = size(x) .ne. 0

end subroutine calculate_add

!===============================================================================

end module cv_add

