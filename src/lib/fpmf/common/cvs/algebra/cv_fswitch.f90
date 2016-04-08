!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_fswitch

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeFSWITCH

    integer     :: arg_cv
    real(PMFDP) :: steepness
    real(PMFDP) :: reference

    contains
        procedure :: load_cv        => load_fswitch
        procedure :: calculate_cv   => calculate_fswitch
end type CVTypeFSWITCH

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_fswitch
!===============================================================================

subroutine load_fswitch(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils
    use pmf_unit

    implicit none
    class(CVTypeFSWITCH)                :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    integer                             :: alloc_failed, i
    character(PMF_MAX_CV_NAME)          :: cv_name
    type(UnitType)                      :: steepnessunit
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'FSWITCH'
    cv_item%gradforanycrd = .true.
    cv_item%isalgebraic   = .true.
    call pmf_unit_init(cv_item%unit)
    call cv_common_read_name(cv_item,prm_fin)

    ! load sub specific data  -----------------------
    allocate(cv_item%algebraicidxs(1), stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for algebraicidxs!')
    endif

    if( .not. prmfile_get_string_by_key(prm_fin,'arg',cv_name) ) then
        call pmf_utils_exit(PMF_OUT,1,'Argument CV (arg) is not specified!')
    else
        write(PMF_OUT,5) cv_name
        cv_item%arg_cv = cv_common_find_cv(cv_name)
        cv_item%algebraicidxs(1) = cv_item%arg_cv
    end if

    ! this CV does not have any groups
    cv_item%ngrps = 0

    ! set atoms as a copy from left and right CVs, this information is required fro proper detection of
    ! SHAKE constraints in collisions
    cv_item%natoms = CVList(cv_item%arg_cv)%cv%natoms

    allocate(cv_item%rindexes(cv_item%natoms), &
             cv_item%lindexes(cv_item%natoms), &
             stat = alloc_failed)

    if( alloc_failed .ne. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'Unable to allocate memory for atoms!')
    endif

    do i=1,CVList(cv_item%arg_cv)%cv%natoms
        cv_item%rindexes(i) = CVList(cv_item%arg_cv)%cv%rindexes(i)
    end do


    if( CVList(cv_item%arg_cv)%cv%gradforanycrd .neqv. .true. ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> Argument CV does not have gradforanycrd == true!')
    end if

    ! load rest of setup
    steepnessunit = pmf_unit_power_unit(CVList(cv_item%arg_cv)%cv%unit,-1)

    if( .not. prmfile_get_real8_by_key(prm_fin,'steepness',cv_item%steepness) ) then
        call pmf_utils_exit(PMF_OUT,1,'steepness is not specified!')
    else
        write(PMF_OUT,10) cv_item%steepness,trim(pmf_unit_label(steepnessunit))
    end if
    call pmf_unit_conv_to_ivalue(steepnessunit,cv_item%steepness)

    if( .not. prmfile_get_real8_by_key(prm_fin,'reference',cv_item%reference) ) then
        call pmf_utils_exit(PMF_OUT,1,'reference is not specified!')
    else
        write(PMF_OUT,20) cv_item%reference,trim(pmf_unit_label(LengthUnit))
    end if
    call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%reference)

    return

     5 format('   ** arg CV             : ',A)
    10 format('   ** steepness          : ',E14.5,' [',A,']')
    20 format('   ** reference distance : ',E14.5,' [',A,']')

end subroutine load_fswitch

!===============================================================================
! Subroutine:  calculate_fswitch
!===============================================================================

subroutine calculate_fswitch(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeFSWITCH)    :: cv_item
    real(PMFDP)             :: x(:,:)
    type(CVContextType)     :: ctx
    ! -----------------------------------------------
    integer                 :: cv_idx
    ! --------------------------------------------------------------------------

!     cv_idx = cv_item%idx
!     ctx%CVsValues(cv_idx) = ctx%CVsValues(cv_item%left_cv) + ctx%CVsValues(cv_item%right_cv)
!     ctx%CVsDrvs(:,:,cv_idx) = ctx%CVsDrvs(:,:,cv_item%left_cv) + ctx%CVsDrvs(:,:,cv_item%right_cv)

    ! disable unused variable warning
    ignored_arg__ = size(x) .ne. 0

end subroutine calculate_fswitch

!===============================================================================

end module cv_fswitch

