!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module pmf_cvs

use pmf_sizes
use pmf_constants
use pmf_unit

implicit none

! CV values and derivatives ----------------------------------------------------
type CVContextType
    real(PMFDP),pointer         :: CVsValues(:)          ! CVs values
    real(PMFDP),pointer         :: CVsDrvs(:,:,:)        ! CVs derivatives
    real(PMFDP),pointer         :: CVsDrvDrvs(:,:,:,:,:) ! CVs second derivatives lam81
    real(PMFDP),pointer         :: CVsFrc(:)             ! PMF derivatives lam81
end type CVContextType

! core definition of collective variables (CVs) --------------------------------
type CVType
    type(CVType),pointer        :: next             ! next CV in the chain
    character(PMF_MAX_TYPE)     :: ctype            ! type of definition (bond,angle,..)
    type(UnitType)              :: unit             ! CV unit
    integer                     :: idx              ! index of CV
    character(PMF_MAX_CV_NAME)  :: name             ! CV unit name
    integer                     :: natoms           ! number of atoms
    integer,pointer             :: rindexes(:)      ! real atom indexes
    integer,pointer             :: lindexes(:)      ! local atom indexes
    integer                     :: ngrps            ! number of groups
    integer,pointer             :: grps(:)          ! groups boundary
    integer                     :: pathidx          ! if CV is part of path
    integer                     :: nindatoms        ! number of individual atoms
    integer,pointer             :: indlindexes(:)   ! individual local atom indexes
    ! CV abilities -------------------------------
    logical                     :: gradforanycrd    ! is gradient available for any coordinates?
    logical                     :: isalgebraic      ! this CV combines other CVs
    integer,pointer             :: algebraicidxs(:) ! which CV indexes are used in algebra
    logical                     :: processed        ! used by CST - see cst_constraints_calc_fdxp

    contains
        ! executive methods
        procedure   :: reset_cv
        procedure   :: load_cv
        procedure   :: calculate_cv
        procedure   :: free_cv
        ! unit methods
        procedure   :: get_ulabel
        procedure   :: conv_to_ivalue
        procedure   :: get_rvalue
        ! PBC related methods
        procedure   :: is_periodic_cv
        procedure   :: get_period_cv_value
        procedure   :: get_min_cv_value
        procedure   :: get_max_cv_value
        procedure   :: get_average_value
        procedure   :: get_deviation
end type CVType

! list of CV -------------------------------------------------------------------
type CVPointer
    class(CVType),pointer   :: cv           ! cv data
end type CVPointer

integer                     :: NumOfCVs     ! number of CVs
integer                     :: NumOfFakeCVs ! number of fake CVs
type(CVPointer),allocatable :: CVList(:)    ! input definition of CVs

contains

!===============================================================================
! Subroutine:  reset_cv
!===============================================================================

subroutine reset_cv(cv_item)

    use pmf_unit

    implicit none
    class(CVType)           :: cv_item
    ! --------------------------------------------------------------------------

    cv_item%ctype           = ''
    call pmf_unit_init(cv_item%unit)
    cv_item%idx             = 0
    cv_item%name            = ''
    cv_item%natoms          = 0
    cv_item%rindexes        => NULL()
    cv_item%lindexes        => NULL()
    cv_item%ngrps           = 0
    cv_item%grps            => NULL()
    cv_item%pathidx         = 0
    cv_item%gradforanycrd   = .false.
    cv_item%nindatoms       = 0
    cv_item%indlindexes     => NULL()
    cv_item%isalgebraic     = .false.
    cv_item%algebraicidxs   => NULL()
    cv_item%processed     = .false.

end subroutine reset_cv

!===============================================================================
! Subroutine:  load_cv
!===============================================================================

subroutine load_cv(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVType)                       :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    call pmf_utils_exit(PMF_OUT,1, &
            'load_cv not implemented for ' // cv_item%ctype)

    ! disable unused variable warning
    ignored_arg__ = same_type_as(prm_fin,prm_fin)

end subroutine load_cv

!===============================================================================
! Subroutine:  calculate_cv
!===============================================================================

subroutine calculate_cv(cv_item,x,ctx)

    use pmf_utils

    implicit none
    class(CVType)       :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! --------------------------------------------------------------------------

    call pmf_utils_exit(PMF_OUT,1, &
            'calculate_cv not implemented for ' // cv_item%ctype)

    ! disable unused variable warning
    ignored_arg__ = size(x) .ne. 0
    ignored_arg__ = same_type_as(ctx,ctx)

end subroutine calculate_cv

!===============================================================================
! Subroutine:  calculate_cv
!===============================================================================

subroutine free_cv(cv_item)

    use pmf_utils

    implicit none
    class(CVType)       :: cv_item
    ! --------------------------------------------------------------------------

    if( associated(cv_item%rindexes) ) then
        deallocate(cv_item%rindexes)
    end if
    if( associated(cv_item%lindexes) ) then
        deallocate(cv_item%lindexes)
    end if
    if( associated(cv_item%indlindexes) ) then
        deallocate(cv_item%indlindexes)
    end if
    if( associated(cv_item%algebraicidxs) ) then
        deallocate(cv_item%algebraicidxs)
    end if

end subroutine free_cv

!===============================================================================
! Function:   get_ulabel
!===============================================================================

character(PMF_MAX_SUNIT) function get_ulabel(cv_item)

    implicit none
    class(CVType)   :: cv_item
    ! --------------------------------------------------------------------------

    get_ulabel = pmf_unit_label(cv_item%unit)

end function get_ulabel


!===============================================================================
! Function:   conv_to_ivalue
! convert value from external to internal unit
!===============================================================================

subroutine conv_to_ivalue(cv_item,rivalue)

    implicit none
    class(CVType)               :: cv_item
    real(PMFDP),intent(inout)   :: rivalue
    ! --------------------------------------------------------------------------

    rivalue = pmf_unit_get_ivalue(cv_item%unit,rivalue)

end subroutine conv_to_ivalue

!===============================================================================
! Function:   get_rvalue
! convert value from internal to external unit
!===============================================================================

real(PMFDP) function get_rvalue(cv_item,ivalue)

    implicit none
    class(CVType)               :: cv_item
    real(PMFDP),intent(in)      :: ivalue
    ! --------------------------------------------------------------------------

    get_rvalue = pmf_unit_get_rvalue(cv_item%unit,ivalue)

end function get_rvalue

!===============================================================================
! Function:  is_periodic_cv
! is CV periodic?
!===============================================================================

logical function is_periodic_cv(cv_item)

    implicit none
    class(CVType)   :: cv_item
    ! --------------------------------------------------------------------------

    is_periodic_cv = .false.

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function is_periodic_cv

!===============================================================================
! Function:  get_period_cv_value
! get periodicity value for periodic CV
!===============================================================================

real(PMFDP) function get_period_cv_value(cv_item)

    implicit none
    class(CVType)   :: cv_item
    ! --------------------------------------------------------------------------

    get_period_cv_value = 0.0d0
    if( .not. cv_item%is_periodic_cv() ) return

    get_period_cv_value = cv_item%get_max_cv_value() - cv_item%get_min_cv_value()

    ! disable unsued variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_period_cv_value

!===============================================================================
! Function:  get_min_cv_value
! get min value for periodic CV
!===============================================================================

real(PMFDP) function get_min_cv_value(cv_item)

    implicit none
    class(CVType)   :: cv_item
    ! --------------------------------------------------------------------------

    get_min_cv_value = 0.0d0

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_min_cv_value

!===============================================================================
! Function:  get_max_cv_value
! get max value for periodic CV
!===============================================================================

real(PMFDP) function get_max_cv_value(cv_item)

    implicit none
    class(CVType)   :: cv_item
    ! --------------------------------------------------------------------------

    get_max_cv_value = 0.0d0

    ! disable unused variable warning
    ignored_arg__ = same_type_as(cv_item,cv_item)

end function get_max_cv_value

!===============================================================================
! Function:  get_average_value
!
! some periodical variables require special care when calculating their
! average values
!===============================================================================

real(PMFDP) function get_average_value(cv_item,value1,value2)

    implicit none
    class(CVType)   :: cv_item
    real(PMFDP)     :: value1
    real(PMFDP)     :: value2
    ! --------------------------------------------
    real(PMFDP)     :: minv,maxv,vec
    ! --------------------------------------------------------------------------

    if( .not. cv_item%is_periodic_cv() ) then
        get_average_value = 0.5d0*(value1+value2)
        return
    end if

    minv = cv_item%get_min_cv_value()
    maxv = cv_item%get_max_cv_value()

    if( abs(value2-value1) .lt. 0.5d0*(maxv-minv) ) then
        get_average_value = 0.5d0*(value1+value2)
        return
    else  
        ! get vector
        vec = value2 - value1
        ! shift to box center
        vec = vec + 0.5d0*(maxv+minv)
        ! image as point
        vec = vec - (maxv-minv)*floor((vec - minv)/(maxv-minv))
        ! return vector back
        vec = vec - 0.5d0*(maxv+minv)
        ! calculate average
        get_average_value = 0.5d0*vec + value1
        ! image average
        get_average_value = get_average_value - (maxv-minv)*floor((get_average_value - minv)/(maxv-minv))

        ! debug
        ! write(PMF_DEBUG,*) 'get_average_value:',value1,value2,get_average_value
    end if

end function get_average_value

!===============================================================================
! Function:  get_deviation
!
! some periodical variables require special care when calculating their
! deviation values
!===============================================================================

real(PMFDP) function get_deviation(cv_item,value1,value2)

    implicit none
    class(CVType)   :: cv_item
    real(PMFDP)     :: value1
    real(PMFDP)     :: value2
    ! --------------------------------------------
    real(PMFDP)     :: minv,maxv,vec
    ! --------------------------------------------------------------------------

    if( .not. cv_item%is_periodic_cv() ) then
        get_deviation = value1 - value2
        return
    end if

    minv = cv_item%get_min_cv_value()
    maxv = cv_item%get_max_cv_value()

    if( abs(value1-value2) .lt. 0.5d0*(maxv-minv) ) then
        get_deviation = value1 - value2
        return
    else
        ! get vector
        vec = value1 - value2
        ! shift to box center
        vec = vec + 0.5d0*(maxv+minv)
        ! image as point
        vec = vec - (maxv-minv)*floor((vec-minv)/(maxv-minv))
        ! return vector back
        get_deviation = vec - 0.5d0*(maxv+minv)

        ! debug
        ! write(PMF_DEBUG,*) 'get_deviation:',value1,value2,get_deviation
    end if

end function get_deviation

!===============================================================================
! Function:  get_imaged_value
!
! return imaged value for periodic CV
! it returns unmodified value for other CVs
!===============================================================================

real(PMFDP) function get_imaged_value(cv_item,value)

    implicit none
    class(CVType)   :: cv_item
    real(PMFDP)     :: value
    ! --------------------------------------------
    real(PMFDP)     :: minv,maxv
    ! --------------------------------------------------------------------------

    if( .not. cv_item%is_periodic_cv() ) then
        get_imaged_value = value
        return
    end if

    minv = cv_item%get_min_cv_value()
    maxv = cv_item%get_max_cv_value()

    if(  (value .gt. minv) .and. (value .le. maxv) ) then
        get_imaged_value = value
        return
    else
        get_imaged_value = value - (maxv-minv)*floor((value - minv)/(maxv-minv))
        ! debug
        ! write(*,*) 'get_imaged_value:',value,get_imaged_value
    end if

end function get_imaged_value

!===============================================================================

end module pmf_cvs
