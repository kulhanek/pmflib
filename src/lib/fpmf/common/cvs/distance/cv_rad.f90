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

module cv_rad

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeRAD

    integer :: direction

    contains
        procedure :: load_cv        => load_rad
        procedure :: calculate_cv   => calculate_rad
end type CVTypeRAD

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_rad
!===============================================================================

subroutine load_rad(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypeRAD)                    :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(1)                        :: cdir
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'RAD'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_read_groups(cv_item,prm_fin)

    ! load orientation ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'direction',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'direction is not specified!')
    else
        write(PMF_OUT,20) cdir
    end if

    select case(cdir)
        case('x')
            cv_item%direction = 1
        case('y')
            cv_item%direction = 2
        case('z')
            cv_item%direction = 3
        case default
            call pmf_utils_exit(PMF_OUT,1,&
                                'direction has to be x, y, or z!')
    end select

    return

    20  format('   ** z-axis direction   : ',A1)

end subroutine load_rad

!===============================================================================
! Subroutine:  calculate_rad
!===============================================================================

subroutine calculate_rad(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeRAD)    :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: m,ai,di
    real(PMFDP)    :: totmass,pos(3),amass,sc
    ! --------------------------------------------------------------------------

    totmass = 0.0d0
    pos(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        pos(:) = pos(:) + x(:,ai)*amass
        totmass = totmass + amass
    end do

    if( totmass .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass is zero in calculate_rad!')
    end if

    pos(:) = pos(:) / totmass

    ! erase coordinate in "z-direction"
    di = cv_item%direction
    pos(di) = 0.0d0

    ! calculate distance
    ctx%CVsValues(cv_item%idx) = sqrt(pos(1)**2 + pos(2)**2 + pos(3)**2)

    ! calculate forces -------------------------------
    if( ctx%CVsValues(cv_item%idx) .gt. 1e-7 ) then
        sc = 1.0d0 / ctx%CVsValues(cv_item%idx)
    else
        sc = 0.0
    end if

    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + sc*pos(:)*amass/totmass
    end do

end subroutine calculate_rad

!===============================================================================

end module cv_rad

