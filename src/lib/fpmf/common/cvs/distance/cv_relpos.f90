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

module cv_relpos

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeRELPOS

    integer :: atom
    integer :: direction
    real(PMFDP) :: length

    contains
        procedure :: load_cv        => load_relpos
        procedure :: calculate_cv   => calculate_relpos
end type CVTypeRELPOS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_relpos
!===============================================================================

subroutine load_relpos(cv_item,prm_fin)

    use prmfile
    use pmf_utils
    use pmf_dat
    use cv_common

    implicit none
    class(CVTypeRELPOS)                 :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    character(1)                        :: cdir
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'RELPOS'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 1
    call cv_common_read_groups(cv_item,prm_fin)

    ! load atom -------------------------------------
    if( .not. prmfile_get_integer_by_key(prm_fin,'atom',cv_item%atom) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: atom is not specified!')
    else
        write(PMF_OUT,10) cv_item%atom
    end if

    ! load length -----------------------------------
    if( .not. prmfile_get_real8_by_key(prm_fin,'atom',cv_item%length) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: atom is not specified!')
    else
        write(PMF_OUT,20) cv_item%length, trim(pmf_unit_label(LengthUnit)) 
    end if

    ! load orientation ------------------------------
    if( .not. prmfile_get_string_by_key(prm_fin,'direction',cdir) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: direction is not specified!')
    else
        write(PMF_OUT,30) cdir
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
                                '>>> ERROR: direction has to be x, y, or z!')
    end select

    return

    10  format('   ** atom : ',I6)
    20  format('   ** length : ',E14.5,' [',A,']')
    30  format('   ** axis direction   : ',A1)

end subroutine load_relpos

!===============================================================================
! Subroutine:  calculate_relpos
!===============================================================================

subroutine calculate_relpos(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeRELPOS) :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: m,ai,di
    real(PMFDP)    :: totmass,rcom(3),amass
    real(PMFDP)    :: dis,dis2,der
    ! --------------------------------------------------------------------------

    totmass = 0.0d0
    rcom(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        rcom(:) = rcom(:) + x(:,ai)*amass
        totmass = totmass + amass
    end do

    if( totmass .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: totmass is zero in calculate_relpos!')
    end if

    rcom(:) = rcom(:) / totmass

    dis = sqrt( (x(1,cv_item%atom)-rcom(1))**2 + (x(2,cv_item%atom)-rcom(2))**2 + (x(3,cv_item%atom)-rcom(3))**2 )
    dis2 = dis*dis

    di = cv_item%direction
    ctx%CVsValues(cv_item%idx) = rcom(di) + ( x(di,cv_item%atom) - rcom(di) ) / dis * cv_item%length

    ! calculate forces -------------------------------
    der = cv_item%length * ( dis - (x(di,cv_item%atom) - rcom(di)) / dis ) / dis2

    ctx%CVsDrvs(di,cv_item%atom,cv_item%idx) = ctx%CVsDrvs(di,cv_item%atom,cv_item%idx) + der

    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        ctx%CVsDrvs(di,ai,cv_item%idx) = ctx%CVsDrvs(di,ai,cv_item%idx) + amass/totmass * (1.0d0 - der)
    end do

end subroutine calculate_relpos

!===============================================================================

end module cv_relpos

