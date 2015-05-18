!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

! MDISG - minimum distance between two groups of atoms

module cv_mdisg

use pmf_sizes
use pmf_constants
use pmf_dat

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeMDISG

    real(PMFDP) :: beta

    contains
        procedure :: load_cv        => load_mdisg
        procedure :: calculate_cv   => calculate_mdisg
end type CVTypeMDISG

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_mdisg
!===============================================================================

subroutine load_mdisg(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use pmf_utils
    use cv_common

    implicit none
    class(CVTypeMDISG)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! simple init and allocation --------------------
    cv_item%ctype         = 'MDISG'
    cv_item%unit          = LengthUnit
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_read_groups(cv_item,prm_fin)

    ! groups cannot overlap because then minimal distance would be always zero
    call cv_common_check_grp_overlap(cv_item,1,2)

    ! read beta ------------------------------------
    cv_item%beta = 1.0d0
    if( prmfile_get_real8_by_key(prm_fin,'beta',cv_item%beta) ) then
        write(PMF_OUT,10) cv_item%beta
    else
        write(PMF_OUT,15) cv_item%beta
    end if

    10 format('   ** beta               : ',F5.2)
    15 format('   ** beta               : ',F5.2,' (default)')

end subroutine load_mdisg

!===============================================================================
! Subroutine:  calculate_mdisg
!===============================================================================

subroutine calculate_mdisg(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils
    use pmf_pbc

    implicit none
    class(CVTypeMDISG)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer        :: m1,m2,ai1,ai2
    real(PMFDP)    :: d1(3),d2(3),dd(3),dx2,dx,ex,beta,sc,lsc
    ! --------------------------------------------------------------------------

    ! calculate value
    beta = cv_item%beta
    ex = 0.0d0
    do  m1 = 1, cv_item%grps(1)
    ai1 = cv_item%lindexes(m1)
    d1(:) = x(:,ai1)
        do  m2 = cv_item%grps(1)+1, cv_item%grps(2)
            ai2 = cv_item%lindexes(m2)
            d2(:) = x(:,ai2)

            dd(:) = d1(:) - d2(:)

            if( fenable_pbc ) then
                call pmf_pbc_image_vector(dd)
            end if

            dx2 = dd(1)**2 + dd(2)**2 + dd(3)**2
            dx  = sqrt(dx2)
            if( dx .lt. 1.0e-7 ) then
                call pmf_utils_exit(PMF_OUT,1, &
                        'Distance is smaller than 1.0e-7 in calculate_mdisg!')
            end if
            ex = ex + exp( - dx / beta)
        end do
    end do

    ctx%CVsValues(cv_item%idx) = - beta * log(ex)

    ! ------------------------------------------------

    ! calculate gradient
    sc = 1.0d0 / ex
    do  m1 = 1, cv_item%grps(1)
        ai1 = cv_item%lindexes(m1)
        d1(:) = x(:,ai1)
        do  m2 =  cv_item%grps(1)+1, cv_item%grps(2)
            ai2 = cv_item%lindexes(m2)
            d2(:) = x(:,ai2)

            dd(:) = d1(:) - d2(:)
            if( fenable_pbc ) then
                call pmf_pbc_image_vector(dd)
            end if

            dx2 = dd(1)**2 + dd(2)**2 + dd(3)**2
            dx  = sqrt(dx2)
            lsc = sc * exp(- dx / beta) / dx
            ctx%CVsDrvs(:,ai1,cv_item%idx) = ctx%CVsDrvs(:,ai1,cv_item%idx) + lsc*dd(:)
            ctx%CVsDrvs(:,ai2,cv_item%idx) = ctx%CVsDrvs(:,ai2,cv_item%idx) - lsc*dd(:)
        end do
    end do

end subroutine calculate_mdisg

!===============================================================================

end module cv_mdisg

