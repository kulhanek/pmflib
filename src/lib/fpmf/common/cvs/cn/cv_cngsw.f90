!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

! coordination number - CNGSW
! switch function
! COM of group_a to any atom of group_b

module cv_cngsw

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeCNGSW

    real(PMFDP) :: cutoff

    contains
        procedure :: load_cv        => load_cngsw
        procedure :: calculate_cv   => calculate_cngsw
end type CVTypeCNGSW

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_cngsw
!===============================================================================

subroutine load_cngsw(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils
    use pmf_unit

    implicit none
    class(CVTypeCNGSW)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'CNGSW'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_read_groups(cv_item,prm_fin)

    ! group_a and group_b cannot overlap
    call cv_common_check_grp_overlap(cv_item,1,2)

    if( .not. prmfile_get_real8_by_key(prm_fin,'cutoff',cv_item%cutoff) ) then
        call pmf_utils_exit(PMF_OUT,1,'cutoff is not specified!')
    else
        write(PMF_OUT,10) cv_item%cutoff,trim(pmf_unit_label(LengthUnit))
    end if
    call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%cutoff)

    return

 10 format('   ** cutoff             : ',E14.5,' [',A,']')

end subroutine load_cngsw

!===============================================================================
! Subroutine:  calculate_cngsw
!===============================================================================

subroutine calculate_cngsw(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils
    use pmf_pbc

    implicit none
    class(CVTypeCNGSW)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,j,aj
    real(PMFDP)         :: com(3),dx(3),amass,totmass,d2,v,dv,cutoff2
    ! --------------------------------------------------------------------------

    ctx%CVsValues(cv_item%idx) = 0.0d0
    cutoff2 = cv_item%cutoff**2

    ! calculate com
    com = 0.0d0
    totmass = 0.0d0
    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        amass = mass(ai)
        totmass = totmass + amass
        com(:) = com(:) + amass*x(:,ai)
    end do
    if( totmass .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass is zero in calculate_cngsw!')
    end if
    com(:) = com(:) / totmass

    do j = cv_item%grps(1)+1, cv_item%grps(2)
        aj = cv_item%lindexes(j)
        dx(:) = x(:,aj) - com(:)

        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dx)
        end if

        d2 = dx(1)**2 + dx(2)**2 + dx(3)**2

        if( d2 .ge. cutoff2 ) cycle

        v = 1.0d0 - 2.0d0*d2/cutoff2 + d2**2/cutoff2**2
        ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + v

        ! derivatives
        dv = -4.0d0/cutoff2 + 4.0d0*d2/cutoff2**2

        ctx%CVsDrvs(:,aj,cv_item%idx) = ctx%CVsDrvs(:,aj,cv_item%idx) + dv*dx(:)

        do i = 1, cv_item%grps(1)
            ai = cv_item%lindexes(i)
            amass = mass(ai)
            ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - amass*dv*dx(:)/totmass
        end do
    end do

end subroutine calculate_cngsw

!===============================================================================

end module cv_cngsw

