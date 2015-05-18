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

! coordination number - CNGFF
! Fermi-like function
! COM of group_a to any atom of group_b

module cv_cngff

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeCNGFF

    real(PMFDP) :: steepness
    real(PMFDP) :: reference

    contains
        procedure :: load_cv        => load_cngff
        procedure :: calculate_cv   => calculate_cngff
end type CVTypeCNGFF

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_cngff
!===============================================================================

subroutine load_cngff(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils
    use pmf_unit

    implicit none
    class(CVTypeCNGFF)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! -----------------------------------------------
    type(UnitType)                      :: steepnessunit
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'CNGFF'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_read_groups(cv_item,prm_fin)

    ! group_a and group_b cannot overlap
    call cv_common_check_grp_overlap(cv_item,1,2)

    ! load rest of setup
    steepnessunit = pmf_unit_power_unit(LengthUnit,-1)

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

 10 format('   ** steepness          : ',E14.5,' [',A,']')
 20 format('   ** reference distance : ',E14.5,' [',A,']')

end subroutine load_cngff

!===============================================================================
! Subroutine:  calculate_cngff
!===============================================================================

subroutine calculate_cngff(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils
    use pmf_pbc

    implicit none
    class(CVTypeCNGFF)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,j,aj
    real(PMFDP)         :: com(3),dx(3),amass,totmass,d2,d,e,e1,v,dv
    ! --------------------------------------------------------------------------

    ctx%CVsValues(cv_item%idx) = 0.0d0

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
        call pmf_utils_exit(PMF_OUT,1,'totmass is zero in calculate_cngff!')
    end if
    com(:) = com(:) / totmass

    do j = cv_item%grps(1)+1, cv_item%grps(2)
        aj = cv_item%lindexes(j)
        dx(:) = x(:,aj) - com(:)

        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dx)
        end if

        d2 = dx(1)**2 + dx(2)**2 + dx(3)**2
        d  = sqrt(d2)
        e  = exp(cv_item%steepness*(d - cv_item%reference))
        e1 = 1.0d0 + e
        ! e1 cannot be zero
        v  = 1.0d0 / e1
        ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + v
        ! for d->0 the derivative should be zero?
        if( d .gt. 1.0e-7 ) then
            dv = - v*v*e*cv_item%steepness/d
            ctx%CVsDrvs(:,aj,cv_item%idx) = ctx%CVsDrvs(:,aj,cv_item%idx) + dv*dx(:)
            do i = 1, cv_item%grps(1)
                ai = cv_item%lindexes(i)
                amass = mass(ai)
                ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) - amass*dv*dx(:)/totmass
            end do
        end if
    end do

end subroutine calculate_cngff

!===============================================================================

end module cv_cngff

