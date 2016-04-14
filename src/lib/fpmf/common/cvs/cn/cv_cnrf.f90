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

! coordination number - CNRF
! rational function
! any atom of group_a to any atom of group_b

module cv_cnrf

use pmf_sizes
use pmf_constants
use pmf_cvs

implicit none

!===============================================================================

type, extends(CVType) :: CVTypeCNRF

    real(PMFDP) :: offset
    real(PMFDP) :: reference
    integer     :: npow
    integer     :: mpow

    contains
        procedure :: load_cv        => load_cnrf
        procedure :: calculate_cv   => calculate_cnrf
end type CVTypeCNRF

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_cnrf
!===============================================================================

subroutine load_cnrf(cv_item,prm_fin)

    use prmfile
    use pmf_dat
    use cv_common
    use pmf_utils
    use pmf_unit

    implicit none
    class(CVTypeCNRF)                   :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------------------------------------

    ! unit and CV name initialization ---------------
    cv_item%ctype         = 'CNRF'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

    ! load groups -----------------------------------
    cv_item%ngrps = 2
    call cv_common_read_groups(cv_item,prm_fin)

    ! group_a and group_b cannot overlap
    call cv_common_check_grp_overlap(cv_item,1,2)

    ! load rest of setup
    if( .not. prmfile_get_real8_by_key(prm_fin,'offset',cv_item%offset) ) then
        call pmf_utils_exit(PMF_OUT,1,'offset is not specified!')
    else
        write(PMF_OUT,10) cv_item%offset,trim(pmf_unit_label(LengthUnit))
    end if
    call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%offset)

    if( .not. prmfile_get_real8_by_key(prm_fin,'reference',cv_item%reference) ) then
        call pmf_utils_exit(PMF_OUT,1,'reference is not specified!')
    else
        write(PMF_OUT,20) cv_item%reference,trim(pmf_unit_label(LengthUnit))
    end if
    call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%reference)

    if( cv_item%reference .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'reference cannot be zero!')
    end if

    if( .not. prmfile_get_integer_by_key(prm_fin,'npower',cv_item%npow) ) then
        call pmf_utils_exit(PMF_OUT,1,'npower is not specified!')
    else
        write(PMF_OUT,30) cv_item%npow
    end if

    if( cv_item%npow .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'npower cannot be zero!')
    end if

    if( .not. prmfile_get_integer_by_key(prm_fin,'mpower',cv_item%mpow) ) then
        call pmf_utils_exit(PMF_OUT,1,'mpower is not specified!')
    else
        write(PMF_OUT,40) cv_item%mpow
    end if

    if( cv_item%mpow .eq. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'mpower cannot be zero!')
    end if

    return

 10 format('   ** distance offset    : ',E14.5,' [',A,']')
 20 format('   ** reference distance : ',E14.5,' [',A,']')
 30 format('   ** n-power            : ',I8)
 40 format('   ** m-power            : ',I8)

end subroutine load_cnrf

!===============================================================================
! Subroutine:  calculate_cnrf
!===============================================================================

subroutine calculate_cnrf(cv_item,x,ctx)

    use pmf_dat
    use pmf_utils

    implicit none
    class(CVTypeCNRF)   :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,j,aj
    real(PMFDP)         :: dx(3),d2,d,dv,iref,up,dn,dm
    ! --------------------------------------------------------------------------

    ctx%CVsValues(cv_item%idx) = 0.0d0

    iref = 1.0d0 / cv_item%reference

    do i = 1, cv_item%grps(1)
        ai = cv_item%lindexes(i)
        do j = cv_item%grps(1)+1, cv_item%grps(2)
            aj = cv_item%lindexes(j)
            dx(:) = x(:,ai) - x(:,aj)
            d2 = dx(1)**2 + dx(2)**2 + dx(3)**2
            d  = sqrt(d2)
            dm = d - cv_item%offset
            if( dm .le. 0 ) then
                ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + 1.0d0
            else
                if( dm .ne. cv_item%reference ) then
                    up = (dm*iref)**cv_item%npow
                    dn = (dm*iref)**cv_item%mpow
                    ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + (1.0d0 - up) / (1.0d0 - dn)
                    dv =    - cv_item%npow*((dm*iref)**(cv_item%npow-1))*iref*(1.0d0 - dn)
                    dv = dv + cv_item%mpow*((dm*iref)**(cv_item%mpow-1))*iref*(1.0d0 - up)
                    dv = dv / (d * (1.0d0 - dn)**2)
                    ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) + dv*dx(:)
                    ctx%CVsDrvs(:,aj,cv_item%idx) = ctx%CVsDrvs(:,aj,cv_item%idx) - dv*dx(:)
                else
                    ctx%CVsValues(cv_item%idx) = ctx%CVsValues(cv_item%idx) + real(cv_item%npow) / real(cv_item%mpow)
                    ! TODO derivatives ?
                    call pmf_utils_exit(PMF_OUT,1,'not implemented!')
                end if
            end if
        end do
    end do

end subroutine calculate_cnrf

!===============================================================================

end module cv_cnrf

