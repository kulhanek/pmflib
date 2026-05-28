!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module cv_ptpaths2pn

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePTPATHS2PN

    integer             :: nrefsp        ! number of reference points - positive
    integer             :: nrefsn        ! number of reference points - negative
    real(PMFDP)         :: alpha
    real(PMFDP)         :: beta

    contains
        procedure :: load_cv        => load_ptpaths2pn
        procedure :: calculate_cv   => calculate_ptpaths2pn
end type CVTypePTPATHS2PN

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_ptpaths2pn
!===============================================================================

subroutine load_ptpaths2pn(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePTPATHS2PN)             :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------
    character(len=PRMFILE_MAX_LINE)     :: mask
    ! --------------------------------------------------------------------------

! simple init and allocation --------------------
    cv_item%ctype         = 'PTPATHS2PN'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

! read nrefs
    cv_item%nrefsp = 0
    if( prmfile_get_integer_by_key(prm_fin,'nrefsp',cv_item%nrefsp) ) then
        write(PMF_OUT,210) cv_item%nrefsp
    else
        call pmf_utils_exit(PMF_OUT,1,'nrefsp is not specified!')
    end if

    if( cv_item%nrefsp .lt. 0 ) then
       call pmf_utils_exit(PMF_OUT,1,'nrefsp must be greater than or equal to 0!')
    end if

    cv_item%nrefsn = 0
    if( prmfile_get_integer_by_key(prm_fin,'nrefsn',cv_item%nrefsn) ) then
        write(PMF_OUT,215) cv_item%nrefsn
    else
        call pmf_utils_exit(PMF_OUT,1,'nrefsn is not specified!')
    end if

    if( cv_item%nrefsn .lt. 0 ) then
       call pmf_utils_exit(PMF_OUT,1,'nrefsn must be greater than or equal to 0!')
    end if

    if( (cv_item%nrefsp + cv_item%nrefsn) .le. 2 ) then
       call pmf_utils_exit(PMF_OUT,1,'nrefsp + nrefsn must be greater than 2!')
    end if

! init groups -----------------------------------
    cv_item%ngrps = 2 + 1 ! plus anchor point
    call cv_common_init_groups_I(cv_item)

    ! anchor
    mask = 'anchor'
    call cv_common_init_groups_II(cv_item,prm_fin,1,mask)
    mask = 'refpoint'
    call cv_common_init_groups_II(cv_item,prm_fin,2,mask)
    mask = 'direction'
    call cv_common_init_groups_II(cv_item,prm_fin,3,mask)

    call cv_common_init_groups_III(cv_item)

! load groups -----------------------------------
    write(PMF_OUT,100)

    ! anchor
    mask = 'anchor'
    call cv_common_read_group_by_name(cv_item,prm_fin,1,mask)

    write(PMF_OUT,200)

    mask = 'refpoint'
    call cv_common_read_group_by_name(cv_item,prm_fin,2,mask)
    mask = 'direction'
    call cv_common_read_group_by_name(cv_item,prm_fin,3,mask)

! read alpha
    cv_item%alpha = 1.0d0
    if( prmfile_get_real8_by_key(prm_fin,'alpha',cv_item%alpha) ) then
        write(PMF_OUT,220) cv_item%alpha,trim(pmf_unit_label(LengthUnit))
        call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%alpha)
    else
        call pmf_utils_exit(PMF_OUT,1,'alpha is not specified!')
    end if

    cv_item%beta = cv_item%alpha
    if( prmfile_get_real8_by_key(prm_fin,'beta',cv_item%beta) ) then
        write(PMF_OUT,230) cv_item%beta,trim(pmf_unit_label(LengthUnit))
        call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%beta)
    else
        call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%beta)
    end if

   100 format('   == Anchor point ===============================')
   200 format('   == Reference points ===========================')
   210 format('   ** Num of ref. pts (+): ',I6)
   215 format('   ** Num of ref. pts (-): ',I6)
   220 format('   ** Alpha              : ',F5.2,' [',A,']')
   230 format('   ** Beta               : ',F5.2,' [',A,']')

end subroutine load_ptpaths2pn

!===============================================================================
! Subroutine:  calculate_ptpaths2pn
!===============================================================================

subroutine calculate_ptpaths2pn(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils
    use cv_math

    implicit none
    class(CVTypePTPATHS2PN) :: cv_item
    real(PMFDP)             :: x(:,:)
    type(CVContextType)     :: ctx
    ! -----------------------------------------------
    integer             :: i
    real(PMFDP)         :: d1(3),d2(3),d3(3),dv(3),n_dv(3),dx(3),dy(3)
    real(PMFDP)         :: totmass1,totmass2,totmass3
    real(PMFDP)         :: cu,cd,ce,r2,sc1,sc2,sce,sci
    ! --------------------------------------------------------------------------

! calculate CV value
    call get_com(cv_item,1,x,d1,totmass1)
    call get_com(cv_item,2,x,d2,totmass2)
    call get_com(cv_item,3,x,d3,totmass3)

    dv(:) = d3(:) - d2(:)
    call norm_vec(dv,n_dv)

    cu = 0.0d0
    cd = 0.0d0

    do i=-cv_item%nrefsn,cv_item%nrefsp

        dx(:) = d1(:) - (real(i,PMFDP)*cv_item%alpha*n_dv(:) + d2(:))

        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dx)
        end if

        r2 = dx(1)**2 + dx(2)**2 + dx(3)**2

        ce = exp(-r2/cv_item%beta**2)
        cu = cu + real(i,PMFDP)*cv_item%alpha*ce
        cd = cd + ce
    end do

    ctx%CVsValues(cv_item%idx) = cu / cd

 !------------------------------------------------
 !calculate derivatives

    ! (a'b - a*b')/b^2
    sc1 = 1.0 / cd     ! cu'
    sc2 = cu / (cd * cd) ! cd'

    cu = 0.0d0

    do i=-cv_item%nrefsn,cv_item%nrefsp

        dx(:) = d1(:) - (real(i,PMFDP)*cv_item%alpha*n_dv(:) + d2(:))

        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dx)
        end if

        r2 = dx(1)**2 + dx(2)**2 + dx(3)**2

        sce = exp(-r2/cv_item%beta**2)

        sci = real(i,PMFDP)*cv_item%alpha*sce

        cu = - 2.0d0*(sc1*sci - sc2*sce) / cv_item%beta**2

        call get_com_der(cv_item,1,dx,totmass1,cu,ctx)
        call get_com_der(cv_item,2,dx,totmass2,-cu,ctx)

        dy(:) = 0.0d0
        call norm_vec_der(dv,dx,dy)

        call get_com_der(cv_item,3,dy,totmass3,-real(i,PMFDP)*cv_item%alpha*cu,ctx)
        call get_com_der(cv_item,2,dy,totmass2,+real(i,PMFDP)*cv_item%alpha*cu,ctx)
    end do

 return

end subroutine calculate_ptpaths2pn

!===============================================================================

end module cv_ptpaths2pn

