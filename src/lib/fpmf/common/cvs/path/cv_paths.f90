!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!
!    This library is free software; you can repathstribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is pathstributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
!===============================================================================

module cv_paths

use pmf_sizes
use pmf_constants
use pmf_dat
use cv_common

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePATHS

    integer             :: nrefs        ! number of reference points
    real(PMFDP)         :: alpha

    contains
        procedure :: load_cv        => load_paths
        procedure :: calculate_cv   => calculate_paths
end type CVTypePATHS

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_paths
!===============================================================================

subroutine load_paths(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePATHS)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------
    integer                             :: i
    character(len=PRMFILE_MAX_LINE)     :: mask,key,title
    logical                             :: ret
    ! --------------------------------------------------------------------------

! simple init and allocation --------------------
    cv_item%ctype         = 'PATHS'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

! determine number of reference points ----------
    ret = prmfile_first_line(prm_fin)
    cv_item%nrefs = 0
    do while( prmfile_get_field_on_line(prm_fin,key) .and. ret)
        if( trim(key) .eq. 'refpoint' )  cv_item%nrefs = cv_item%nrefs + 1
        ret = prmfile_next_line(prm_fin)
    end do

! init groups -----------------------------------
    cv_item%ngrps = cv_item%nrefs + 1 ! plus anchor point
    call cv_common_init_groups_I(cv_item)

    ! anchor
    mask = 'anchor'
    call cv_common_init_groups_II(cv_item,prm_fin,1,mask)

    ! refopints
    ret = prmfile_first_line(prm_fin)
    i = 1
    do while( prmfile_get_field_on_line(prm_fin,key) .and. ret)
        if( trim(key) .eq. 'refpoint' ) then
            ret = prmfile_get_string_value_on_line(prm_fin,mask)
            if( .not. ret ) then
                call pmf_utils_exit(PMF_OUT,1,'refpoint error')
            end if
            call cv_common_init_groups_II_bymask(cv_item,1+i,mask)
            i = i + 1
        end if
        ret = prmfile_next_line(prm_fin)
    end do

    call cv_common_init_groups_III(cv_item)

! load groups -----------------------------------
    write(PMF_OUT,100)

    ! anchor
    mask = 'anchor'
    call cv_common_read_group_by_name(cv_item,prm_fin,1,mask)

    ! refpoints
    write(PMF_OUT,200)
    write(PMF_OUT,210) cv_item%nrefs
    ret = prmfile_first_line(prm_fin)
    i = 1
    do while( prmfile_get_field_on_line(prm_fin,key) .and. ret)
        if( trim(key) .eq. 'refpoint' ) then
            ret = prmfile_get_string_value_on_line(prm_fin,mask)
            if( .not. ret ) then
                call pmf_utils_exit(PMF_OUT,1,'refpoint error')
            end if
            write(title,205) i
            call cv_common_read_group_by_mask(cv_item,i+1,title,mask)
            i = i + 1
        end if
        ret = prmfile_next_line(prm_fin)
    end do

! read alpha
    cv_item%alpha = 1.0d0
    if( prmfile_get_real8_by_key(prm_fin,'alpha',cv_item%alpha) ) then
        write(PMF_OUT,220) cv_item%alpha,trim(pmf_unit_label(LengthUnit))
        call pmf_unit_conv_to_ivalue(LengthUnit,cv_item%alpha)
    else
        call pmf_utils_exit(PMF_OUT,1,'alpha is not specified!')
    end if

   100 format('   == Anchor point ===============================')
   200 format('   == Reference points ===========================')
   205 format('refpoint #',I2.2)
   210 format('   ** Num of ref. pts    : ',I6)
   220 format('   ** Alpha              : ',F5.2,' [',A,']')

end subroutine load_paths

!===============================================================================
! Subroutine:  calculate_paths
!===============================================================================

subroutine calculate_paths(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypePATHS)  :: cv_item
    real(PMFDP)         :: x(:,:)
    type(CVContextType) :: ctx
    ! -----------------------------------------------
    integer             :: i,ai,m
    real(PMFDP)         :: d1(3),d2(3),dx(3)
    real(PMFDP)         :: totmass1,totmass2,amass
    real(PMFDP)         :: cu,cd,ce,r2,sc1,sc2,sce,sci
    ! --------------------------------------------------------------------------

! calculate CV value
    totmass1 = 0.0d0
    d1(:) = 0.0
    do  m = 1, cv_item%grps(1)
        ai = cv_item%lindexes(m)
        amass = mass(ai)
        d1(:) = d1(:) + x(:,ai)*amass
        totmass1 = totmass1 + amass
    end do
    if( totmass1 .le. 0 ) then
        call pmf_utils_exit(PMF_OUT,1,'totmass1 is zero in calculate_paths!')
    end if
    d1(:) = d1(:) / totmass1

    cu = 0.0d0
    cd = 0.0d0

    do i=1,cv_item%nrefs
        totmass2 = 0.0d0
        d2(:) = 0.0d0
        do  m = cv_item%grps(i) + 1 , cv_item%grps(i+1)
            ai = cv_item%lindexes(m)
            amass = mass(ai)
            d2(:) = d2(:) + x(:,ai)*amass
            totmass2 = totmass2 + amass
        end do
        if( totmass2 .le. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_paths!')
        end if
        d2(:) = d2(:) / totmass2

        dx(:) = d1(:) - d2(:)
        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dx)
        end if

        r2 = dx(1)**2 + dx(2)**2 + dx(3)**2

        ce = exp(-r2/cv_item%alpha**2)
        cu = cu + real(i-0,PMFDP)*ce
        cd = cd + ce
    end do

    ctx%CVsValues(cv_item%idx) = cu / ( cd *real(cv_item%nrefs-1) )

    ! (a'b - a*b')/b^2
    sc1 = 1.0 / ( cd * real(cv_item%nrefs-1) )    ! cu'
    sc2 = cu / (cd * cd * real(cv_item%nrefs-1) ) ! cd'

! ------------------------------------------------
! calculate derivatives

    cu = 0.0d0

    do i=1,cv_item%nrefs
        totmass2 = 0.0d0
        d2(:) = 0.0d0
        do  m = cv_item%grps(i) + 1 , cv_item%grps(i+1)
            ai = cv_item%lindexes(m)
            amass = mass(ai)
            d2(:) = d2(:) + x(:,ai)*amass
            totmass2 = totmass2 + amass
        end do
        if( totmass2 .le. 0 ) then
            call pmf_utils_exit(PMF_OUT,1,'totmass2 is zero in calculate_pathz!')
        end if
        d2(:) = d2(:) / totmass2

        dx(:) = d1(:) - d2(:)
        if( fenable_pbc ) then
            call pmf_pbc_image_vector(dx)
        end if

        r2 = dx(1)**2 + dx(2)**2 + dx(3)**2

        sce = exp(-r2/cv_item%alpha**2)

        sci = real(i-0,PMFDP)*sce

        do  m = 1, cv_item%grps(1)
            ai = cv_item%lindexes(m)
            amass = mass(ai)
            ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) &
                                          - 2.0d0*(sc1*sci - sc2*sce)*dx(:)*amass/(totmass1 * cv_item%alpha**2)
        end do

        do  m = cv_item%grps(i) + 1 , cv_item%grps(i+1)
            ai = cv_item%lindexes(m)
            amass = mass(ai)
            ctx%CVsDrvs(:,ai,cv_item%idx) = ctx%CVsDrvs(:,ai,cv_item%idx) &
                                          + 2.0d0*(sc1*sci - sc2*sce)*dx(:)*amass/(totmass2 * cv_item%alpha**2)
        end do

    end do

 return

end subroutine calculate_paths

!===============================================================================

end module cv_paths

