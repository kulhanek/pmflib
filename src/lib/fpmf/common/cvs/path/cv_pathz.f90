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

! PTPATHZ - path-distance coordinate for point-based paths

module cv_pathz

use pmf_constants
use pmf_dat
use cv_common
use pmf_paths
use pmf_sizes

implicit none

!===============================================================================

type, extends(CVType) :: CVTypePATHZ

    integer                 :: pathindx         ! path index
    class(PathType),pointer :: path             ! path data
    real(PMFDP)             :: alpha            ! unit less

    contains
        procedure :: load_cv        => load_pathz
        procedure :: calculate_cv   => calculate_pathz
end type CVTypePATHZ

!===============================================================================

contains

!===============================================================================
! Subroutine:  load_pathz
!===============================================================================

subroutine load_pathz(cv_item,prm_fin)

    use prmfile
    use pmf_utils

    implicit none
    class(CVTypePATHZ)                  :: cv_item
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    ! --------------------------------------------
    character(len=PRMFILE_MAX_LINE)     :: pathname
    ! --------------------------------------------------------------------------

! simple init and allocation --------------------
    cv_item%ctype         = 'PATHS'
    call pmf_unit_init(cv_item%unit)
    cv_item%gradforanycrd = .true.
    call cv_common_read_name(cv_item,prm_fin)

! determine number of reference points ----------
    if( .not. prmfile_get_string_by_key(prm_fin,'path',pathname)) then
        call pmf_utils_exit(PMF_OUT,1,'[PATHS] The PATH name (path) is not provided!')
    end if

    write(PMF_OUT,10) trim(pathname)

    cv_item%pathindx = pmf_paths_find_path(pathname)
    cv_item%path => PathList(cv_item%pathindx)%path

! read alpha
    cv_item%alpha = 0.0d0
    if( prmfile_get_real8_by_key(prm_fin,'alpha',cv_item%alpha) ) then
        write(PMF_OUT,20) cv_item%alpha,trim(pmf_unit_label(LengthUnit))
    end if

    if( cv_item%alpha .eq. 0.0d0 ) then
        ! calculate
        write(PMF_OUT,20) cv_item%alpha,trim(pmf_unit_label(LengthUnit))
    end if

10 format('   ** Path               : ',A)
20 format('   ** Alpha              : ',F5.3)

end subroutine load_pathz

!===============================================================================
! Subroutine:  calculate_paths
!===============================================================================

subroutine calculate_pathz(cv_item,x,ctx)

    use pmf_dat
    use pmf_pbc
    use pmf_utils

    implicit none
    class(CVTypePATHZ)      :: cv_item
    real(PMFDP)             :: x(:,:)
    type(CVContextType)     :: ctx
    ! -----------------------------------------------
    ! --------------------------------------------------------------------------



 return

end subroutine calculate_pathz

!===============================================================================

end module cv_pathz

