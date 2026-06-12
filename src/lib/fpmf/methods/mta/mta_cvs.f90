!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module mtc_cvs

implicit none
contains

!===============================================================================
! Subroutine:  mtc_cvs_reset_cv
!===============================================================================

subroutine mtc_cvs_reset_cv(mtc_item)

    use mtc_dat

    implicit none
    type(CVTypeMTC) :: mtc_item
    ! --------------------------------------------------------------------------

    mtc_item%cvindx         = 0         ! CV index
    mtc_item%cv             => null()

    mtc_item%min_value      = 0.0       ! left range
    mtc_item%max_value      = 0.0       ! right range
    mtc_item%nbins          = 0         ! number of bins

end subroutine mtc_cvs_reset_cv

!===============================================================================
! Subroutine:  mtc_cvs_read_cv
!===============================================================================

subroutine mtc_cvs_read_cv(prm_fin,mtc_item)

    use prmfile
    use mtc_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeMTC)                     :: mtc_item
    ! --------------------------------------------------------------------------

! used CV cannot be controlled by the path subsystem
    if( mtc_item%cv%pathidx .gt. 0 ) then
        if( PathList(mtc_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in a driven mode!')
        end if
    end if

! main CV setup
    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',mtc_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) mtc_item%min_value, trim(mtc_item%cv%get_ulabel())
    call mtc_item%cv%conv_to_ivalue(mtc_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',mtc_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) mtc_item%max_value, trim(mtc_item%cv%get_ulabel())
    call mtc_item%cv%conv_to_ivalue(mtc_item%max_value)

    if( mtc_item%max_value .le. mtc_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',mtc_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    if( mtc_item%nbins .lt. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins has to be greater then zero!')
    end if
    write(PMF_OUT,125) mtc_item%nbins

    return

110 format('    ** Min value         : ',F16.7,' [',A,']')
120 format('    ** Max value         : ',F16.7,' [',A,']')
125 format('    ** Number of bins    : ',I8)

end subroutine mtc_cvs_read_cv

!===============================================================================
! Subroutine:  mtc_cvs_cv_info
!===============================================================================

subroutine mtc_cvs_cv_info(mtc_item)

    use mtc_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit
    use prmfile

    implicit none
    type(CVTypeMTC) :: mtc_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,145) trim(mtc_item%cv%name)
    write(PMF_OUT,146) trim(mtc_item%cv%ctype)
    write(PMF_OUT,150) mtc_item%cv%get_rvalue(CVContext%CVsValues(mtc_item%cvindx)), &
                    trim(mtc_item%cv%get_ulabel())

    write(PMF_OUT,155) mtc_item%cv%get_rvalue(mtc_item%min_value), &
                    trim(mtc_item%cv%get_ulabel())
    write(PMF_OUT,160) mtc_item%cv%get_rvalue(mtc_item%max_value), &
                    trim(mtc_item%cv%get_ulabel())
    write(PMF_OUT,165) mtc_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)

end subroutine mtc_cvs_cv_info

!===============================================================================

end module mtc_cvs
