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

module mta_cvs

implicit none
contains

!===============================================================================
! Subroutine:  mta_cvs_reset_cv
!===============================================================================

subroutine mta_cvs_reset_cv(mta_item)

    use mta_dat

    implicit none
    type(CVTypeMTA) :: mta_item
    ! --------------------------------------------------------------------------

    mta_item%cvindx         = 0         ! CV index
    mta_item%cv             => null()

    mta_item%min_value      = 0.0d0     ! left range
    mta_item%max_value      = 0.0d0     ! right range
    mta_item%nbins          = 0         ! number of bins

end subroutine mta_cvs_reset_cv

!===============================================================================
! Subroutine:  mta_cvs_read_cv
!===============================================================================

subroutine mta_cvs_read_cv(prm_fin,mta_item)

    use prmfile
    use mta_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths
    use pmf_control_utils

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeMTA)                     :: mta_item
    ! --------------------------------------------------------------------------

! used CV cannot be controlled by the path subsystem
    if( mta_item%cv%pathidx .gt. 0 ) then
        if( PathList(mta_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in a driven mode!')
        end if
    end if

! main CV setup
    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',mta_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) mta_item%min_value, trim(mta_item%cv%get_ulabel())
    call mta_item%cv%conv_to_ivalue(mta_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',mta_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) mta_item%max_value, trim(mta_item%cv%get_ulabel())
    call mta_item%cv%conv_to_ivalue(mta_item%max_value)

    if( mta_item%max_value .le. mta_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater than min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',mta_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    if( mta_item%nbins .lt. 1 ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins has to be greater than zero!')
    end if
    write(PMF_OUT,125) mta_item%nbins

    return

110 format('    ** Min value         : ',F16.7,' [',A,']')
120 format('    ** Max value         : ',F16.7,' [',A,']')
125 format('    ** Number of bins    : ',I8)

end subroutine mta_cvs_read_cv

!===============================================================================
! Subroutine:  mta_cvs_cv_info
!===============================================================================

subroutine mta_cvs_cv_info(mta_item)

    use mta_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit
    use prmfile

    implicit none
    type(CVTypeMTA) :: mta_item
    ! --------------------------------------------------------------------------

    write(PMF_OUT,145) trim(mta_item%cv%name)
    write(PMF_OUT,146) trim(mta_item%cv%ctype)
    write(PMF_OUT,150) mta_item%cv%get_rvalue(CVContext%CVsValues(mta_item%cvindx)), &
                    trim(mta_item%cv%get_ulabel())

    write(PMF_OUT,155) mta_item%cv%get_rvalue(mta_item%min_value), &
                    trim(mta_item%cv%get_ulabel())
    write(PMF_OUT,160) mta_item%cv%get_rvalue(mta_item%max_value), &
                    trim(mta_item%cv%get_ulabel())
    write(PMF_OUT,165) mta_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)

end subroutine mta_cvs_cv_info

!===============================================================================

end module mta_cvs
