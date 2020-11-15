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

module abp_cvs

implicit none
contains

!===============================================================================
! Subroutine:  abp_cvs_reset_cv
!===============================================================================

subroutine abp_cvs_reset_cv(abp_item)

    use abp_dat

    implicit none
    type(CVTypeABP)       :: abp_item
    ! --------------------------------------------------------------------------

    abp_item%cvindx          = 0   ! CV index
    abp_item%cv              => null()
    abp_item%set             = 0
    abp_item%min_value       = 0.0 ! left range
    abp_item%max_value       = 0.0 ! right range
    abp_item%nbins           = 0.0 ! number of bins
    abp_item%alpha           = 0.1 ! alpha factor

end subroutine abp_cvs_reset_cv

!===============================================================================
! Subroutine:  abp_cvs_read_cv
!===============================================================================

subroutine abp_cvs_read_cv(prm_fin,abp_item)

    use prmfile
    use abp_dat
    use pmf_cvs
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeABP)                    :: abp_item

    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( abp_item%cv%pathidx .gt. 0 ) then
        if( PathList(abp_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
        end if
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',abp_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) abp_item%min_value, trim(abp_item%cv%get_ulabel())
    call abp_item%cv%conv_to_ivalue(abp_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',abp_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) abp_item%max_value, trim(abp_item%cv%get_ulabel())
    call abp_item%cv%conv_to_ivalue(abp_item%max_value)

    if( abp_item%max_value .le. abp_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'alpha',abp_item%alpha) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: alpha is not specified!')
    end if

    write(PMF_OUT,130) abp_item%alpha, trim(abp_item%cv%get_ulabel())
    call abp_item%cv%conv_to_ivalue(abp_item%alpha)

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',abp_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: nbins is not specified!')
    end if
    write(PMF_OUT,140) abp_item%nbins

    return

110 format('   ** Min value          : ',F16.6,' [',A,']')
120 format('   ** Max value          : ',F16.6,' [',A,']')
130 format('   ** Alpha              : ',F16.6,' [',A,']')
140 format('   ** Number of bins     : ',I9)

end subroutine abp_cvs_read_cv

!===============================================================================
! Subroutine:  abp_cvs_cv_info
!===============================================================================

subroutine abp_cvs_cv_info(abp_item)

    use abp_dat
    use pmf_dat
    use pmf_cvs

    implicit none
    type(CVTypeABP)    :: abp_item
    ! -----------------------------------------------
    integer            :: ci
    ! --------------------------------------------------------------------------

    ci = abp_item%cvindx

    write(PMF_OUT,145) trim(abp_item%cv%name)
    write(PMF_OUT,146) trim(abp_item%cv%ctype)
    write(PMF_OUT,150) abp_item%cv%get_rvalue(CVContext%CVsValues(abp_item%cvindx)), &
                    trim(abp_item%cv%get_ulabel())
    write(PMF_OUT,155) abp_item%cv%get_rvalue(abp_item%min_value), &
                    trim(abp_item%cv%get_ulabel())
    write(PMF_OUT,160) abp_item%cv%get_rvalue(abp_item%max_value), &
                    trim(abp_item%cv%get_ulabel())
    write(PMF_OUT,170) abp_item%cv%get_rvalue(abp_item%alpha), &
                    trim(abp_item%cv%get_ulabel())
    write(PMF_OUT,180) abp_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
170 format('    ** Alpha             : ',E16.7,' [',A,']')
180 format('    ** Number of bins    : ',I8)

end subroutine abp_cvs_cv_info

!===============================================================================

end module abp_cvs
