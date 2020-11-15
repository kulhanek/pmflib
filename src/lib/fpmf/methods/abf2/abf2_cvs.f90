!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module abf2_cvs

implicit none
contains

!===============================================================================
! Subroutine:  abf2_cvs_reset_cv
!===============================================================================

subroutine abf2_cvs_reset_cv(abf2_item)

    use abf2_dat

    implicit none
    type(CVTypeABF) :: abf2_item
    ! --------------------------------------------------------------------------

    abf2_item%cvindx         = 0     ! CV index
    abf2_item%cv             => null()
    abf2_item%set            = 0
    abf2_item%min_value      = 0.0   ! left range
    abf2_item%max_value      = 0.0   ! right range
    abf2_item%nbins          = 0     ! number of bins

end subroutine abf2_cvs_reset_cv

!===============================================================================
! Subroutine:  abf2_cvs_read_cv
!===============================================================================

subroutine abf2_cvs_read_cv(prm_fin,abf2_item)

    use prmfile
    use abf2_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeABF)                     :: abf2_item
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( abf2_item%cv%pathidx .gt. 0 ) then
        if( PathList(abf2_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
        end if
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',abf2_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) abf2_item%min_value, trim(abf2_item%cv%get_ulabel())
    call abf2_item%cv%conv_to_ivalue(abf2_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',abf2_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) abf2_item%max_value, trim(abf2_item%cv%get_ulabel())
    call abf2_item%cv%conv_to_ivalue(abf2_item%max_value)

    if( abf2_item%max_value .le. abf2_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',abf2_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    write(PMF_OUT,125) abf2_item%nbins

    return

110 format('   ** Min value          : ',F16.7,' [',A,']')
120 format('   ** Max value          : ',F16.7,' [',A,']')
125 format('   ** Number of bins     : ',I8)

end subroutine abf2_cvs_read_cv

!===============================================================================
! Subroutine:  abf2_cvs_cv_info
!===============================================================================

subroutine abf2_cvs_cv_info(abf2_item)

    use abf2_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit

    implicit none
    type(CVTypeABF) :: abf2_item
    ! -----------------------------------------------
    type(UnitType)  :: forceunit
    type(UnitType)  :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, abf2_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    write(PMF_OUT,145) trim(abf2_item%cv%name)
    write(PMF_OUT,146) trim(abf2_item%cv%ctype)
    write(PMF_OUT,150) abf2_item%cv%get_rvalue(CVContext%CVsValues(abf2_item%cvindx)), &
                    trim(abf2_item%cv%get_ulabel())
    write(PMF_OUT,155) abf2_item%cv%get_rvalue(abf2_item%min_value), &
                    trim(abf2_item%cv%get_ulabel())
    write(PMF_OUT,160) abf2_item%cv%get_rvalue(abf2_item%max_value), &
                    trim(abf2_item%cv%get_ulabel())
    write(PMF_OUT,165) abf2_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)

end subroutine abf2_cvs_cv_info

!===============================================================================

end module abf2_cvs
