!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2007 Martin Petrek, petrek@chemi.muni.cz &
!                       Petr Kulhanek, kulhanek@enzim.hu
!    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
!                       Martin Petrek, petrek@chemi.muni.cz
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

module usabf_cvs

implicit none
contains

!===============================================================================
! Subroutine:  usabf_cvs_reset_cv
!===============================================================================

subroutine usabf_cvs_reset_cv(usabf_item)

    use usabf_dat

    implicit none
    type(CVTypeUSABF)   :: usabf_item
    ! --------------------------------------------------------------------------

    usabf_item%cvindx           = 0     ! CV index
    usabf_item%cv               => null()
    usabf_item%set              = 0
    usabf_item%min_value        = 0.0   ! left range
    usabf_item%max_value        = 0.0   ! right range
    usabf_item%nbins            = 0     ! number of bins
    usabf_item%target_value     = 0
    usabf_item%force_constant   = 0
    usabf_item%deviation        = 0

end subroutine usabf_cvs_reset_cv

!===============================================================================
! Subroutine:  usabf_cvs_read_cv
!===============================================================================

subroutine usabf_cvs_read_cv(prm_fin,usabf_item)

    use prmfile
    use usabf_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeUSABF)                   :: usabf_item
    ! -----------------------------------------------
    type(UnitType)                      :: forceunit
    type(UnitType)                      :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, usabf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    ! used CV cannot be controlled by the path subsystem
    if( usabf_item%cv%pathidx .gt. 0 ) then
        if( PathList(usabf_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
        end if
    end if

    if( prmfile_get_real8_by_key(prm_fin,'value',usabf_item%target_value) ) then
        call usabf_item%cv%conv_to_ivalue(usabf_item%target_value)
        write(PMF_OUT,200) usabf_item%cv%get_rvalue(usabf_item%target_value), trim(usabf_item%cv%get_ulabel())
    else
        call pmf_utils_exit(PMF_OUT,1,'Target value is not specified (required by US-ABF)!')
    end if

    ! force constant =========
    if( .not. prmfile_get_real8_by_key(prm_fin,'force_constant',usabf_item%force_constant) ) then
        call pmf_utils_exit(PMF_OUT,1,'force_constant is not specified!')
    end if
    call pmf_unit_conv_to_ivalue(forceunit,usabf_item%force_constant)
    write(PMF_OUT,210) pmf_unit_get_rvalue(forceunit,usabf_item%force_constant), trim(pmf_unit_label(forceunit))

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',usabf_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) usabf_item%min_value, trim(usabf_item%cv%get_ulabel())
    call usabf_item%cv%conv_to_ivalue(usabf_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',usabf_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) usabf_item%max_value, trim(usabf_item%cv%get_ulabel())
    call usabf_item%cv%conv_to_ivalue(usabf_item%max_value)

    if( usabf_item%max_value .le. usabf_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',usabf_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    write(PMF_OUT,125) usabf_item%nbins

    ! ============================================
    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, usabf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    return

110 format('    ** Min value         : ',F16.7,' [',A,']')
120 format('    ** Max value         : ',F16.7,' [',A,']')
125 format('    ** Number of bins    : ',I8)

200 format('    ** Value             : ',F16.7,' [',A,']')
210 format('    ** Force constant    : ',F16.7,' [',A,']')

end subroutine usabf_cvs_read_cv

!===============================================================================
! Subroutine:  usabf_cvs_cv_info
!===============================================================================

subroutine usabf_cvs_cv_info(usabf_item)

    use usabf_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit

    implicit none
    type(CVTypeUSABF)   :: usabf_item
    ! -----------------------------------------------
    type(UnitType)      :: forceunit
    type(UnitType)      :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, usabf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    write(PMF_OUT,145) trim(usabf_item%cv%name)
    write(PMF_OUT,146) trim(usabf_item%cv%ctype)
    write(PMF_OUT,150) usabf_item%cv%get_rvalue(CVContext%CVsValues(usabf_item%cvindx)), &
                    trim(usabf_item%cv%get_ulabel())

    write(PMF_OUT,130) usabf_item%cv%get_rvalue(usabf_item%target_value), &
                        trim(usabf_item%cv%get_ulabel())
    write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,usabf_item%force_constant), &
                        trim(pmf_unit_label(forceunit))

    write(PMF_OUT,155) usabf_item%cv%get_rvalue(usabf_item%min_value), &
                    trim(usabf_item%cv%get_ulabel())
    write(PMF_OUT,160) usabf_item%cv%get_rvalue(usabf_item%max_value), &
                    trim(usabf_item%cv%get_ulabel())
    write(PMF_OUT,165) usabf_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')

130 format('    ** Target value      : ',E16.9,' [',A,']')
180 format('    ** Force constant    : ',E16.9,' [',A,']')

155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)

end subroutine usabf_cvs_cv_info

!===============================================================================

end module usabf_cvs
