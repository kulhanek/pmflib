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

module tabf_cvs

implicit none
contains

!===============================================================================
! Subroutine:  tabf_cvs_reset_cv
!===============================================================================

subroutine tabf_cvs_reset_cv(tabf_item)

    use tabf_dat

    implicit none
    type(CVTypeABF) :: tabf_item
    ! --------------------------------------------------------------------------

    tabf_item%cvindx         = 0     ! CV index
    tabf_item%cv             => null()
    tabf_item%set            = 0
    tabf_item%min_value      = 0.0   ! left range
    tabf_item%max_value      = 0.0   ! right range
    tabf_item%nbins          = 0     ! number of bins

end subroutine tabf_cvs_reset_cv

!===============================================================================
! Subroutine:  tabf_cvs_read_cv
!===============================================================================

subroutine tabf_cvs_read_cv(prm_fin,tabf_item)

    use prmfile
    use tabf_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)    :: prm_fin
    type(CVTypeABF)                     :: tabf_item
    ! -----------------------------------------------
    type(UnitType)                      :: forceunit
    type(UnitType)                      :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( tabf_item%cv%pathidx .gt. 0 ) then
        if( PathList(tabf_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
        end if
    end if

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'min_value',tabf_item%min_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'min_value is not specified!')
    end if
    write(PMF_OUT,110) tabf_item%min_value, trim(tabf_item%cv%get_ulabel())
    call tabf_item%cv%conv_to_ivalue(tabf_item%min_value)

    ! ========================
    if( .not. prmfile_get_real8_by_key(prm_fin,'max_value',tabf_item%max_value) ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value is not specified!')
    end if
    write(PMF_OUT,120) tabf_item%max_value, trim(tabf_item%cv%get_ulabel())
    call tabf_item%cv%conv_to_ivalue(tabf_item%max_value)

    if( tabf_item%max_value .le. tabf_item%min_value ) then
        call pmf_utils_exit(PMF_OUT,1,'max_value has to be greater then min_value!')
    end if

    ! ========================
    if( .not. prmfile_get_integer_by_key(prm_fin,'nbins',tabf_item%nbins) ) then
        call pmf_utils_exit(PMF_OUT,1,'nbins is not specified!')
    end if
    write(PMF_OUT,125) tabf_item%nbins

    ! ============================================
    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, tabf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    return

110 format('   ** Min value          : ',F16.7,' [',A,']')
120 format('   ** Max value          : ',F16.7,' [',A,']')
125 format('   ** Number of bins     : ',I8)

end subroutine tabf_cvs_read_cv

!===============================================================================
! Subroutine:  tabf_cvs_cv_info
!===============================================================================

subroutine tabf_cvs_cv_info(tabf_item)

    use tabf_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit

    implicit none
    type(CVTypeABF) :: tabf_item
    ! -----------------------------------------------
    type(UnitType)  :: forceunit
    type(UnitType)  :: forceunit_recip
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit       = pmf_unit_div_units( EnergyUnit, tabf_item%cv%unit )
    forceunit_recip = pmf_unit_power_unit( forceunit, -1 )

    write(PMF_OUT,145) trim(tabf_item%cv%name)
    write(PMF_OUT,146) trim(tabf_item%cv%ctype)
    write(PMF_OUT,150) tabf_item%cv%get_rvalue(CVContext%CVsValues(tabf_item%cvindx)), &
                    trim(tabf_item%cv%get_ulabel())
    write(PMF_OUT,155) tabf_item%cv%get_rvalue(tabf_item%min_value), &
                    trim(tabf_item%cv%get_ulabel())
    write(PMF_OUT,160) tabf_item%cv%get_rvalue(tabf_item%max_value), &
                    trim(tabf_item%cv%get_ulabel())
    write(PMF_OUT,165) tabf_item%nbins

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
155 format('    ** Min value         : ',E16.7,' [',A,']')
160 format('    ** Max value         : ',E16.7,' [',A,']')
165 format('    ** Number of bins    : ',I9)

end subroutine tabf_cvs_cv_info

!===============================================================================

end module tabf_cvs
