!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module stm_cvs

implicit none
contains

!===============================================================================
! Subroutine:  stm_cvs_reset_cv
!===============================================================================

subroutine stm_cvs_reset_cv(stm_item)

    use stm_dat

    implicit none
    type(CVTypeSTM)       :: stm_item
    ! --------------------------------------------------------------------------

    stm_item%cvindx          = 0   ! CV index
    stm_item%force_constant  = 50.0 ! force constant

end subroutine stm_cvs_reset_cv

!===============================================================================
! Subroutine:  stm_cvs_read_cv
!===============================================================================

subroutine stm_cvs_read_cv(prm_fin,stm_item)

    use prmfile
    use stm_dat
    use pmf_cvs
    use pmf_unit
    use pmf_utils
    use pmf_paths

    implicit none
    type(PRMFILE_TYPE),intent(inout)   :: prm_fin
    type(CVTypeSTM)                    :: stm_item
    ! -----------------------------------------------
    type(UnitType)                     :: forceunit
    ! --------------------------------------------------------------------------

    ! used CV cannot be controlled by the path subsystem
    if( stm_item%cv%pathidx .gt. 0 ) then
        if( PathList(stm_item%cv%pathidx)%path%driven_mode ) then
            call pmf_utils_exit(PMF_OUT,1,'Requested CV is connected with the path that is in driven mode!')
        end if
    end if

    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, pmf_unit_power_unit(stm_item%cv%unit,2) )

    ! force_constant ==============================================================
    if( .not. prmfile_get_real8_by_key(prm_fin,'force_constant',stm_item%force_constant) ) then
        call pmf_utils_exit(PMF_OUT,1,'force_constant is not specified!')
    else
        write(PMF_OUT,110) stm_item%force_constant, trim(pmf_unit_label(forceunit))
        call pmf_unit_conv_to_ivalue(forceunit,stm_item%force_constant)
    end if

    return

    110 format('   ** Force constant :',F16.6,' [',A,']')

end subroutine stm_cvs_read_cv

!===============================================================================
! Subroutine:  stm_cvs_cv_info
!===============================================================================

subroutine stm_cvs_cv_info(stm_item)

    use stm_dat
    use pmf_dat
    use pmf_cvs
    use pmf_unit

    implicit none
    type(CVTypeSTM)     :: stm_item
    ! -----------------------------------------------
    type(UnitType)      :: forceunit
    ! --------------------------------------------------------------------------

    ! prepare force unit
    forceunit = pmf_unit_div_units( EnergyUnit, pmf_unit_power_unit(stm_item%cv%unit,2) )

    write(PMF_OUT,145) trim(stm_item%cv%name)
    write(PMF_OUT,146) trim(stm_item%cv%ctype)
    write(PMF_OUT,150) stm_item%cv%get_rvalue(CVContext%CVsValues(stm_item%cvindx)), &
                    trim(stm_item%cv%get_ulabel())
    write(PMF_OUT,180) pmf_unit_get_rvalue(forceunit,stm_item%force_constant), &
                    trim(pmf_unit_label(forceunit))

    return

145 format('    ** Name              : ',a)
146 format('    ** Type              : ',a)
150 format('    ** Current value     : ',E16.7,' [',A,']')
180 format('    ** Force constant   : ',E16.9,' [',A,']')

end subroutine stm_cvs_cv_info

!===============================================================================
! Subroutine:  stm_restraints_increment
!===============================================================================

subroutine stm_restraints_increment(stm_item)

    use stm_dat
    use pmf_dat

    implicit none
    type(CVTypeSTM)     :: stm_item
    ! --------------------------------------------------------------------------

    if( stmmode .eq. BMO_ACCUMULATION ) return  ! no change of constant restraint
    if( stmmode .eq. BMO_PRODUCTION ) return  ! no change of constant restraint
    if( stmmode .eq. BMO_TERMINATE ) return  ! no change of constant restraint

    ! incrementation is in linear mode
    stm_item%target_value = stm_item%startvalue + &
                  (stm_item%stopvalue-stm_item%startvalue) * curstep / stmsteps

    return

end subroutine stm_restraints_increment

!===============================================================================

end module stm_cvs
