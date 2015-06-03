!===============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
!-------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
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

module pmf_unit

use pmf_sizes
use pmf_constants

implicit none

!===============================================================================

! conversion factors -----------------------------------------------------------
! this is used for the conversion of data from external MD codes to
! internal PMFLib units
real(PMFDP) :: MassConv     = 1.0d0
real(PMFDP) :: LengthConv   = 1.0d0
real(PMFDP) :: AngleConv    = 1.0d0
real(PMFDP) :: TimeConv     = 1.0d0
real(PMFDP) :: VelocityConv = 1.0d0
real(PMFDP) :: EnergyConv   = 1.0d0
real(PMFDP) :: ForceConv    = 1.0d0

! internal PMFLib units are
! time   - fs
! length - A
! angle  - rad
! energy - kcal/mol

! quantities -------------------------------------------------------------------
integer,parameter   :: QUANTITY_MASS            = 1     ! mass
integer,parameter   :: QUANTITY_TIME            = 2     ! time
integer,parameter   :: QUANTITY_LENGTH          = 3     ! length
integer,parameter   :: QUANTITY_ANGLE           = 4     ! angle
integer,parameter   :: QUANTITY_ENERGY          = 5     ! energy
integer,parameter   :: QUANTITY_AMOUNT          = 6     ! amount
integer,parameter   :: QUANTITY_TEMPERATURE     = 7     ! temperature
integer,parameter   :: MAX_QUANTITIES   = 7

! unit object ------------------------------------------------------------------
type UnitType
    real(PMFDP) :: ConversionFac            ! conversion factor
    integer     :: Units(MAX_QUANTITIES)    ! unit types
    integer     :: Expns(MAX_QUANTITIES)    ! unit exponents
end type UnitType

type(UnitType)  :: TimeUnit                 ! for input/output conversions
type(UnitType)  :: EnergyUnit
type(UnitType)  :: LengthUnit
type(UnitType)  :: AngleUnit
type(UnitType)  :: MassUnit
type(UnitType)  :: TemperatureUnit

! unit types -------------------------------------------------------------------
integer,parameter   :: AMOUNT_MOL  = 1          ! mol

integer,parameter   :: MASS_AU  = 1             ! au
integer,parameter   :: MASS_G   = 2             ! g

integer,parameter   :: TIME_AU = 1              ! au
integer,parameter   :: TIME_FS = 2              ! fs
integer,parameter   :: TIME_PS = 3              ! ps

integer,parameter   :: LENGTH_AU  = 1           ! au
integer,parameter   :: LENGTH_ANG = 2           ! ang
integer,parameter   :: LENGTH_PM  = 3           ! pm
integer,parameter   :: LENGTH_NM  = 4           ! nm

integer,parameter   :: ANGLE_RAD = 1            ! rad
integer,parameter   :: ANGLE_DEG = 2            ! deg

integer,parameter   :: ENERGY_AU    = 1         ! au
integer,parameter   :: ENERGY_KCAL  = 2         ! kcal
integer,parameter   :: ENERGY_KJ    = 3         ! kJ
integer,parameter   :: ENERGY_eV    = 4         ! eV

integer,parameter   :: TEMPERATURE_K    = 1     ! K

!===============================================================================

contains

!===============================================================================
! Subroutine:   pmf_unit_init
!===============================================================================

subroutine pmf_unit_init(unit)

    implicit none
    type(UnitType)      :: unit
    ! --------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    unit%ConversionFac  = 1.0d0
    do i=1,MAX_QUANTITIES
        unit%Units(i) = 0
        unit%Expns(i) = 0
    end do

end subroutine pmf_unit_init

!===============================================================================
! Function:   pmf_unit_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_label(unit)

    implicit none
    type(UnitType)              :: unit
    ! --------------------------------------------
    integer                     :: i, any
    character(PMF_MAX_SUNIT)    :: sexp
    ! --------------------------------------------------------------------------

    pmf_unit_label = ''
    any = 0

    ! positive exponents
    do i=1,MAX_QUANTITIES
        if( unit%Expns(i) .gt. 0 ) then
            pmf_unit_label = trim(pmf_unit_label) // ' ' // trim(pmf_unit_quantity_label(i,unit%Units(i)))
            if( unit%Expns(i) .gt. 1 ) then
                write(sexp,'(I1)') unit%Expns(i)
                pmf_unit_label = trim(pmf_unit_label) // '^' // trim(sexp)
            end if
            any = any + 1
        end if
    end do

    ! negative exponents
    do i=1,MAX_QUANTITIES
        if( unit%Expns(i) .lt. 0 ) then
            pmf_unit_label = trim(pmf_unit_label) // ' ' // trim(pmf_unit_quantity_label(i,unit%Units(i)))
            write(sexp,'(I2)') unit%Expns(i)
            pmf_unit_label = trim(pmf_unit_label) // '^' // trim(sexp)
            any = any + 1
        end if
    end do

    if( any .eq. 0 ) then
        ! unit less quantity
        pmf_unit_label = '1'
    end if

    ! remove leading space
    pmf_unit_label = adjustl(pmf_unit_label)

end function pmf_unit_label

!===============================================================================
! Function:   pmf_unit_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_quantity_label(quantity,unit)

    use pmf_utils

    implicit none
    integer      :: quantity
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(quantity)
        case(QUANTITY_MASS)
            pmf_unit_quantity_label =  pmf_unit_mass_label(unit)
        case(QUANTITY_AMOUNT)
            pmf_unit_quantity_label =  pmf_unit_amount_label(unit)
        case(QUANTITY_TIME)
            pmf_unit_quantity_label =  pmf_unit_time_label(unit)
        case(QUANTITY_LENGTH)
            pmf_unit_quantity_label =  pmf_unit_length_label(unit)
        case(QUANTITY_ANGLE)
            pmf_unit_quantity_label =  pmf_unit_angle_label(unit)
        case(QUANTITY_ENERGY)
            pmf_unit_quantity_label =  pmf_unit_energy_label(unit)
        case(QUANTITY_TEMPERATURE)
            pmf_unit_quantity_label =  pmf_unit_temperature_label(unit)
        case default
            call pmf_utils_exit(PMF_OUT,1,'Unsupported quantity kind in pmf_unit_quantity_label!')
    end select

end function pmf_unit_quantity_label

!===============================================================================
! Function:   pmf_unit_mass_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_mass_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(MASS_AU)
           pmf_unit_mass_label = 'au'
       case(MASS_G)
           pmf_unit_mass_label = 'g'
       case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported mass unit in pmf_unit_mass_label!')
    end select

end function pmf_unit_mass_label

!===============================================================================
! Function:   pmf_unit_amount_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_amount_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(AMOUNT_MOL)
           pmf_unit_amount_label = 'mol'
       case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported amount unit in pmf_unit_amount_label!')
    end select

end function pmf_unit_amount_label

!===============================================================================
! Function:   pmf_unit_time_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_time_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
        case(TIME_AU)
           pmf_unit_time_label = 'au'
        case(TIME_FS)
           pmf_unit_time_label = 'fs'
        case(TIME_PS)
           pmf_unit_time_label = 'ps'
        case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported time unit in pmf_unit_time_label!')
    end select

end function pmf_unit_time_label

!===============================================================================
! Function:   pmf_unit_length_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_length_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(LENGTH_AU)
           pmf_unit_length_label = 'au'
       case(LENGTH_ANG)
           pmf_unit_length_label = 'A'
       case(LENGTH_PM)
           pmf_unit_length_label = 'pm'
       case(LENGTH_NM)
           pmf_unit_length_label = 'nm'
       case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported lenght unit in pmf_unit_length_label!')
    end select

end function pmf_unit_length_label

!===============================================================================
! Function:   pmf_unit_angle_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_angle_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(ANGLE_RAD)
           pmf_unit_angle_label = 'rad'
       case(ANGLE_DEG)
           pmf_unit_angle_label = 'deg'
       case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported angle unit in pmf_unit_angle_label!')
    end select

end function pmf_unit_angle_label

!===============================================================================
! Function:   pmf_unit_energy_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_energy_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(ENERGY_AU)
           pmf_unit_energy_label = 'au'
       case(ENERGY_KCAL)
           pmf_unit_energy_label = 'kcal'
       case(ENERGY_KJ)
           pmf_unit_energy_label = 'kJ'
       case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported energy unit in pmf_unit_energy_label!')
    end select

end function pmf_unit_energy_label

!===============================================================================
! Function:   pmf_unit_temperature_label
!===============================================================================

character(PMF_MAX_SUNIT) function pmf_unit_temperature_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(TEMPERATURE_K)
           pmf_unit_temperature_label = 'K'
       case default
           call pmf_utils_exit(PMF_OUT,1,'Unsupported temperature unit in pmf_unit_temperature_label!')
    end select

end function pmf_unit_temperature_label

!===============================================================================
! #############################################################################
!===============================================================================

!===============================================================================
! Subroutine:   pmf_unit_decode_time_unit
!===============================================================================

subroutine pmf_unit_decode_time_unit(unit_label,unit)

    use pmf_utils

    implicit none
    character(*)        :: unit_label
    type(UnitType)      :: unit
    ! --------------------------------------------------------------------------

    call pmf_unit_init(unit)

    select case(unit_label)
        case('au','a.u.')
            unit%ConversionFac = PMF_FS2AU
            unit%Units(QUANTITY_TIME) = TIME_AU
            unit%Expns(QUANTITY_TIME) = 1
        case('fs')
            unit%ConversionFac = 1.0d0
            unit%Units(QUANTITY_TIME) = TIME_FS
            unit%Expns(QUANTITY_TIME) = 1
        case('ps')
            unit%ConversionFac = PMF_FS2PS
            unit%Units(QUANTITY_TIME) = TIME_PS
            unit%Expns(QUANTITY_TIME) = 1
        case default
            write(PMF_OUT,'(A)')   ''//trim(unit_label)//' is not supported time unit!'
            write(PMF_OUT,'(A,/)') '           Supported units are: au (a.u.); fs; ps'
            call pmf_utils_exit(PMF_OUT,1,'Unable to decode time unit!')
    end select

end subroutine pmf_unit_decode_time_unit

!===============================================================================
! Subroutine:   pmf_unit_decode_energy_unit
!===============================================================================

subroutine pmf_unit_decode_energy_unit(unit_label,unit)

    use pmf_utils

    implicit none
    character(*)        :: unit_label
    type(UnitType)      :: unit
    ! --------------------------------------------------------------------------

    call pmf_unit_init(unit)

    select case(unit_label)
        case('au','a.u.')
            unit%ConversionFac = PMF_KCL2HARTREE
            unit%Units(QUANTITY_ENERGY) = ENERGY_AU
            unit%Expns(QUANTITY_ENERGY) = 1
        case('eV')
            unit%ConversionFac = PMF_KCL2eV
            unit%Units(QUANTITY_ENERGY) = ENERGY_eV
            unit%Expns(QUANTITY_ENERGY) = 1
        case('kcal/mol')
            unit%ConversionFac = 1.0d0
            unit%Units(QUANTITY_ENERGY) = ENERGY_KCAL
            unit%Expns(QUANTITY_ENERGY) = 1
            unit%Units(QUANTITY_AMOUNT) = AMOUNT_MOL
            unit%Expns(QUANTITY_AMOUNT) = -1
        case('kJ/mol')
            unit%ConversionFac = PMF_KCL2KJ
            unit%Units(QUANTITY_ENERGY) = ENERGY_KJ
            unit%Expns(QUANTITY_ENERGY) = 1
            unit%Units(QUANTITY_AMOUNT) = AMOUNT_MOL
            unit%Expns(QUANTITY_AMOUNT) = -1
        case default
            write(PMF_OUT,'(A)')   ''//trim(unit_label)//' is not supported energy unit!'
            write(PMF_OUT,'(A,/)') '           Supported units are: au (a.u.); eV; kcal/mol; kJ/mol'
            call pmf_utils_exit(PMF_OUT,1,'Unable to decode energy unit!')
    end select

end subroutine pmf_unit_decode_energy_unit

!===============================================================================
! Subroutine:   pmf_unit_decode_length_unit
!===============================================================================

subroutine pmf_unit_decode_length_unit(unit_label,unit)

    use pmf_utils

    implicit none
    character(*)        :: unit_label
    type(UnitType)      :: unit
    ! --------------------------------------------------------------------------

    call pmf_unit_init(unit)

    select case(unit_label)
        case('au','a.u.')
            unit%ConversionFac = PMF_A2AU
            unit%Units(QUANTITY_LENGTH) = LENGTH_AU
            unit%Expns(QUANTITY_LENGTH) = 1
        case('A')
            unit%ConversionFac = 1.0d0
            unit%Units(QUANTITY_LENGTH) = LENGTH_ANG
            unit%Expns(QUANTITY_LENGTH) = 1
        case('pm')
            unit%ConversionFac = PMF_A2PM
            unit%Units(QUANTITY_LENGTH) = LENGTH_PM
            unit%Expns(QUANTITY_LENGTH) = 1
        case('nm')
            unit%ConversionFac = PMF_A2NM
            unit%Units(QUANTITY_LENGTH) = LENGTH_NM
            unit%Expns(QUANTITY_LENGTH) = 1
        case default
            write(PMF_OUT,'(A)')   ''//trim(unit_label)//' is not supported length unit'
            write(PMF_OUT,'(A,/)') '           Supported units are: au (a.u.); A; pm; nm'
            call pmf_utils_exit(PMF_OUT,1,'Unable to decode length unit!')
    end select

end subroutine pmf_unit_decode_length_unit

!===============================================================================
! Subroutine:   pmf_unit_decode_angle_unit
!===============================================================================

subroutine pmf_unit_decode_angle_unit(unit_label,unit)

    use pmf_utils

    implicit none
    character(*)        :: unit_label
    type(UnitType)      :: unit
    ! --------------------------------------------------------------------------

    call pmf_unit_init(unit)

    select case(unit_label)
        case('rad')
            unit%ConversionFac = 1.0d0
            unit%Units(QUANTITY_ANGLE) = ANGLE_RAD
            unit%Expns(QUANTITY_ANGLE) = 1
        case('deg')
            unit%ConversionFac = PMF_R2D
            unit%Units(QUANTITY_ANGLE) = ANGLE_DEG
            unit%Expns(QUANTITY_ANGLE) = 1
        case default
            write(PMF_OUT,'(A)')   ''//trim(unit_label)//' is not supported angle unit'
            write(PMF_OUT,'(A,/)') '           Supported units are: rad; deg'
            call pmf_utils_exit(PMF_OUT,1,'Unable to decode angle unit!')
    end select

end subroutine pmf_unit_decode_angle_unit

!===============================================================================
! Subroutine:   pmf_unit_decode_mass_unit
!===============================================================================

subroutine pmf_unit_decode_mass_unit(unit_label,unit)

    use pmf_utils

    implicit none
    character(*)        :: unit_label
    type(UnitType)      :: unit
    ! --------------------------------------------------------------------------

    call pmf_unit_init(unit)

    select case(unit_label)
        case('au','a.u.')
            unit%ConversionFac = PMF_AMU2AU
            unit%Units(QUANTITY_MASS) = MASS_AU
            unit%Expns(QUANTITY_MASS) = 1
        case('g/mol','amu')
            unit%ConversionFac = 1.0d0
            unit%Units(QUANTITY_MASS) = MASS_G
            unit%Expns(QUANTITY_MASS) = 1
            unit%Units(QUANTITY_AMOUNT) = AMOUNT_MOL
            unit%Expns(QUANTITY_AMOUNT) = -1
        case default
            write(PMF_OUT,'(A)')   ''//trim(unit_label)//' is not supported mass unit'
            write(PMF_OUT,'(A,/)') '           Supported units are: au (a.u.); g/mol (amu)'
            call pmf_utils_exit(PMF_OUT,1,'Unable to decode mass unit!')
    end select

end subroutine pmf_unit_decode_mass_unit

!===============================================================================
! Subroutine:   pmf_unit_decode_temp_unit
!===============================================================================

subroutine pmf_unit_decode_temp_unit(unit_label,unit)

    use pmf_utils

    implicit none
    character(*)        :: unit_label
    type(UnitType)      :: unit
    ! --------------------------------------------------------------------------

    call pmf_unit_init(unit)

    select case(unit_label)
        case('K')
            unit%ConversionFac = 1
            unit%Units(QUANTITY_TEMPERATURE) = TEMPERATURE_K
            unit%Expns(QUANTITY_TEMPERATURE) = 1
        case default
            write(PMF_OUT,'(A)')   ''//trim(unit_label)//' is not supported temperature unit'
            write(PMF_OUT,'(A,/)') '           Supported units are: K'
            call pmf_utils_exit(PMF_OUT,1,'Unable to decode temperature unit!')
    end select

end subroutine pmf_unit_decode_temp_unit

!===============================================================================
! #############################################################################
!===============================================================================

!===============================================================================
! Function:   pmf_unit_get_rvalue
! convert ivalue from internal unit to user specific unit
!===============================================================================

real(PMFDP) function pmf_unit_get_rvalue(unit,ivalue)

    implicit none
    type(UnitType),intent(in)   :: unit
    real(PMFDP),intent(in)      :: ivalue
    ! --------------------------------------------------------------------------

    pmf_unit_get_rvalue = ivalue * unit%ConversionFac

end function pmf_unit_get_rvalue

!===============================================================================
! Function:   pmf_unit_get_ivalue
! convert rvalue from user specific unit to value in internal unit
!===============================================================================

real(PMFDP) function pmf_unit_get_ivalue(unit,rvalue)

    implicit none
    type(UnitType),intent(in)   :: unit
    real(PMFDP),intent(in)      :: rvalue
    ! --------------------------------------------------------------------------

    pmf_unit_get_ivalue = rvalue / unit%ConversionFac

end function pmf_unit_get_ivalue

!===============================================================================
! Function:   pmf_unit_conv_to_ivalue
!===============================================================================

subroutine pmf_unit_conv_to_ivalue(unit,rivalue)

    implicit none
    type(UnitType),intent(in)   :: unit
    real(PMFDP),intent(inout)   :: rivalue
    ! --------------------------------------------------------------------------

    rivalue = pmf_unit_get_ivalue(unit,rivalue)

end subroutine pmf_unit_conv_to_ivalue

!===============================================================================
! #############################################################################
!===============================================================================

!===============================================================================
! Function:   pmf_unit_power_unit
!===============================================================================

type(UnitType) function pmf_unit_power_unit(unit,power)

    implicit none
    type(UnitType),intent(in)   :: unit
    integer,intent(in)          :: power
    ! --------------------------------------------
    integer                     :: i
    ! --------------------------------------------------------------------------

    pmf_unit_power_unit = unit

    pmf_unit_power_unit%ConversionFac = pmf_unit_power_unit%ConversionFac ** power

    do i=1,MAX_QUANTITIES
        pmf_unit_power_unit%Expns(i) = pmf_unit_power_unit%Expns(i) * power
    end do

end function pmf_unit_power_unit

!===============================================================================
! Function:   pmf_unit_mult_units
!===============================================================================

type(UnitType) function pmf_unit_mult_units(unit1,unit2)

    use pmf_utils

    implicit none
    type(UnitType),intent(in)   :: unit1
    type(UnitType),intent(in)   :: unit2
    ! --------------------------------------------
    integer                     :: i
    ! --------------------------------------------------------------------------

    pmf_unit_mult_units = unit1

    pmf_unit_mult_units%ConversionFac = pmf_unit_mult_units%ConversionFac * unit2%ConversionFac

    do i=1,MAX_QUANTITIES
        if( (pmf_unit_mult_units%Units(i) .gt. 0) .and. (unit2%Units(i) .gt. 0) ) then
            if( pmf_unit_mult_units%Units(i) .ne. unit2%Units(i) ) then
                call pmf_utils_exit(PMF_OUT,1,'Unit mismatch in units multiplication!')
            end if
        end if
        if( pmf_unit_mult_units%Units(i) .eq. 0 ) then
            pmf_unit_mult_units%Units(i) = unit2%Units(i)
        end if
        pmf_unit_mult_units%Expns(i) = pmf_unit_mult_units%Expns(i) + unit2%Expns(i)
    end do

end function pmf_unit_mult_units

!===============================================================================
! Function:   pmf_unit_div_units
!===============================================================================

type(UnitType) function pmf_unit_div_units(unit1,unit2)

    use pmf_utils

    implicit none
    type(UnitType),intent(in)   :: unit1
    type(UnitType),intent(in)   :: unit2
    ! --------------------------------------------
    integer                     :: i
    ! --------------------------------------------------------------------------

    pmf_unit_div_units = unit1

    pmf_unit_div_units%ConversionFac = pmf_unit_div_units%ConversionFac / unit2%ConversionFac

    do i=1,MAX_QUANTITIES
        if( (pmf_unit_div_units%Units(i) .gt. 0) .and. (unit2%Units(i) .gt. 0) ) then
            if( pmf_unit_div_units%Units(i) .ne. unit2%Units(i) ) then
                call pmf_utils_exit(PMF_OUT,1,'Unit mismatch in units division!')
            end if
        end if
        if( pmf_unit_div_units%Units(i) .eq. 0 ) then
            pmf_unit_div_units%Units(i) = unit2%Units(i)
        end if
        pmf_unit_div_units%Expns(i) = pmf_unit_div_units%Expns(i) - unit2%Expns(i)
    end do

end function pmf_unit_div_units

!===============================================================================
! Function:   pmf_unit_compare_units
!===============================================================================

logical function pmf_unit_compare_units(unit1,unit2)

    implicit none
    type(UnitType),intent(in)   :: unit1
    type(UnitType),intent(in)   :: unit2
    ! --------------------------------------------
    integer                     :: i
    ! --------------------------------------------------------------------------

    pmf_unit_compare_units = .true.

    if( unit1%ConversionFac .ne. unit2%ConversionFac ) then
        pmf_unit_compare_units = .false.
        return
    end if

    do i=1,MAX_QUANTITIES
        if( unit1%Units(i) .ne. unit2%Units(i) ) then
            pmf_unit_compare_units = .false.
            return
        end if
        if( unit1%Expns(i) .ne. unit2%Expns(i) ) then
            pmf_unit_compare_units = .false.
            return
        end if
    end do

end function pmf_unit_compare_units

!===============================================================================

end module pmf_unit


