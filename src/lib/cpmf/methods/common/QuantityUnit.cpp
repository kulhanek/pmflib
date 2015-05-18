// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// =============================================================================

#include <QuantityUnit.hpp>
#include <ErrorSystem.hpp>

// unit types ------------------------------------------------------------------

#define AMOUNT_MOL  = 1         // mol

#define MASS_AU     = 1         // au
#define MASS_G      = 2         // g

#define TIME_AU     = 1         // au
#define TIME_FS     = 2         // fs
#define TIME_PS     = 3         // ps

#define LENGTH_AU   = 1         // au
#define LENGTH_ANG  = 2         // ang
#define LENGTH_PM   = 3         // pm
#define LENGTH_NM   = 4         // nm

#define ANGLE_RAD   = 1         // rad
#define ANGLE_DEG   = 2         // deg

#define ENERGY_AU   = 1         // au
#define ENERGY_KCAL = 2         // kcal
#define ENERGY_KJ   = 3         // kJ
#define ENERGY_eV   = 4         // eV

#define TEMPERATURE_K   = 1     // K

//------------------------------------------------------------------------------

CUnitType TimeUnit;
CUnitType EnergyUnit;
CUnitType LengthUnit;
CUnitType AngleUnit;
CUnitType MassUnit;
CUnitType TemperatureUnit;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CUnitType::CUnitType(void)
{
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
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSmallString CUnitType::GetUnitLabel(void)
{
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
}

//------------------------------------------------------------------------------

    /// convert from internal to real value
    double GetRValue(double ivalue);

    /// convert from real to internal value
    double GetIValue(double rvalue);

// operators -------------------------------------------------------------------
    /// assignement
    CUnitType& operator = (const CUnitType& right);
    /// multiplication
    const CUnitType operator * (const CUnitType& right) const;
    /// division
    const CUnitType operator / (const CUnitType& right) const;
    /// power
    const CUnitType operator ^ (double power) const;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSmallString CUnitType::GetQuantityLabel(int quantity)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported quantity kind in pmf_unit_quantity_label!')
    end select

end function pmf_unit_quantity_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetMassLabel(int unit)
{
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
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported mass unit in pmf_unit_mass_label!')
    end select

end function pmf_unit_mass_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetAmoutLabel(int unit)
{
character(PMF_MAX_SUNIT) function pmf_unit_amount_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(AMOUNT_MOL)
           pmf_unit_amount_label = 'mol'
       case default
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported amount unit in pmf_unit_amount_label!')
    end select

end function pmf_unit_amount_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetTimeLabel(int unit)
{
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
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported time unit in pmf_unit_time_label!')
    end select

end function pmf_unit_time_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetLengthLabel(int unit)
{
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
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported lenght unit in pmf_unit_length_label!')
    end select

end function pmf_unit_length_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetAngleLabel(int unit)
{
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
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported angle unit in pmf_unit_angle_label!')
    end select

end function pmf_unit_angle_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetEnergyLabel(int unit)
{
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
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported energy unit in pmf_unit_energy_label!')
    end select

end function pmf_unit_energy_label
}

//------------------------------------------------------------------------------

CSmallString CUnitType::GetTemperatureLabel(int unit)
{
character(PMF_MAX_SUNIT) function pmf_unit_temperature_label(unit)

    use pmf_utils

    implicit none
    integer      :: unit
    ! --------------------------------------------------------------------------

    select case(unit)
       case(TEMPERATURE_K)
           pmf_unit_temperature_label = 'K'
       case default
           call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unsupported temperature unit in pmf_unit_temperature_label!')
    end select

end function pmf_unit_temperature_label
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CUnitType::DecodeTimeUnit(const CSmallString& unit_label)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unable to decode time unit!')
    end select

end subroutine pmf_unit_decode_time_unit
}

//------------------------------------------------------------------------------

void CUnitType::DecodeEnergyUnit(const CSmallString& unit_label)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unable to decode energy unit!')
    end select

end subroutine pmf_unit_decode_energy_unit
}

//------------------------------------------------------------------------------

void CUnitType::DecodeLengthUnit(const CSmallString& unit_label)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unable to decode length unit!')
    end select

end subroutine pmf_unit_decode_length_unit
}

//------------------------------------------------------------------------------

void CUnitType::DecodeAngleUnit(const CSmallString& unit_label)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unable to decode angle unit!')
    end select

end subroutine pmf_unit_decode_angle_unit
}

//------------------------------------------------------------------------------

void CUnitType::DecodeMassUnit(const CSmallString& unit_label)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unable to decode mass unit!')
    end select

end subroutine pmf_unit_decode_mass_unit
}

//------------------------------------------------------------------------------

void CUnitType::DecodeTemperatureUnit(const CSmallString& unit_label)
{
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
            call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unable to decode temperature unit!')
    end select

end subroutine pmf_unit_decode_temp_unit
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


