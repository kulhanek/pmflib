#ifndef QuantityUnitH
#define QuantityUnitH
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

//! internal PMFLib units are
//! time   - fs
//! length - A
//! angle  - rad
//! energy - kcal/mol

#define MAX_QUANTITIES

// -----------------------------------------------------------------------------

class CUnitType {
public:
    CUnitType(void);

// information methods ---------------------------------------------------------
    CSmallString    GetUnitLabel(void);

// executive methods -----------------------------------------------------------
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

// section of private data -----------------------------------------------------
private:
    double  ConversionFac;          // conversion factor
    int     Units[MAX_QUANTITIES];  //  unit types
    int     Expns[MAX_QUANTITIES];  // unit exponents

private:
    CSmallString    GetQuantityLabel(int quantity);

    CSmallString    GetMassLabel(int unit);
    CSmallString    GetAmoutLabel(int unit);
    CSmallString    GetTimeLabel(int unit);
    CSmallString    GetLengthLabel(int unit);
    CSmallString    GetAngleLabel(int unit);
    CSmallString    GetEnergyLabel(int unit);
    CSmallString    GetTemperatureLabel(int unit);

    void DecodeTimeUnit(const CSmallString& unit_label);
    void DecodeEnergyUnit(const CSmallString& unit_label);
    void DecodeLengthUnit(const CSmallString& unit_label);
    void DecodeAngleUnit(const CSmallString& unit_label);
    void DecodeMassUnit(const CSmallString& unit_label);
    void DecodeTemperatureUnit(const CSmallString& unit_label);
};

// -----------------------------------------------------------------------------

// predefined units for I/O conversions
extern CUnitType TimeUnit;
extern CUnitType EnergyUnit;
extern CUnitType LengthUnit;
extern CUnitType AngleUnit;
extern CUnitType MassUnit;
extern CUnitType TemperatureUnit;

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
                call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unit mismatch in units multiplication!')
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
                call pmf_utils_exit(PMF_OUT,1,'>>> ERROR: Unit mismatch in units division!')
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


