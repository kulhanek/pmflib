// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//                       Martin Petrek, petrek@chemi.muni.cz
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

#include <math.h>
#include <string.h>
#include <ColVariable.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CColVariable::CColVariable(void)
{
    ID          = 0;
    NBins       = 0;
    MinValue    = 0.0;
    MaxValue    = 0.0;
    BinWidth    = 0.0;
    Width       = 0.0;
    MaxMovement = 0.0;
    FConv       = 1.0;
    Alpha       = 0.0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CColVariable::SetName(const CSmallString& lname)
{
    Name = lname;
}

//------------------------------------------------------------------------------

void CColVariable::SetType(const CSmallString& ltype)
{
    Type = ltype;
}

//------------------------------------------------------------------------------

void CColVariable::SetMinValue(double lmin)
{
    MinValue = lmin;

    Width = MaxValue - MinValue;
    BinWidth = 0.0;
    if(NBins > 0){
        BinWidth = Width / ((double)NBins);
    } else {
        NBins = -1;
    }

}

//------------------------------------------------------------------------------

void CColVariable::SetMaxValue(double lmax)
{
    MaxValue = lmax;

    Width = MaxValue - MinValue;
    BinWidth = 0.0;

    if(NBins > 0){
        BinWidth = Width / ((double)NBins);
    } else {
        NBins = -1;
    }
}

//------------------------------------------------------------------------------

void CColVariable::SetMaxMovement(double lmax)
{
    MaxMovement = lmax;
}

//------------------------------------------------------------------------------

void CColVariable::SetCV(int lid,const CSmallString& lname,const CSmallString& ltype,
                            double lmin_value,double lmax_value,unsigned int lnbins)
{
    ID = lid;
    Name = lname;
    Type = ltype;
    MinValue = lmin_value;
    MaxValue = lmax_value;
    NBins = lnbins;

    Width = MaxValue - MinValue;
    BinWidth = 0.0;
    if(NBins > 0) BinWidth = Width / ((double)NBins);
}

//------------------------------------------------------------------------------

void CColVariable::SetCV(int lid,const CSmallString& lname,const CSmallString& ltype)
{
    ID = lid;
    Name = lname;
    Type = ltype;
    MinValue = 0.0;
    MaxValue = 0.0;
    NBins = 0;
    Width = 0.0;
    BinWidth = 0.0;
}

//------------------------------------------------------------------------------

void CColVariable::SetCV(int lid,const CSmallString& lname,const CSmallString& ltype,
                  double lmin_value,double lmax_value,unsigned int lnbins,
                  double lfconv,const CSmallString& lunit)
{
    ID = lid;
    Name = lname;
    Type = ltype;
    MinValue = lmin_value;
    MaxValue = lmax_value;
    NBins = lnbins;

    Width = MaxValue - MinValue;
    BinWidth = 0.0;
    if(NBins > 0) BinWidth = Width / ((double)NBins);

    FConv = lfconv;
    Unit = lunit;
}

//------------------------------------------------------------------------------

void CColVariable::SetCV(int lid,const CSmallString& lname,const CSmallString& ltype,
              double lmin_value,double lmax_value,unsigned int lnbins,
              double alpha,
              double lfconv,const CSmallString& lunit)
{
    ID = lid;
    Name = lname;
    Type = ltype;
    MinValue = lmin_value;
    MaxValue = lmax_value;
    NBins = lnbins;

    Width = MaxValue - MinValue;
    BinWidth = 0.0;
    if(NBins > 0) BinWidth = Width / ((double)NBins);

    Alpha = alpha;
    FConv = lfconv;
    Unit = lunit;
}

//------------------------------------------------------------------------------

// for energy surface
void CColVariable::CopyFrom(const CColVariable* p_coord)
{
    ID        = p_coord->ID;
    Name      = p_coord->Name;
    Type      = p_coord->Type;
    Unit      = p_coord->Unit;
    NBins     = p_coord->NBins;
    FConv     = p_coord->FConv;
    MinValue  = p_coord->MinValue;
    MaxValue  = p_coord->MaxValue;
    BinWidth  = p_coord->BinWidth;
    Width     = p_coord->Width;
    Alpha     = p_coord->Alpha;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CColVariable::GetNumOfBins(void) const
{
    return(NBins);
}

//------------------------------------------------------------------------------

double CColVariable::GetMinValue(void) const
{
    return(MinValue);
}

//------------------------------------------------------------------------------

double CColVariable::GetMaxValue(void) const
{
    return(MaxValue);
}

//------------------------------------------------------------------------------

double CColVariable::GetMaxMovement(void) const
{
    return(MaxMovement);
}

//------------------------------------------------------------------------------

double CColVariable::GetRange(void) const
{
    return(Width);
}

//------------------------------------------------------------------------------

double CColVariable::GetBinWidth(void) const
{
    return(BinWidth);
}

//------------------------------------------------------------------------------

double CColVariable::GetValue(unsigned int bin) const
{
// bin  = 0,...,NBins-1
    return(MinValue + ((double)bin + 0.5)*(MaxValue-MinValue)/((double)NBins));
}

//------------------------------------------------------------------------------

double CColVariable::GetRValue(unsigned int bin) const
{
  return( GetValue(bin)*FConv );
}

//------------------------------------------------------------------------------

const CSmallString& CColVariable::GetType(void) const
{
    return(Type);
}

//------------------------------------------------------------------------------

const CSmallString& CColVariable::GetName(void) const
{
    return(Name);
}

//------------------------------------------------------------------------------

const CSmallString& CColVariable::GetUnit(void) const
{
    return(Unit);
}

//------------------------------------------------------------------------------

double CColVariable::GetDifference(double left,double right) const
{
    // FIXME
    return(left-right);
}

//------------------------------------------------------------------------------

bool CColVariable::IsPeriodic(void) const
{
    if(strstr(Type,"DIH") == NULL) return(false);
    // FIXME - TODO !!!
// currently only dihedral angle can be periodic
    if(fabs(MaxValue - MinValue - 2*M_PI) <= 0.1) {
        return(true);
    }
    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CColVariable::PrintInfo(std::ostream& vout)
{
//    vout << "ID P Type       Unit  Name                       Min value   Max value   NBins  " << endl;
//    vout << "-- - ---------- ----- -------------------------- ----------- ----------- -------" << endl;

    char periodic = 'F';
    if( IsPeriodic() ) periodic = 'T';

    vout << setw(2) << ID+1 << " " << setw(1) << periodic << " " << left << setw(10) << Type << " " << left << setw(5) << Unit << " " << setw(26) << Name << " ";
    vout << right << fixed << setw(11) << setprecision(4) << MinValue*FConv << " ";
    vout << setw(11) << setprecision(4) << MaxValue*FConv << " ";
    vout << setw(7) << NBins << endl;
}

//------------------------------------------------------------------------------

void CColVariable::LoadInfo(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    bool result = true;

    result &= p_ele->GetAttribute("ID",ID);
    result &= p_ele->GetAttribute("Name",Name);
    result &= p_ele->GetAttribute("Type",Type);
    result &= p_ele->GetAttribute("Unit",Unit);
    result &= p_ele->GetAttribute("FConv",FConv);

    if( p_ele->GetAttribute("NBins",NBins) ){
        result &= p_ele->GetAttribute("MinValue",MinValue);
        result &= p_ele->GetAttribute("MaxValue",MaxValue);
    }

    // optional
    p_ele->GetAttribute("MaxMovement",MaxMovement);
    p_ele->GetAttribute("Alpha",Alpha);

    if(result == false) {
        LOGIC_ERROR("unable to get attribute(s)");
    }

    Width = MaxValue - MinValue;
    BinWidth = 0.0;
    if(NBins > 0) BinWidth = Width/((double)NBins);
}

//------------------------------------------------------------------------------

bool CColVariable::CheckInfo(CXMLElement* p_ele) const
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    bool result = true;

    int             lid;
    int             lnbins;
    CSmallString    ltype,lname,lunit;
    double          lmin_value,lmax_value,lunit_fac;

    result &= p_ele->GetAttribute("ID",lid);
    result &= p_ele->GetAttribute("Name",lname);
    result &= p_ele->GetAttribute("Type",ltype);
    result &= p_ele->GetAttribute("Unit",lunit);
    result &= p_ele->GetAttribute("FConv",lunit_fac);

    if( NBins > 0 ){
        result &= p_ele->GetAttribute("NBins",lnbins);
        result &= p_ele->GetAttribute("MinValue",lmin_value);
        result &= p_ele->GetAttribute("MaxValue",lmax_value);
    }

    if(result == false) {
        ES_ERROR("unable to get attribute(s)");
        return(false);
    }

    if(lid != ID) {
        ES_ERROR("mismatch in ID");
        return(false);
    }

    if(lname != Name) {
        ES_ERROR("mismatch in Name");
        return(false);
    }

    if(ltype != Type) {
        ES_ERROR("mismatch in Type");
        return(false);
    }

// DO NOT CHECK - we work in internal units
//    if(lunit != Unit) {
//        ES_ERROR("mismatch in Unit");
//        return(false);
//    }
//
//    if(fabs(lunit_fac-FConv) > fabs(FConv/100000.0)) {
//        ES_ERROR("mismatch in FConv");
//        return(false);
//    }

    if( NBins > 0 ){
        if(fabs(lmin_value-MinValue) > fabs(MinValue/100000.0)) {
            ES_ERROR("mismatch in MinValue");
            return(false);
        }

        if(fabs(lmax_value-MaxValue) > fabs(MaxValue/100000.0)) {
            ES_ERROR("mismatch in MaxValue");
            return(false);
        }

        if(lnbins != NBins) {
            ES_ERROR("mismatch in NBins");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CColVariable::CheckInfo(const CColVariable* p_coord) const
{
    if(p_coord == NULL) {
        INVALID_ARGUMENT("p_coord is NULL");
    }

    if(p_coord->ID != ID) {
        ES_ERROR("mismatch in ID");
        return(false);
    }

    if(p_coord->Name != Name) {
        ES_ERROR("mismatch in Name");
        return(false);
    }

    if(p_coord->Type != Type) {
        ES_ERROR("mismatch in Type");
        return(false);
    }

// DO NOT CHECK - we work in internal units
//    if(p_coord->Unit != Unit) {
//        ES_ERROR("mismatch in Unit");
//        return(false);
//    }
//
//    if(fabs(p_coord->FConv - FConv) > fabs(FConv/100000.0)) {
//        ES_ERROR("mismatch in FConv");
//        return(false);
//    }

    if( NBins > 0 ){

        if(fabs(p_coord->MinValue-MinValue) > fabs(MinValue/100000.0)) {
            ES_ERROR("mismatch in MinValue");
            return(false);
        }

        if(fabs(p_coord->MaxValue-MaxValue) > fabs(MaxValue/100000.0)) {
            ES_ERROR("mismatch in MaxValue");
            return(false);
        }

        if(p_coord->NBins != NBins) {
            ES_ERROR("mismatch in NBins");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CColVariable::SaveInfo(CXMLElement* p_ele) const
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    p_ele->SetAttribute("ID",ID);
    p_ele->SetAttribute("Name",Name);
    p_ele->SetAttribute("Type",Type);
    p_ele->SetAttribute("Unit",Unit);
    p_ele->SetAttribute("FConv",FConv);
    if( NBins != 0 ){
        p_ele->SetAttribute("MinValue",MinValue);
        p_ele->SetAttribute("MaxValue",MaxValue);
        p_ele->SetAttribute("NBins",NBins);
    }
    // optional
    p_ele->SetAttribute("MaxMovement",MaxMovement);
    p_ele->SetAttribute("Alpha",Alpha);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

