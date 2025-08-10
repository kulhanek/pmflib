// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <EnergyDerProxy.hpp>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyDerProxy::CEnergyDerProxy(void)
{
    Realm = "NONE";
    Description = "NONE";
}

//------------------------------------------------------------------------------

CEnergyDerProxy::~CEnergyDerProxy(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergyDerProxy::Init(CPMFAccumulatorPtr accu)
{
    if( std::find(Requires.begin(),Requires.end(),string(accu->GetMethod())) == Requires.end() ) {
        CSmallString error;
        error << "PMF accumulator '" << accu->GetMethod() << "' is inconsistent with EnergyDerProxy requirements '" << join(Requires,",") << "'";
        RUNTIME_ERROR(error)
    }
    Accu = accu;
}

//------------------------------------------------------------------------------

bool CEnergyDerProxy::IsCompatible(CPMFAccumulatorPtr accu)
{
    return( std::find(Requires.begin(), Requires.end(), std::string(accu->GetMethod())) != Requires.end());
}

//------------------------------------------------------------------------------

bool CEnergyDerProxy::SetType(const CSmallString& realm)
{
    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccumulatorPtr CEnergyDerProxy::GetAccu(void)
{
    return(Accu);
}

//------------------------------------------------------------------------------

CEnergyProxyPtr CEnergyDerProxy::GetEnergyCorrection(void)
{
    return(CEnergyProxyPtr());
}

//------------------------------------------------------------------------------

CSmallString CEnergyDerProxy::GetRealm(void)
{
    return(Realm);
}

//------------------------------------------------------------------------------

CSmallString CEnergyDerProxy::GetDescription(void)
{
    return(Description);
}

//------------------------------------------------------------------------------

int CEnergyDerProxy::GetNumOfCVs(void) const
{
    if( Accu == NULL ) return(0);
    return(Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

int CEnergyDerProxy::GetNumOfBins(void) const
{
    if( Accu == NULL ) return(0);
    return(Accu->GetNumOfBins());
}

//------------------------------------------------------------------------------

int CEnergyDerProxy::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

void CEnergyDerProxy::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CEnergyDerProxy::GetValue(int ibin,int cv,EProxyRealm realm) const
{
    RUNTIME_ERROR("need to be overloaded");
    return(0.0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
