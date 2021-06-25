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

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyDerProxy::CEnergyDerProxy(void)
{
    Method = "NONE";
}

//------------------------------------------------------------------------------

CEnergyDerProxy::~CEnergyDerProxy(void)
{

}

//------------------------------------------------------------------------------

int CEnergyDerProxy::GetNumOfSamples(int ibin) const
{
    RUNTIME_ERROR("need to be overloaded");
    return(0);
}

//------------------------------------------------------------------------------

void CEnergyDerProxy::SetNumOfSamples(int ibin,int nsamples)
{
    RUNTIME_ERROR("need to be overloaded");
}

//------------------------------------------------------------------------------

double CEnergyDerProxy::GetValue(int ibin,int cv,EProxyRealm realm) const
{
    RUNTIME_ERROR("need to be overloaded");
    return(0.0);
}

//------------------------------------------------------------------------------

void CEnergyDerProxy::SetNCorr(double ncorr)
{
    RUNTIME_ERROR("need to be overloaded");
}


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
