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

#include <EnergyProxy.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyProxy::CEnergyProxy(void)
{
    Method = "NONE";
}

//------------------------------------------------------------------------------

CEnergyProxy::~CEnergyProxy(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CEnergyProxy::GetNumOfSamples(int ibin) const
{
    RUNTIME_ERROR("need to be overloaded");
    return(0);
}


//------------------------------------------------------------------------------

void CEnergyProxy::SetNumOfSamples(int ibin,int nsamples)
{
    RUNTIME_ERROR("need to be overloaded");
}

//------------------------------------------------------------------------------

double CEnergyProxy::GetValue(int ibin,EProxyRealm realm) const
{
    RUNTIME_ERROR("need to be overloaded");
    return(0.0);
}

//------------------------------------------------------------------------------

void CEnergyProxy::SetNCorr(double ncorr)
{
    RUNTIME_ERROR("need to be overloaded");
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
