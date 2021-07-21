// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <ABPProxy_dG.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABPProxy_dG::CABPProxy_dG(void)
{
    Requires.push_back("ABP");
    Provide = "ABP dG(x)";
}

//------------------------------------------------------------------------------

CABPProxy_dG::~CABPProxy_dG(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPProxy_dG::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    return( Accu->GetData("NSAMPLES",ibin) );
}

//------------------------------------------------------------------------------

double CABPProxy_dG::GetValue(int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double pop  = Accu->GetData("POP",ibin);
    double temp = Accu->GetTemperature();
    double ene  = 0.0;

    switch(realm){
        // -------------------
        case(E_PROXY_VALUE):
            ene = - temp*PMF_Rgas*log(pop);
            return( ene );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return( ene );
}


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



