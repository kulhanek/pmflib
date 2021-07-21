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

#include <MTDProxy_dG.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDProxy_dG::CMTDProxy_dG(void)
{
    Requires.push_back("MTD");
    Provide = "MTD dG(x)";
}

//------------------------------------------------------------------------------

CMTDProxy_dG::~CMTDProxy_dG(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDProxy_dG::IsWTMeta(void)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    return( Accu->HasSectionData("MTD-WT") );
}

//------------------------------------------------------------------------------

int CMTDProxy_dG::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    return( Accu->GetData("NSAMPLES",ibin) );
}

//------------------------------------------------------------------------------

double CMTDProxy_dG::GetValue(int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double mtdpot = Accu->GetData("MTDPOT",ibin);
    double fact = 1.0;

    if( Accu->HasSectionData("MTD-WT") ){
        // well-tempered metadynamics
        double temp = Accu->GetTemperature();
        double wtem = Accu->GetData("MTD-WT",0);
        if( wtem > 0 ){
            fact = (temp + wtem) / wtem;
        } else {
            RUNTIME_ERROR("MTD-WT temerature is not greater than zero");
        }
    }

    switch(realm){
        // -------------------
        case(E_PROXY_VALUE):
            return( - mtdpot * fact );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return( mtdpot );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



