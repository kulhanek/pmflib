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

    double mtdpot = Accu->GetData("MTDPOT",ibin);

    switch(realm){
        // -------------------
        case(E_PROXY_VALUE):
            return( - mtdpot );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return( mtdpot );
}

//    const double Rgas = 0.0019872065;
//
//    double m = 1.0;
//    for(int i=0; i < Accu->GetNumOfBins(); i++){
//        double pop = Accu->GetPop(i);
//        if( pop > m ) m = pop;
//    }
//    for(int i=0; i < Accu->GetNumOfBins(); i++){
//        double pop = Accu->GetPop(i);
//        double ene = - Accu->GetTemperature()*Rgas*log(pop/m);
//        FES->SetEnergy(i,ene);
//        FES->SetNumOfSamples(i,Accu->GetNumberOfABPSamples(i));
//    }

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



