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

#include <ABFProxy_PP.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFProxy_PP::CABFProxy_PP(void)
{
    SetType(ABF_PPNP_HH);
    Requires.push_back("ABF");
}

//------------------------------------------------------------------------------

CABFProxy_PP::~CABFProxy_PP(void)
{
}

//------------------------------------------------------------------------------

bool CABFProxy_PP::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "ABF" ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

void CABFProxy_PP::SetType(EABFPPNPType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(ABF_PPNP_HH):
            Provide = "ABF PP-HH";
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFProxy_PP::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

double CABFProxy_PP::GetValue(int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    int icv = 0;

    double  nsamples = Accu->GetData("NSAMPLES",ibin);

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    double  micf     = 0.0;
    double  m2icf    = 0.0;

    switch(Type){
    // -------------------
        case(ABF_PPNP_HH):
            micf     = Accu->GetData("MPP",ibin,icv);
            m2icf    = Accu->GetData("M2PP",ibin,icv);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    switch(realm){
    // -------------------
        // -------------------
        case(E_PROXY_VALUE):
            return( micf );
        // -------------------
        case(E_PROXY_SIGMA):
        case(E_PROXY_ERROR):
            return( sqrt(m2icf / nsamples) );
        // -------------------

           // return( sqrt(m2icf ) / nsamples );
    // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



