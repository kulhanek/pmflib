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

#include <ABFProxy_dH.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFProxy_dH::CABFProxy_dH(void)
{
    SetType(ABF_MICFP);

    Requires.push_back("ABF");
}

//------------------------------------------------------------------------------

CABFProxy_dH::~CABFProxy_dH(void)
{
}

//------------------------------------------------------------------------------

bool CABFProxy_dH::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "ABF" ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

void CABFProxy_dH::SetType(EABFdHType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(ABF_MICFP):
            Provide = "ABF ICFP";
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFProxy_dH::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CABFProxy_dH::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CABFProxy_dH::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = 0.0;
    double  micf     = 0.0;
    double  m2icf    = 0.0;
    double  ncorr    = Accu->GetNCorr();

    switch(Type){
    // -------------------
        case(ABF_MICFP):
            nsamples = Accu->GetData("NTDS",ibin);
            micf     = Accu->GetData("MHICFP",ibin,icv);
            m2icf    = Accu->GetData("M2HICFP",ibin,icv);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
// mean force
        // -------------------
        case(E_PROXY_VALUE):
            return( micf );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(m2icf / nsamples) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(m2icf * ncorr) / nsamples );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



