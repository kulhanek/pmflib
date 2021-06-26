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

#include <CSTProxy_dH.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dH::CCSTProxy_dH(void)
{
    Require = "CST";
    Provide = "CST dH(x)";
    NCorr = 1.0;
}

//------------------------------------------------------------------------------

CCSTProxy_dH::~CCSTProxy_dH(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dH::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dH::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
}

//------------------------------------------------------------------------------

double CCSTProxy_dH::GetValue( int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = Accu->GetData("NSAMPLES",ibin);
    double  mepot    = Accu->GetData("MEPOT",ibin);
    double  m2epot   = Accu->GetData("M2EPOT",ibin);

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
// mean force
        // -------------------
        case(E_PROXY_VALUE):
            return( mepot );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(m2epot / nsamples) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(m2epot * NCorr) / nsamples );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//------------------------------------------------------------------------------

void CCSTProxy_dH::SetNCorr(double ncorr)
{
    NCorr  = ncorr;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



