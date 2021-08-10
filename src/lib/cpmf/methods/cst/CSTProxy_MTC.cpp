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

#include <CSTProxy_MTC.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_MTC::CCSTProxy_MTC(void)
{
    Requires.push_back("CST");
    Provide = "CST dG(x)^{MTC}";      // metric tensor correction
}

//------------------------------------------------------------------------------

CCSTProxy_MTC::~CCSTProxy_MTC(void)
{
}

//------------------------------------------------------------------------------

bool CCSTProxy_MTC::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "CST" ) return(true);
    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_MTC::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_MTC::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NSAMPLES",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_MTC::GetValue(int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = Accu->GetData("NSAMPLES",ibin);
    double  misrz    = Accu->GetData("MISRZ",ibin);
    double  temp     = Accu->GetTemperature();

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
// MTC correction
        // -------------------
        case(E_PROXY_VALUE): {
            double value = - PMF_Rgas * temp * log(misrz) ;
            return( value );
        }
        // -------------------
        case(E_PROXY_SIGMA):
            return( 0.0 );      // FIXME
        // -------------------
        case(E_PROXY_ERROR):
            return( 0.0 );      // FIXME
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



