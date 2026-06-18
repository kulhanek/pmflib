// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <CSTProxy_dU.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dU::CCSTProxy_dU(void)
{
    RegisterRealm(CST_dU,       "dU",       "CST", "dU=<Etot>_FW");
}

//------------------------------------------------------------------------------

CCSTProxy_dU::~CCSTProxy_dU(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dU::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dU::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_dU::GetValue( int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  value   = 0.0; // sample mean
    double  sd      = 0.0; // unbiased sample standard deviation
    double  sem     = 0.0; // standard error of the sample mean

    switch(RealmID){
        case(CST_dU):
            GetWMeanValue("ETOTFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
        break;
    }

// return result
    switch(realm){
        // -------------------
        case(E_PROXY_MEAN):
            return( value );
        // -------------------
        case(E_PROXY_SD):
            return( sd );
        // -------------------
        case(E_PROXY_SEM):
            return( sem );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



