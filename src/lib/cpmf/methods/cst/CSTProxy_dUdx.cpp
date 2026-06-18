// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <CSTProxy_dUdx.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dUdx::CCSTProxy_dUdx(void)
{
    RegisterRealm(CST_dUdx_VF,  "dU/dx",    "CST",  "dU=|[<ICFP>_FW - cov(ICF,Eint)_FW/(k_B*T)] dx|");
}

//------------------------------------------------------------------------------

CCSTProxy_dUdx::~CCSTProxy_dUdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dUdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dUdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_dUdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double temp  = Accu->GetTemperature();

    double  value   = 0.0;  // result
    double  sd      = 0.0;  // unbiased sample standard deviation
    double  sem     = 0.0;  // standard error of the sample result

    switch(RealmID){
    // -------------------
        case(CST_dUdx_VF): { 
            double icfpfw_mean = 0.0;
            double icfpfw_sd = 0.0;
            double icfpfw_sem = 0.0;

            GetWMeanValue("ICFPFW","FW",icfpfw_mean,icfpfw_sd,icfpfw_sem,realm==E_PROXY_MEAN,ibin,icv);

            double c11iifw_cval = 0.0;
            double c11iifw_sd = 0.0;
            double c11iifw_sem = 0.0;

            GetWCovarianceValue("C11IIFW","FW",c11iifw_cval,c11iifw_sd,c11iifw_sem,realm==E_PROXY_MEAN,ibin,icv);

            value  =  icfpfw_mean - c11iifw_cval / (temp * PMF_Rgas);

            // https://en.wikipedia.org/wiki/Propagation_of_uncertainty
            // approximative estimates - icfkfw_mean and c11iifw_cval are considered as independent
            sd = sqrt( icfpfw_sd*icfpfw_sd + (c11iifw_sd / (temp * PMF_Rgas)) * (c11iifw_sd / (temp * PMF_Rgas)) );
            sem = sqrt( icfpfw_sem*icfpfw_sem + (c11iifw_sem / (temp * PMF_Rgas)) * (c11iifw_sem / (temp * PMF_Rgas)) );
        }        
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
