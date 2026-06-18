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

#include <CSTProxy_dCdx.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dCdx::CCSTProxy_dCdx(void)
{
    RegisterRealm(CST_TdS_LT,       "C11LT",            "CST", "|cov(lambda,Etot)/(k_B*T) dx| (TDS subsystem)");
    RegisterRealm(CST_TdS_LI,       "C11LI",            "CST", "|cov(lambda,Eint)/(k_B*T) dx| (TDS subsystem)");
    RegisterRealm(CST_TdS_LP,       "C11LP",            "CST", "|cov(lambda,Epot)/(k_B*T) dx| (TDS subsystem)");
    RegisterRealm(CST_TdS_LR,       "C11LR",            "CST", "|cov(lambda,Erst)/(k_B*T) dx| (TDS subsystem)");
    RegisterRealm(CST_TdS_LK,       "C11LK",            "CST", "|cov(lambda,Ekin)/(k_B*T) dx| (TDS subsystem)");

    RegisterRealm(CST_TdS_LTFW,     "C11LTFW",          "CST", "|cov(lambda,Etot)_FW/(k_B*T) dx| (TDS subsystem)");

    RegisterRealm(CST_TdS_II,       "C11II",            "CST", "|cov(ICF,Eint)/(k_B*T) dx| (VF subsystem)");
    RegisterRealm(CST_TdS_IIFW,     "C11IIFW",          "CST", "|cov(ICF,Eint)_FW/(k_B*T) dx| (VF subsystem)");
    RegisterRealm(CST_TdS_PIFW,     "C11PIFW",          "CST", "|cov(ICFP,Eint)_FW/(k_B*T) dx| (VF subsystem)");
    RegisterRealm(CST_TdS_KIFW,     "C11KIFW",          "CST", "|cov(ICFK,Eint)_FW/(k_B*T) dx| (VF subsystem)");
}

//------------------------------------------------------------------------------

CCSTProxy_dCdx::~CCSTProxy_dCdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dCdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dCdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_dCdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double temp  = Accu->GetTemperature();

    double  value   = 0.0;  // result
    double  sd      = 0.0;  // unbiased sample standard deviation
    double  sem     = 0.0;  // standard error of the sample result

    switch(RealmID){
        case(CST_TdS_LT):
            GetCovarianceValue("C11LT",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_LI):
            GetCovarianceValue("C11LI",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_LP):
            GetCovarianceValue("C11LP",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_LR):
            GetCovarianceValue("C11LR",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_LK):
            GetCovarianceValue("C11LK",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_LTFW):
            GetWCovarianceValue("C11LTFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_II):
            GetCovarianceValue("C11II",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_IIFW):
            GetWCovarianceValue("C11IIFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_PIFW):
            GetWCovarianceValue("C11PIFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_TdS_KIFW):
            GetWCovarianceValue("C11KIFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    value = value / (temp * PMF_Rgas);
    sd = sd / (temp * PMF_Rgas);
    sem = sem / (temp * PMF_Rgas);

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
