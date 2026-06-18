// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <CSTProxy_dAdx.hpp>
#include <CSTProxy_Ecorr.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dAdx::CCSTProxy_dAdx(void)
{
    RegisterRealm(CST_dAdx,         "dA/dx",        "CST", "dA(x)=|<lam> dx| + dA{CST}corr (MD subsystem)");
    RegisterRealm(CST_dAdx_TdS,     "dA/dx(TDS)",   "CST", "dA(x)=|<lam> dx| + dA{CST}corr (TDS subsystem)");
    RegisterRealm(CST_dAdx_VF,      "dA/dx(VF)",    "CST", "dA(x)=|<ICFFW> dx| (VF subsystem)");
}

//------------------------------------------------------------------------------

CCSTProxy_dAdx::~CCSTProxy_dAdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dAdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    switch(RealmID){
    // -------------------
        case(CST_dAdx):
            return(Accu->GetData("NSAMPLES",ibin));
    // -------------------
        default:
            return(Accu->GetData("NTDS",ibin));
    }
}

//------------------------------------------------------------------------------

void CCSTProxy_dAdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    switch(RealmID){
    // -------------------
        case(CST_dAdx):
            Accu->SetData("NSAMPLES",ibin,nsamples);
        break;
    // -------------------  
        default:
            Accu->SetData("NTDS",ibin,nsamples);
        break;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyProxyPtr CCSTProxy_dAdx::GetEnergyCorrection(void)
{
    CEnergyProxyPtr ene_proxy;
    if( RealmID == CST_dAdx ){
        ene_proxy = CCSTProxy_Ecorr_Ptr(new CCSTProxy_Ecorr);
        ene_proxy->SetRealm(CST_dA_corr);
        ene_proxy->Init(Accu);
    }
    if( RealmID == CST_dAdx_TdS ){
        ene_proxy = CCSTProxy_Ecorr_Ptr(new CCSTProxy_Ecorr);
        ene_proxy->SetRealm(CST_dA_corr_TdS);
        ene_proxy->Init(Accu);
    }
    return(ene_proxy);
}

//------------------------------------------------------------------------------

double CCSTProxy_dAdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  value   = 0.0;  // result
    double  sd      = 0.0;  // unbiased sample standard deviation
    double  sem     = 0.0;  // standard error of the sample result

// get requested data
    switch(RealmID){
    // -------------------
        case(CST_dAdx):     // this adds MTC correction
            GetMeanValue("LAMBDA",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dAdx_TdS): // this adds MTC correction
            GetMeanValue("LAMTDS",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dAdx_VF):
            GetWMeanValue("ICFFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
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



