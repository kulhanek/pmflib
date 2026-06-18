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

#include <CSTProxy_dYdx.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dYdx::CCSTProxy_dYdx(void)
{
    RegisterRealm(CST_dLAMBDAdx,    "dLAMBDA/dx(MD)",   "CST", "|<lambda> dx| (MD subsystem)");
    RegisterRealm(CST_dLAMTDSdx,    "dLAMBDA/dx(TDS)",  "CST", "|<lambda> dx| (TDS subsystem)");

    RegisterRealm(CST_dMICFdx,      "dMICF/dx",         "CST", "|<ICF> dx| (VF subsystem)");
    RegisterRealm(CST_dMICFFWdx,    "dMICFFW/dx",       "CST", "|<ICF>_FW dx| (VF subsystem)");
    RegisterRealm(CST_dMICFPFWdx,   "dMICFPFW/dx",      "CST", "|<ICFP>_FW dx| (VF subsystem)");
    RegisterRealm(CST_dMICFKFWdx,   "dMICFKFW/dx",      "CST", "|<ICFK>_FW dx| (VF subsystem)");
}

//------------------------------------------------------------------------------

CCSTProxy_dYdx::~CCSTProxy_dYdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dYdx::GetNumOfSamples(int ibin) const
{
    switch(RealmID){
    // -------------------
        case(CST_dLAMBDAdx):
            return(Accu->GetData("NSAMPLES",ibin));
    // -------------------
        default:
            return(Accu->GetData("NTDS",ibin));
    }
}

//------------------------------------------------------------------------------

void CCSTProxy_dYdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    switch(RealmID){
    // -------------------
        case(CST_dLAMBDAdx):
            Accu->SetData("NSAMPLES",ibin,nsamples);
        break;
    // -------------------
        default:
            Accu->SetData("NTDS",ibin,nsamples);
        break;
    }
}

//------------------------------------------------------------------------------

double CCSTProxy_dYdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  value   = 0.0;  // result
    double  sd      = 0.0;  // unbiased sample standard deviation
    double  sem     = 0.0;  // standard error of the sample result

    switch(RealmID){
    // -------------------
        case(CST_dLAMBDAdx):
            GetMeanValue("LAMBDA",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dLAMTDSdx):
            GetMeanValue("LAMTDS",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dMICFdx):
            GetMeanValue("ICF",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dMICFFWdx):
            GetWMeanValue("ICFFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dMICFPFWdx):
            GetWMeanValue("ICFPFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
        break;
    // -------------------
        case(CST_dMICFKFWdx):
            GetWMeanValue("ICFKFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin,icv);
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
