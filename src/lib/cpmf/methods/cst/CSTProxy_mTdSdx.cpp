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

#include <CSTProxy_mTdSdx.hpp>
#include <CSTProxy_MTC.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_mTdSdx::CCSTProxy_mTdSdx(void)
{
//    Requires.push_back("CST");
//
//    SupportedRealms["-TdS/dx"]      = CST_TdS;
//    SupportedRealms["mTdS/dx"]      = CST_TdS;
//
//    SupportedRealms["-TdS_LT/dx"]   = CST_TdS_LT;
//    SupportedRealms["mTdS_LT/dx"]   = CST_TdS_LT;
//
//    SupportedRealms["-TdS_LTFW/dx"] = CST_TdS_LTFW;
//    SupportedRealms["mTdS_LTFW/dx"] = CST_TdS_LTFW;
//
//    SupportedRealms["-TdS_LIFW/dx"] = CST_TdS_LIFW;
//    SupportedRealms["mTdS_LIFW/dx"] = CST_TdS_LIFW;
//
//    SupportedRealms["-TdS_LPFW/dx"] = CST_TdS_LPFW;
//    SupportedRealms["mTdS_LPFW/dx"] = CST_TdS_LPFW;
//
//    SupportedRealms["-TdS_LRFW/dx"] = CST_TdS_LRFW;
//    SupportedRealms["mTdS_LRFW/dx"] = CST_TdS_LRFW;
//
//    SupportedRealms["-TdS_LKFW/dx"] = CST_TdS_LKFW;
//    SupportedRealms["mTdS_LKFW/dx"] = CST_TdS_LKFW;
//
//    SupportedRealms["-TdS_II/dx"]   = CST_TdS_II;
//    SupportedRealms["mTdS_II/dx"]   = CST_TdS_II;
//
//    SupportedRealms["-TdS_IIFW/dx"] = CST_TdS_IIFW;
//    SupportedRealms["mTdS_IIFW/dx"] = CST_TdS_IIFW;

//    switch(type){
//    // -------------------
//        case(CST_TdS):
//            return("CST -TdS(x)");
//    // -------------------
//        case(CST_TdS_LT):
//            return("CST -TdS_LT(x)   cov(lambda,Epot)");
//    // -------------------
//        case(CST_TdS_LTFW):
//            return("CST -TdS_LTFW(x) cov(lambda,Epot) - Fixman weighted");
//    // -------------------
//        case(CST_TdS_II):
//            return("CST -TdS_II(x)   cov(ICF,Eint)");
//    // -------------------
//        case(CST_TdS_IIFW):
//            return("CST -TdS_IIFW(x) cov(ICF,Eint) - Fixman weighted");
//    // -------------------
//        default:
//            RUNTIME_ERROR("unsupported type");
}

//------------------------------------------------------------------------------

CCSTProxy_mTdSdx::~CCSTProxy_mTdSdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyProxyPtr CCSTProxy_mTdSdx::GetEnergyCorrection(void)
{
    CEnergyProxyPtr ene_proxy;
    if( RealmID == CST_TdS ){
        ene_proxy = CCSTProxy_MTC_Ptr(new CCSTProxy_MTC);
        ene_proxy->Init(Accu);
    }
    return(ene_proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_mTdSdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_mTdSdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_mTdSdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  ncorr    = Accu->GetNCorr();
    double  temp     = Accu->GetTemperature();
    double  mean     = 0.0; // sample mean
    double  samvar   = 0.0; // sample variance
    double  meanvar  = 0.0; // variance of sample mean

// do we have enough samples?
    double nsamples    = GetNumOfSamples(ibin);
    if( nsamples <= 0 ) return(mean);

// get requested data
    switch(RealmID){
    // -------------------
        case(CST_TdS):      // plus MTC
        case(CST_TdS_LT):{
            double C        = Accu->GetData("C11LT",ibin,icv);
            mean            = C / nsamples;
            samvar          = 0.0;  // FIXME
            meanvar         = 0.0;
        }
        break;
    // -------------------
        case(CST_TdS_LTFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double C        = Accu->GetData("C11LTFW",ibin,icv);
            mean            = C / fwsum;
            samvar          = 0.0;  // FIXME
            meanvar         = 0.0;
        }
        break;
    // -------------------
        case(CST_TdS_II):{
            double C        = Accu->GetData("C11II",ibin,icv);
            mean            = C / nsamples;
            samvar          = 0.0;  // FIXME
            meanvar         = 0.0;
        }
        break;
    // -------------------
        case(CST_TdS_IIFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double C        = Accu->GetData("C11IIFW",ibin,icv);
            mean            = C / fwsum;
            samvar          = 0.0;  // FIXME
            meanvar         = 0.0;
        }
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

// return result
    switch(realm){
        // -------------------
        case(E_PROXY_VALUE):
            return( mean / (temp * PMF_Rgas) );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(samvar) / (temp * PMF_Rgas) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(ncorr * meanvar) / (temp * PMF_Rgas) );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



