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

#include <CSTProxy_mTdS.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_mTdS::CCSTProxy_mTdS(void)
{
    SetType(CST_TdS_LTFW);
    Requires.push_back("CST");
}

//------------------------------------------------------------------------------

CCSTProxy_mTdS::~CCSTProxy_mTdS(void)
{
}

//------------------------------------------------------------------------------

bool CCSTProxy_mTdS::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "CST" ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

void CCSTProxy_mTdS::SetType(ECSTTdSType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(CST_TdS_LT):
            Provide = "CST -TdS(x)^{c}";    // entropy of the constrained system
        break;
    // -------------------
        case(CST_TdS_LTFW):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Epot)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        case(CST_TdS_LIFW):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Erst)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        case(CST_TdS_LPFW):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Ekin)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        case(CST_TdS_LRFW):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Ekin)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        case(CST_TdS_LKFW):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Ekin)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_mTdS::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_mTdS::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_mTdS::GetValue(int ibin,int icv,EProxyRealm realm) const
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
    switch(Type){
    // -------------------
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



