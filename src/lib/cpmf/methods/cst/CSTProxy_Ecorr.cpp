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

#include <CSTProxy_Ecorr.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_Ecorr::CCSTProxy_Ecorr(void)
{
    RegisterRealm(CST_dA_corr,      "dA_corr",     "CST", "dA{CST}corr");
    RegisterRealm(CST_dA_corr_TdS,  "dA_corr_TdS", "CST", "dA{CST}corr - TdS source");
    RegisterRealm(CST_mTdS_corr,    "mTdS_corr",   "CST", "-TdS{CST}corr");
    RegisterRealm(CST_mTdS_corr,    "-TdS_corr",   "CST", "-TdS{CST}corr");
}

//------------------------------------------------------------------------------

CCSTProxy_Ecorr::~CCSTProxy_Ecorr(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_Ecorr::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    switch(RealmID){
    // -------------------
        case(CST_dA_corr):
            return(Accu->GetData("NSAMPLES",ibin));
    // -------------------
        case(CST_dA_corr_TdS):
        case(CST_mTdS_corr):
            return(Accu->GetData("NTDS",ibin));
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

void CCSTProxy_Ecorr::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    switch(RealmID){
    // -------------------
        case(CST_dA_corr):
            Accu->SetData("NSAMPLES",ibin,nsamples);
        break;
    // -------------------
        case(CST_dA_corr_TdS):
        case(CST_mTdS_corr):
            Accu->SetData("NTDS",ibin,nsamples);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

double CCSTProxy_Ecorr::GetValue(int ibin,EProxyRealm realm) const
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
        case(CST_dA_corr): {
            double mfw  = Accu->GetData("MFW",ibin);
            double m2fw = Accu->GetData("M2FW",ibin);

            mean        = - PMF_Rgas * temp * log(mfw);

            samvar      = m2fw / nsamples;
            samvar      = (PMF_Rgas*temp/mfw)*(PMF_Rgas*temp/mfw) * samvar;

            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_dA_corr_TdS): {
            double mfw  = Accu->GetData("MFWTDS",ibin);
            double m2fw = Accu->GetData("M2FWTDS",ibin);

            mean        = - PMF_Rgas * temp * log(mfw);

            samvar      = m2fw / nsamples;
            samvar      = (PMF_Rgas*temp/mfw)*(PMF_Rgas*temp/mfw) * samvar;

            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_mTdS_corr): {
            double mfw   = Accu->GetData("MFWTDS",ibin);
            double corr1 = PMF_Rgas * temp * log(mfw);

            double C     = Accu->GetData("C11ZH",ibin);
            double corr2 = C / nsamples / mfw;

            mean         = - (corr1 + corr2);
            samvar       = 0.0; // FIXME
            meanvar      = 0.0;
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
            return( mean );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(samvar) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(ncorr * meanvar) );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



