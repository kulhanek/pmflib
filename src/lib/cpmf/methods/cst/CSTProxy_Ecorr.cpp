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

#include <CSTProxy_Ecorr.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_Ecorr::CCSTProxy_Ecorr(void)
{
    RegisterRealm(CST_dA_corr,      "dA_corr",      "CST", "dA{CST}corr");
    RegisterRealm(CST_dA_corr_TdS,  "dA_corr(TdS)", "CST", "dA{CST}corr (TDS subsystem)");
    RegisterRealm(CST_mTdS_corr,    "mTdS_corr",    "CST", "-TdS{CST}corr");
    RegisterRealm(CST_mTdS_corr,    "-TdS_corr",    "CST", "-TdS{CST}corr");
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
        break;
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
        break;
    }
}

//------------------------------------------------------------------------------

double CCSTProxy_Ecorr::GetValue(int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  temp     = Accu->GetTemperature();

    double  value   = 0.0;  // result
    double  sd      = 0.0;  // unbiased sample standard deviation
    double  sem     = 0.0;  // standard error of the sample result

// get requested data
    switch(RealmID){
    // -------------------
        case(CST_dA_corr): {
            double fw_mean = 0.0;
            double fw_sd = 0.0;
            double fw_sem = 0.0;

            GetMeanValue("FW",fw_mean,fw_sd,fw_sem,realm==E_PROXY_MEAN,ibin);

            value = - PMF_Rgas * temp * log(fw_mean);

            // https://en.wikipedia.org/wiki/Propagation_of_uncertainty
            if( fw_mean != 0.0 ){
                sd = fabs(- PMF_Rgas * temp * fw_sd / fw_mean);
                sem = fabs( - PMF_Rgas * temp * fw_sem / fw_mean );
            }
        }
        break;
    // -------------------
        case(CST_dA_corr_TdS): {
            double fw_mean = 0.0;
            double fw_sd = 0.0;
            double fw_sem = 0.0;

            GetMeanValue("FWTDS",fw_mean,fw_sd,fw_sem,realm==E_PROXY_MEAN,ibin);

            value = - PMF_Rgas * temp * log(fw_mean);

            // https://en.wikipedia.org/wiki/Propagation_of_uncertainty
            if( fw_mean != 0.0 ){
                sd = fabs(- PMF_Rgas * temp * fw_sd / fw_mean);
                sem = fabs( - PMF_Rgas * temp * fw_sem / fw_mean );
            }
        }
        break;
    // -------------------
        case(CST_mTdS_corr): {
            double fw_mean = 0.0;
            double fw_sd = 0.0;
            double fw_sem = 0.0;

            GetMeanValue("FWTDS",fw_mean,fw_sd,fw_sem,realm==E_PROXY_MEAN,ibin);

            double corr1_mean = - PMF_Rgas * temp * log(fw_mean);
            double corr1_sd = 0.0;
            double corr1_sem = 0.0;

            // https://en.wikipedia.org/wiki/Propagation_of_uncertainty
            if( fw_mean != 0.0 ){
                corr1_sd = fabs(- PMF_Rgas * temp * fw_sd / fw_mean);
                corr1_sem = fabs( - PMF_Rgas * temp * fw_sem / fw_mean );
            }

            double c11_cval = 0.0;
            double c11_sd = 0.0;
            double c11_sem = 0.0;

            GetCovarianceValue("C11ZH",c11_cval,c11_sd,c11_sem,realm==E_PROXY_MEAN,ibin);

            value  = corr1_mean - c11_cval / fw_mean;

            // https://en.wikipedia.org/wiki/Propagation_of_uncertainty
            // approximative estimates - terms are considered as independent
            if( fw_mean != 0.0 ){
                const double f  = fw_mean;
                const double c  = c11_cval;
                const double f2 = f*f;
                const double f4 = f2*f2;

                sd = sqrt( corr1_sd*corr1_sd
                        + c11_sd*c11_sd / f2
                        + c*c * fw_sd*fw_sd / f4 );

                sem = sqrt( corr1_sem*corr1_sem
                        + c11_sem*c11_sem / f2
                        + c*c * fw_sem*fw_sem / f4 );
            }
        }
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
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



