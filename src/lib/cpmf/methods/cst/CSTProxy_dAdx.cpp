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
    RegisterRealm(CST_dAdx,     "dA/dx",    "CST", "dA(x)=|<lam> dx| + dA{CST}corr");
    RegisterRealm(CST_ICF,      "ICF",      "CST", "|ICF dx|");
    RegisterRealm(CST_ICFFW,    "ICFFW",    "CST", "|ICFFW dx|");
    RegisterRealm(CST_ICFPFW,   "ICFPFW",   "CST", "|ICFPFW dx|");
    RegisterRealm(CST_ICFKFW,   "ICFKFW",   "CST", "|ICFKFW dx|");
}

//------------------------------------------------------------------------------

CCSTProxy_dAdx::~CCSTProxy_dAdx(void)
{
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
    return(ene_proxy);
}

//------------------------------------------------------------------------------

double CCSTProxy_dAdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  ncorr    = Accu->GetNCorr();
    double  mean     = 0.0; // sample mean
    double  samvar   = 0.0; // sample variance
    double  meanvar  = 0.0; // variance of sample mean

// do we have enough samples?
    double nsamples    = GetNumOfSamples(ibin);
    if( nsamples <= 0 ) return(mean);

// get requested data
    switch(RealmID){
    // -------------------
        case(CST_dAdx): {  // this requires MTC correction
            mean        = Accu->GetData("MLAMBDA",ibin,icv);
            double M2   = Accu->GetData("M2LAMBDA",ibin,icv);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_ICF): {
            mean        = Accu->GetData("MICF",ibin,icv);
            double M2   = Accu->GetData("M2MICF",ibin,icv);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_ICFFW): {
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MICFFW",ibin);
            double M2       = Accu->GetData("M2ICFFW",ibin);

            // number of effective measurements
            double neff = fwsum2 / (fwsum * fwsum);

            // unbiased weighted sample variance
            samvar          = M2 / fwsum * neff / (neff - 1.0);

            // variance of the weighted mean
            // unbiased importance weights
            meanvar         = samvar / neff;
        }
        break;
    // -------------------
        case(CST_ICFPFW): {
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MICFPFW",ibin);
            double M2       = Accu->GetData("M2ICFPFW",ibin);

            // number of effective measurements
            double neff = fwsum2 / (fwsum * fwsum);

            // unbiased weighted sample variance
            samvar          = M2 / fwsum * neff / (neff - 1.0);

            // variance of the weighted mean
            // unbiased importance weights
            meanvar         = samvar / neff;
        }
        break;
    // -------------------
        case(CST_ICFKFW): {
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MICFKFW",ibin);
            double M2       = Accu->GetData("M2ICFKFW",ibin);

            // number of effective measurements
            double neff = fwsum2 / (fwsum * fwsum);

            // unbiased weighted sample variance
            samvar          = M2 / fwsum * neff / (neff - 1.0);

            // variance of the weighted mean
            // unbiased importance weights
            meanvar         = samvar / neff;
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



