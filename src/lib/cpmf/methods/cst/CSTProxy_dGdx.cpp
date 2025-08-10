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

#include <CSTProxy_dGdx.hpp>
#include <CSTProxy_MTC.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dGdx::CCSTProxy_dGdx(void)
{
    Requires.push_back("CST");

    SupportedRealms["dG/dx"]        = CST_dGdx;
    SupportedRealms["LAMBDA/dx"]    = CST_LAMBDAdx;
    SupportedRealms["MICF/dx"]      = CST_MICFdx;
    SupportedRealms["MICFFW/dx"]    = CST_MICFFWdx;
    SupportedRealms["MICFPFW/dx"]   = CST_MICFPFWdx;
    SupportedRealms["MICFKFW/dx"]   = CST_MICFKFWdx;
}

//------------------------------------------------------------------------------

CCSTProxy_dGdx::~CCSTProxy_dGdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CCSTProxy_dGdx::SetType(const CSmallString& realm)
{
    if( SupportedRealms.count(realm) == 0 ) return(false);
    Realm = realm;
    SetType(SupportedRealms[realm]);
    return(true);
}

//------------------------------------------------------------------------------

void CCSTProxy_dGdx::SetType(ECSTdGdxType type)
{
    Type = type;
    Description = GetTypeDescription(Type);
}

//------------------------------------------------------------------------------

const CSmallString CCSTProxy_dGdx::GetTypeDescription(ECSTdGdxType type)
{
    switch(type){
    // -------------------
        case(CST_dGdx):
            return("CST dG(x) (|<l> dx| + MTC)");
    // -------------------
        case(CST_LAMBDAdx):
            return("CST |<lambda> dx|");
    // -------------------
        case(CST_MICFdx):
            return("CST |MICF(x) dx|");
    // -------------------
        case(CST_MICFFWdx):
            return("CST |MICF(x)FW dx|");
    // -------------------
        case(CST_MICFPFWdx):
            return("CST |MICFP(x)FW dx|");
    // -------------------
        case(CST_MICFKFWdx):
            return("CST |MICFK(x)FW dx|");
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

CEnergyProxyPtr CCSTProxy_dGdx::GetEnergyCorrection(void)
{
    CEnergyProxyPtr ene_proxy;
    if( Type == CST_dGdx ){
        ene_proxy = CCSTProxy_MTC_Ptr(new CCSTProxy_MTC);
        ene_proxy->Init(Accu);
    }
    return(ene_proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dGdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    switch(Type){
    // -------------------
        case(CST_dGdx):
        case(CST_LAMBDAdx):
            return(Accu->GetData("NSAMPLES",ibin));
    // -------------------
        case(CST_MICFdx):
        case(CST_MICFFWdx):
        case(CST_MICFPFWdx):
        case(CST_MICFKFWdx):
            return(Accu->GetData("NTDS",ibin));
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

void CCSTProxy_dGdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    switch(Type){
    // -------------------
        case(CST_dGdx):
        case(CST_LAMBDAdx):
            Accu->SetData("NSAMPLES",ibin,nsamples);
    // -------------------
        case(CST_MICFdx):
        case(CST_MICFFWdx):
        case(CST_MICFPFWdx):
        case(CST_MICFKFWdx):
            Accu->SetData("NTDS",ibin,nsamples);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

double CCSTProxy_dGdx::GetValue(int ibin,int icv,EProxyRealm realm) const
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
    switch(Type){
    // -------------------
        case(CST_dGdx):       // this requires MTC correction
        case(CST_LAMBDAdx): {
            mean        = Accu->GetData("MLAMBDA",ibin,icv);
            double M2   = Accu->GetData("M2LAMBDA",ibin,icv);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // ------------------
        case(CST_MICFdx): {
            mean        = Accu->GetData("MICF",ibin,icv);
            double M2   = Accu->GetData("M2ICF",ibin,icv);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_MICFFWdx): {
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MICFFW",ibin,icv);
            double M2       = Accu->GetData("M2ICFFW",ibin,icv);

            // https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf

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
        case(CST_MICFPFWdx):  {
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MICFPFW",ibin,icv);
            double M2       = Accu->GetData("M2ICFPFW",ibin,icv);

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
        case(CST_MICFKFWdx):  {
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MICFKFW",ibin,icv);
            double M2       = Accu->GetData("M2ICFKFW",ibin,icv);

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



