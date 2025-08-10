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

#include <CSTProxy_dG.hpp>
#include <CSTProxy_MTC.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dG::CCSTProxy_dG(void)
{
    Requires.push_back("CST");

    SupportedRealms["dG/dx"]        = CST_dG;
    SupportedRealms["LAMBDA/dx"]    = CST_LAMBDA;
    SupportedRealms["MICF/dx"]      = CST_MICF;
    SupportedRealms["MICFFW/dx"]    = CST_MICFFW;
    SupportedRealms["MICFPFW/dx"]   = CST_MICFPFW;
    SupportedRealms["MICFKFW/dx"]   = CST_MICFKFW;
}

//------------------------------------------------------------------------------

CCSTProxy_dG::~CCSTProxy_dG(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CCSTProxy_dG::SetType(const CSmallString& realm)
{
    if( SupportedRealms.count(realm) == 0 ) return(false);
    Realm = realm;
    SetType(SupportedRealms[realm]);
    return(true);
}

//------------------------------------------------------------------------------

void CCSTProxy_dG::SetType(ECSTdGType type)
{
    Type = type;
    Description = GetTypeDescription(Type);
}

//------------------------------------------------------------------------------

const CSmallString CCSTProxy_dG::GetTypeDescription(ECSTdGType type)
{
    switch(type){
    // -------------------
        case(CST_dG):
            return("CST dG(x) (|<l> dx| + MTC)");
    // -------------------
        case(CST_LAMBDA):
            return("CST |<lambda> dx|");
    // -------------------
        case(CST_MICF):
            return("CST |MICF(x) dx|");
    // -------------------
        case(CST_MICFFW):
            return("CST |MICF(x)FW dx|");
    // -------------------
        case(CST_MICFPFW):
            return("CST |MICFP(x)FW dx|");
    // -------------------
        case(CST_MICFKFW):
            return("CST |MICFK(x)FW dx|");
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

CEnergyProxyPtr CCSTProxy_dG::GetEnergyCorrection(void)
{
    CEnergyProxyPtr ene_proxy;
    if( Type == CST_dG ){
        ene_proxy = CCSTProxy_MTC_Ptr(new CCSTProxy_MTC);
        ene_proxy->Init(Accu);
    }
    return(ene_proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dG::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    switch(Type){
    // -------------------
        case(CST_dG):
        case(CST_LAMBDA):
            return(Accu->GetData("NSAMPLES",ibin));
    // -------------------
        case(CST_MICF):
        case(CST_MICFFW):
        case(CST_MICFPFW):
        case(CST_MICFKFW):
            return(Accu->GetData("NTDS",ibin));
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

void CCSTProxy_dG::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    switch(Type){
    // -------------------
        case(CST_dG):
        case(CST_LAMBDA):
            Accu->SetData("NSAMPLES",ibin,nsamples);
    // -------------------
        case(CST_MICF):
        case(CST_MICFFW):
        case(CST_MICFPFW):
        case(CST_MICFKFW):
            Accu->SetData("NTDS",ibin,nsamples);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

double CCSTProxy_dG::GetValue(int ibin,int icv,EProxyRealm realm) const
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
        case(CST_dG):       // this requires MTC correction
        case(CST_LAMBDA): {
            mean        = Accu->GetData("MLAMBDA",ibin,icv);
            double M2   = Accu->GetData("M2LAMBDA",ibin,icv);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // ------------------
        case(CST_MICF): {
            mean        = Accu->GetData("MICF",ibin,icv);
            double M2   = Accu->GetData("M2ICF",ibin,icv);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_MICFFW): {
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
        case(CST_MICFPFW):  {
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
        case(CST_MICFKFW):  {
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



