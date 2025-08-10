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

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dG::CCSTProxy_dG(void)
{
    Requires.push_back("CST");
    SetType(CST_dG);
}

//------------------------------------------------------------------------------

CCSTProxy_dG::~CCSTProxy_dG(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCSTProxy_dG::SetType(ECSTdGType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(CST_dG):
            Provide = "CST dG(x) (|<l> dx| + MTC)";
    // -------------------
        case(CST_MICF):
            Provide = "CST |MICF(x) dx|";
        break;
    // -------------------
        case(CST_MICFFW):
            Provide = "CST |MICF(x)FW dx|";
        break;
    // -------------------
        case(CST_MICFPFW):
            Provide = "CST |MICFP(x)FW dx|";
        break;
    // -------------------
        case(CST_MICFKFW):
            Provide = "CST |MICFK(x)FW dx|";
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//------------------------------------------------------------------------------

bool CCSTProxy_dG::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "CST" ) return(true);
    return(false);
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
        case(CST_dG): {
            // this requires MTC correction
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



