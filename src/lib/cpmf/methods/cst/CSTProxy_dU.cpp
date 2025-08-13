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

#include <CSTProxy_dU.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dU::CCSTProxy_dU(void)
{
//    Requires.push_back("CST");
//
//    SupportedRealms["dH"]       = CST_dH;
//    SupportedRealms["EINT"]     = CST_EINT;
//    SupportedRealms["EINTFW"]   = CST_EINTFW;

//    // -------------------
//        case(CST_dH):
//            return("dH(x)=<Eint>");
//    // -------------------
//        case(CST_EINT):
//            return("dH(x)=<Eint>");
//    // -------------------
//        case(CST_EINTFW):
//            return("dH(x)=<EintFW>");
//    // -------------------
}

//------------------------------------------------------------------------------

CCSTProxy_dU::~CCSTProxy_dU(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dU::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dU::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_dU::GetValue( int ibin,EProxyRealm realm) const
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

    switch(RealmID){
    // -------------------
        case(CST_EINT):{
            mean        = Accu->GetData("MEINT",ibin);
            double M2   = Accu->GetData("M2EINT",ibin);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_dH):
        case(CST_EINTFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MEINTFW",ibin);
            double M2       = Accu->GetData("M2EINTFW",ibin);

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



