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
    RegisterRealm(CST_dU,       "dU",       "CST", "dU=<Etot>FW");
    RegisterRealm(CST_ETOT,     "<Etot>",   "CST", "<Etot>");
    RegisterRealm(CST_ETOTFW,   "<Etot>FW", "CST", "<Etot>FW");
    RegisterRealm(CST_EINT,     "<Eint>",   "CST", "<Eint>");
    RegisterRealm(CST_EINTFW,   "<Eint>FW", "CST", "<Eint>FW");
    RegisterRealm(CST_EPOT,     "<Epot>",   "CST", "<Epot>");
    RegisterRealm(CST_EPOTFW,   "<Epot>FW", "CST", "<Epot>FW");
    RegisterRealm(CST_ERST,     "<Erst>",   "CST", "<Erst>");
    RegisterRealm(CST_ERSTFW,   "<Erst>FW", "CST", "<Erst>FW");
    RegisterRealm(CST_EKIN,     "<Ekin>",   "CST", "<Ekin>");
    RegisterRealm(CST_EKINFW,   "<Ekin>FW", "CST", "<Ekin>FW");
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
        case(CST_ETOT):{
            mean        = Accu->GetData("METOT",ibin);
            double M2   = Accu->GetData("M2ETOT",ibin);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_dU):
        case(CST_ETOTFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("METOTFW",ibin);
            double M2       = Accu->GetData("M2ETOTFW",ibin);

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
        case(CST_EINT):{
            mean        = Accu->GetData("MEINT",ibin);
            double M2   = Accu->GetData("M2EINT",ibin);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
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
        case(CST_EPOT):{
            mean        = Accu->GetData("MEPOT",ibin);
            double M2   = Accu->GetData("M2EPOT",ibin);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_EPOTFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MEPOTFW",ibin);
            double M2       = Accu->GetData("M2EPOTFW",ibin);

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
        case(CST_ERST):{
            mean        = Accu->GetData("MERST",ibin);
            double M2   = Accu->GetData("M2ERST",ibin);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_ERSTFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MERSTFW",ibin);
            double M2       = Accu->GetData("M2ERSTFW",ibin);

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
        case(CST_EKIN):{
            mean        = Accu->GetData("MEKIN",ibin);
            double M2   = Accu->GetData("M2EKIN",ibin);
            samvar      = M2 / nsamples;
            meanvar     = samvar / nsamples;
        }
        break;
    // -------------------
        case(CST_EKINFW):{
            double fwsum    = Accu->GetData("FWSUM",ibin);
            double fwsum2   = Accu->GetData("FWSUM2",ibin);

            mean            = Accu->GetData("MEKINFW",ibin);
            double M2       = Accu->GetData("M2EKINFW",ibin);

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



