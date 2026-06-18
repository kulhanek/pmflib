// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <ABFProxy_dUdx.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFProxy_dUdx::CABFProxy_dUdx(void)
{
//    Requires.push_back("ABF");

//    SupportedRealms["dG/dx"]        = CST_dG;
//    SupportedRealms["MICF/dx"]      = CST_MICF;
//    SupportedRealms["MICFFW/dx"]    = CST_MICFFW;
//    SupportedRealms["MICFPFW/dx"]   = CST_MICFPFW;
//    SupportedRealms["MICFKFW/dx"]   = CST_MICFKFW;
//
//    // -------------------
//        case(ABF_dH):
//            return("ABF dH(x) (based on derivatives)");
//    // -------------------
//        case(ABF_MICFP):
//            return("ABF ICFP(x)");
}

//------------------------------------------------------------------------------

CABFProxy_dUdx::~CABFProxy_dUdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFProxy_dUdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CABFProxy_dUdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CABFProxy_dUdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double value = 0.0;
    double temp  = Accu->GetTemperature();

    switch(RealmID){
    // -------------------
        case(ABF_dH): {
            double  nsamples    = Accu->GetData("NTDS",ibin);
            double  micfp       = Accu->GetData("MICFP",ibin,icv);
            double  m2icfp      = Accu->GetData("M2ICFP",ibin,icv);

            double  chp         = Accu->GetData("C11PP",ibin,icv) / nsamples;
            double  m2eint      = Accu->GetData("M2EINT",ibin);

            if( nsamples <= 0 ) return(value);

            double value = micfp - chp / (temp * PMF_Rgas);
            double sicfp = sqrt(m2icfp / nsamples);
            double shp  = sqrt(m2icfp / nsamples) * sqrt( m2eint / nsamples )  / (temp * PMF_Rgas);

            // approximation
            double sigma = sqrt( sicfp*sicfp + shp*shp );

            switch(realm){
                // -------------------
                case(E_PROXY_MEAN):
                    return( value );
                // -------------------
                case(E_PROXY_SD):
                    return( sigma );
                // -------------------
                case(E_PROXY_SEM):
                    return( sigma / sqrt(nsamples) );
                // -------------------
                default:
                    RUNTIME_ERROR("unsupported realm");
            }
        }
        break;
    // -------------------
        case(ABF_MICFP): {
            double  nsamples = Accu->GetData("NTDS",ibin);
            double  micf     = Accu->GetData("MICFP",ibin,icv);
            double  m2icf    = Accu->GetData("M2ICFP",ibin,icv);

            if( nsamples <= 0 ) return(value);

            switch(realm){
                // -------------------
                case(E_PROXY_MEAN):
                    return( micf );
                // -------------------
                case(E_PROXY_SD):
                    return( sqrt(m2icf / nsamples) );
                // -------------------
                case(E_PROXY_SEM):
                    return( sqrt(m2icf) / nsamples );
                // -------------------
                default:
                    RUNTIME_ERROR("unsupported realm");
            }
        }
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
