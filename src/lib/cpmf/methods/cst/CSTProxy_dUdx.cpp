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

#include <CSTProxy_dUdx.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dUdx::CCSTProxy_dUdx(void)
{
//    SetType(CST_dH);
//
//    Requires.push_back("CST");

//        case(CST_dH):
//            return("CST dH(x) (based on derivatives)");

//    RegisterRealm(CST_dLAMBDAdx,    "dLAMBDA/dx",   "CST", "|<lambda> dx|");
//    RegisterRealm(CST_dMICFdx,      "dMICF/dx",     "CST", "|<lambda> dx|");
//    RegisterRealm(CST_dMICFPFWdx,   "dMICFPFW/dx",  "CST", "|<lambda> dx|");
//    RegisterRealm(CST_dMICFKFWdx,   "dMICFKFW/dx",  "CST", "|<lambda> dx|");

}

//------------------------------------------------------------------------------

CCSTProxy_dUdx::~CCSTProxy_dUdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dUdx::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dUdx::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_dUdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double value = 0.0;
    double ncorr = Accu->GetNCorr();
    double temp  = Accu->GetTemperature();

    switch(RealmID){
    // -------------------
        case(CST_dH): {
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
                case(E_PROXY_VALUE):
                    return( value );
                // -------------------
                case(E_PROXY_SIGMA):
                    return( sigma );
                // -------------------
                case(E_PROXY_ERROR):
                    return( sqrt(ncorr) * sigma / sqrt(nsamples) );
                // -------------------
                default:
                    RUNTIME_ERROR("unsupported realm");
            }
        }
        break;
//       case(CST_dGdx):       // this requires MTC correction
//        case(CST_dLAMBDAdx): {
//            mean        = Accu->GetData("MLAMBDA",ibin,icv);
//            double M2   = Accu->GetData("M2LAMBDA",ibin,icv);
//            samvar      = M2 / nsamples;
//            meanvar     = samvar / nsamples;
//        }
//        break;
//    // ------------------
//        case(CST_dMICFdx): {
//            mean        = Accu->GetData("MICF",ibin,icv);
//            double M2   = Accu->GetData("M2ICF",ibin,icv);
//            samvar      = M2 / nsamples;
//            meanvar     = samvar / nsamples;
//        }
//        break;
//    // -------------------
//        case(CST_dMICFFWdx): {
//            double fwsum    = Accu->GetData("FWSUM",ibin);
//            double fwsum2   = Accu->GetData("FWSUM2",ibin);
//
//            mean            = Accu->GetData("MICFFW",ibin,icv);
//            double M2       = Accu->GetData("M2ICFFW",ibin,icv);
//
//            // https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf
//
//            // number of effective measurements
//            double neff = fwsum2 / (fwsum * fwsum);
//
//            // unbiased weighted sample variance
//            samvar          = M2 / fwsum * neff / (neff - 1.0);
//
//            // variance of the weighted mean
//            // unbiased importance weights
//            meanvar         = samvar / neff;
//        }
//        break;
//    // -------------------
//        case(CST_dMICFPFWdx):  {
//            double fwsum    = Accu->GetData("FWSUM",ibin);
//            double fwsum2   = Accu->GetData("FWSUM2",ibin);
//
//            mean            = Accu->GetData("MICFPFW",ibin,icv);
//            double M2       = Accu->GetData("M2ICFPFW",ibin,icv);
//
//            // number of effective measurements
//            double neff = fwsum2 / (fwsum * fwsum);
//
//            // unbiased weighted sample variance
//            samvar          = M2 / fwsum * neff / (neff - 1.0);
//
//            // variance of the weighted mean
//            // unbiased importance weights
//            meanvar         = samvar / neff;
//        }
//        break;
//    // -------------------
//        case(CST_dMICFKFWdx):  {
//            double fwsum    = Accu->GetData("FWSUM",ibin);
//            double fwsum2   = Accu->GetData("FWSUM2",ibin);
//
//            mean            = Accu->GetData("MICFKFW",ibin,icv);
//            double M2       = Accu->GetData("M2ICFKFW",ibin,icv);
//
//            // number of effective measurements
//            double neff = fwsum2 / (fwsum * fwsum);
//
//            // unbiased weighted sample variance
//            samvar          = M2 / fwsum * neff / (neff - 1.0);
//
//            // variance of the weighted mean
//            // unbiased importance weights
//            meanvar         = samvar / neff;
//        }
//        break;
//    // -------------------
//        case(CST_MICFP): {
//            double  nsamples = Accu->GetData("NTDS",ibin);
//            double  micf     = Accu->GetData("MICFP",ibin,icv);
//            double  m2icf    = Accu->GetData("M2ICFP",ibin,icv);
//
//            if( nsamples <= 0 ) return(value);
//
//            switch(realm){
//                // -------------------
//                case(E_PROXY_VALUE):
//                    return( micf );
//                // -------------------
//                case(E_PROXY_SIGMA):
//                    return( sqrt(m2icf / nsamples) );
//                // -------------------
//                case(E_PROXY_ERROR):
//                    return( sqrt(m2icf * ncorr) / nsamples );
//                // -------------------
//                default:
//                    RUNTIME_ERROR("unsupported realm");
//            }
//        }
//        break;
//    // -------------------
//        case(CST_MICFPFW): {
//            double  nsamples = Accu->GetData("NTDS",ibin);
//
//            double  mfixmanw    = Accu->GetData("MFW",ibin,icv);
////            double  m2fixmanw   = Accu->GetData("M2FIXW",ibin,icv);
//
//            double  micfz       = Accu->GetData("MICFPFW",ibin,icv);
//      //      double  m2icfz      = Accu->GetData("M2ICFPZ",ibin,icv);
//
//            if( nsamples <= 0 ) return(value);
//
//            switch(realm){
//                // -------------------
//                case(E_PROXY_VALUE):
//                    return( micfz / mfixmanw );
//                // -------------------
//                case(E_PROXY_SIGMA):
//                    return(0.0); // FIXME
//                    // return( sqrt(m2icf / nsamples) );
//                // -------------------
//                case(E_PROXY_ERROR):
//                    return(0.0); // FIXME
//                    // return( sqrt(m2icf * ncorr) / nsamples );
//                // -------------------
//                default:
//                    RUNTIME_ERROR("unsupported realm");
//            }
//        }
//        break;
//    // -------------------
//        case(CST_MICFK): {
//            double  nsamples = Accu->GetData("NTDS",ibin);
//            double  micf     = Accu->GetData("MICFK",ibin,icv);
//            double  m2icf    = Accu->GetData("M2ICFK",ibin,icv);
//
//            if( nsamples <= 0 ) return(value);
//
//            switch(realm){
//                // -------------------
//                case(E_PROXY_VALUE):
//                    return( micf );
//                // -------------------
//                case(E_PROXY_SIGMA):
//                    return( sqrt(m2icf / nsamples) );
//                // -------------------
//                case(E_PROXY_ERROR):
//                    return( sqrt(m2icf * ncorr) / nsamples );
//                // -------------------
//                default:
//                    RUNTIME_ERROR("unsupported realm");
//            }
//        }
//        break;
//    // -------------------
//        case(CST_MICFKFW): {
//            double  nsamples = Accu->GetData("NTDS",ibin);
//
//            double  mfixmanw    = Accu->GetData("MFW",ibin,icv);
////            double  m2fixmanw   = Accu->GetData("M2FIXW",ibin,icv);
//
//            double  micfz       = Accu->GetData("MICFKFW",ibin,icv);
//      //      double  m2icfz      = Accu->GetData("M2ICFPZ",ibin,icv);
//
//            if( nsamples <= 0 ) return(value);
//
//            switch(realm){
//                // -------------------
//                case(E_PROXY_VALUE):
//                    return( micfz / mfixmanw );
//                // -------------------
//                case(E_PROXY_SIGMA):
//                    return(0.0); // FIXME
//                    // return( sqrt(m2icf / nsamples) );
//                // -------------------
//                case(E_PROXY_ERROR):
//                    return(0.0); // FIXME
//                    // return( sqrt(m2icf * ncorr) / nsamples );
//                // -------------------
//                default:
//                    RUNTIME_ERROR("unsupported realm");
//            }
//        }
//        break;
//    // -------------------
//        case(CST_C11PP): {
//            double  nsamples    = Accu->GetData("NTDS",ibin);
//            double  fw          = Accu->GetData("MFW",ibin,icv);
//            double  micfpeintfw = Accu->GetData("MICFPEINTFW",ibin,icv);
//            double  micfpfw     = Accu->GetData("MICFPFW",ibin,icv);
//            double  meintfw     = Accu->GetData("MEINTFW",ibin);
//
//            if( nsamples <= 0 ) return(value);
//
//            double value = (micfpeintfw/fw - micfpfw*meintfw/(fw*fw)) / (temp * PMF_Rgas);
//
//            // FIXME
//            double sigma = 0.0;
//
//            switch(realm){
//                // -------------------
//                case(E_PROXY_VALUE):
//                    return( value );
//                // -------------------
//                case(E_PROXY_SIGMA):
//                    return( sigma );
//                // -------------------
//                case(E_PROXY_ERROR):
//                    return( sqrt(ncorr) * sigma / sqrt(nsamples) );
//                // -------------------
//                default:
//                    RUNTIME_ERROR("unsupported realm");
//            }
//        }
//        break;
//// -------------------
//        case(CST_C11PPFW): {
//            double  nsamples    = Accu->GetData("NTDS",ibin);
//            double  m2icfp      = Accu->GetData("M2ICFP",ibin,icv);
//
//            double  chp         = Accu->GetData("C11PP",ibin,icv) / nsamples;
//            double  m2eint      = Accu->GetData("M2EINT",ibin);
//
//            if( nsamples <= 0 ) return(value);
//
//            double value = chp / (temp * PMF_Rgas);
//  //          double sicfp = sqrt(m2icfp / nsamples);
//            double shp  = sqrt(m2icfp / nsamples) * sqrt( m2eint / nsamples )  / (temp * PMF_Rgas);
//
//            // approximation
//            double sigma = sqrt(  shp*shp );
//
//            switch(realm){
//                // -------------------
//                case(E_PROXY_VALUE):
//                    return( value );
//                // -------------------
//                case(E_PROXY_SIGMA):
//                    return( sigma );
//                // -------------------
//                case(E_PROXY_ERROR):
//                    return( sqrt(ncorr) * sigma / sqrt(nsamples) );
//                // -------------------
//                default:
//                    RUNTIME_ERROR("unsupported realm");
//            }
//        }
//        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
