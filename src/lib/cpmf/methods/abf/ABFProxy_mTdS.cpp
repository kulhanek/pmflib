// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <ABFProxy_mTdS.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFProxy_mTdS::CABFProxy_mTdS(void)
{
    SetType(ABF_TdS_HH);
    Requires.push_back("TABF");
    Requires.push_back("ABF");
    Requires.push_back("US-ABF");
}

//------------------------------------------------------------------------------

CABFProxy_mTdS::~CABFProxy_mTdS(void)
{
}

//------------------------------------------------------------------------------

bool CABFProxy_mTdS::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "ABF" ) return(true);
    if( accu->GetMethod() == "TABF" ) return(true);
    if( accu->GetMethod() == "US-ABF" ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

void CABFProxy_mTdS::SetType(EABFTdSType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(ABF_TdS_HH):
            Provide = "ABF -TdS(x)";
        break;

    // -------------------
        case(ABF_TdS_HH_MTC):
            Provide = "ABF -TdS(x) MTC";
        break;

    // -------------------
        case(ABF_TdS_FP):
            Provide = "ABF -TdS(x) - FP";
        break;
    // -------------------
        case(ABF_TdS_FR):
            Provide = "ABF -TdS(x) - FR";
        break;
    // -------------------
        case(ABF_TdS_FK):
            Provide = "ABF -TdS(x) - FK";
        break;

    // -------------------
        case(ABF_TdS_VP):
            Provide = "ABF -TdS(x) - VP";
        break;
    // -------------------
        case(ABF_TdS_VR):
            Provide = "ABF -TdS(x) - VR";
        break;
    // -------------------
        case(ABF_TdS_VK):
            Provide = "ABF -TdS(x) - VK";
        break;

    // -------------------
        case(ABF_TdS_BP):
            Provide = "ABF -TdS(x) - BP";
        break;
    // -------------------
        case(ABF_TdS_BR):
            Provide = "ABF -TdS(x) - BR";
        break;
    // -------------------
        case(ABF_TdS_BK):
            Provide = "ABF -TdS(x) - BK";
        break;

    // -------------------
        case(ABF_TdS_SP):
            Provide = "ABF -TdS(x) - SP";
        break;
    // -------------------
        case(ABF_TdS_SR):
            Provide = "ABF -TdS(x) - SR";
        break;
    // -------------------
        case(ABF_TdS_SK):
            Provide = "ABF -TdS(x) - SK";
        break;

    // -------------------
        case(ABF_TdS_LP):
            Provide = "ABF -TdS(x) - LP";
        break;
    // -------------------
        case(ABF_TdS_LR):
            Provide = "ABF -TdS(x) - LR";
        break;
    // -------------------
        case(ABF_TdS_LK):
            Provide = "ABF -TdS(x) - LK";
        break;

    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFProxy_mTdS::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

void CABFProxy_mTdS::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NSAMPLES",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CABFProxy_mTdS::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = 0.0;

    if( Type == ABF_TdS_HH ) {
        nsamples = Accu->GetData("NSAMPLES",ibin);
    } else {
        nsamples = Accu->GetData("NTDS",ibin);
    }

    double  ncorr    = Accu->GetNCorr();
    double  temp     = Accu->GetTemperature();

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    double  c11     = 0.0;
    double  m2icf   = 0.0;
    double  m2ene   = 0.0;

    switch(Type){
    // -------------------
        case(ABF_TdS_HH):{
            double m2pp = Accu->GetData("M2PP",ibin,icv);
            double m2pn = Accu->GetData("M2PN",ibin,icv);
            c11 = 0.25*(m2pp-m2pn)/nsamples;
            m2icf   = Accu->GetData("M2ICF",ibin,icv);
            m2ene   = Accu->GetData("M2ETOT",ibin);
        }
        break;

    // -------------------
        case(ABF_TdS_HH_MTC):{
            double micfetot = Accu->GetData("MICFETOT_MTC",ibin,icv);
            double micf     = Accu->GetData("MICF_MTC",ibin,icv);
            double metot    = Accu->GetData("METOT_MTC",ibin);
            double mmtc     = 1.0;
            if( Accu->HasSectionData("MMTC") ){
                mmtc     = Accu->GetData("MMTC",ibin);
            }

            c11     = micfetot/mmtc - metot*micf/(mmtc*mmtc);
            m2icf   = Accu->GetData("M2ICF_MTC",ibin,icv);
            m2ene   = Accu->GetData("M2ETOT_MTC",ibin);
        }
        break;

    // -------------------
        case(ABF_TdS_FP):
            c11     = Accu->GetData("C11TDSFP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSFX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_FR):
            c11     = Accu->GetData("C11TDSFR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSFX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_FK):
            c11     = Accu->GetData("C11TDSFK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSFX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEKIN",ibin);
        break;

    // -------------------
        case(ABF_TdS_VP):
            c11     = Accu->GetData("C11TDSVP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSVX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_VR):
            c11     = Accu->GetData("C11TDSVR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSVX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_VK):
            c11     = Accu->GetData("C11TDSVK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSVX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEKIN",ibin);
        break;

    // -------------------
        case(ABF_TdS_BP):
            c11     = Accu->GetData("C11TDSBP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSBX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_BR):
            c11     = Accu->GetData("C11TDSBR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSBX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_BK):
            c11     = Accu->GetData("C11TDSBK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSBX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEKIN",ibin);
        break;

    // -------------------
        case(ABF_TdS_SP):
            c11     = Accu->GetData("C11TDSSP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSSX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_SR):
            c11     = Accu->GetData("C11TDSSR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSSX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_SK):
            c11     = Accu->GetData("C11TDSSK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSSX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEKIN",ibin);
        break;

                    // -------------------
        case(ABF_TdS_LP):
            c11     = Accu->GetData("C11TDSLP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSLX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_LR):
            c11     = Accu->GetData("C11TDSLR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSLX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_LK):
            c11     = Accu->GetData("C11TDSLK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2TDSLX",ibin,icv);
            m2ene   = Accu->GetData("M2TDSEKIN",ibin);
        break;

    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    switch(realm){
    // -------------------
        case(E_PROXY_VALUE): {
            return( c11  / (temp * PMF_Rgas) );
        }
    // -------------------
        case(E_PROXY_SIGMA): {
            // approximation
            return( sqrt(m2icf / nsamples) * sqrt( m2ene / nsamples )  / (temp * PMF_Rgas) );
        }
    // -------------------
        case(E_PROXY_ERROR): {
            // approximation
            return( sqrt(ncorr) * sqrt(m2icf / nsamples) * sqrt( m2ene / nsamples ) / sqrt(nsamples) / (temp * PMF_Rgas) );
        }
    // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



