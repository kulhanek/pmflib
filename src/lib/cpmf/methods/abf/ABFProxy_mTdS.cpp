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
    Requires.push_back("ABF");
}

//------------------------------------------------------------------------------

CABFProxy_mTdS::~CABFProxy_mTdS(void)
{
}

//------------------------------------------------------------------------------

bool CABFProxy_mTdS::IsCompatible(CPMFAccumulatorPtr accu)
{
    if( accu->GetMethod() == "ABF" ) return(true);
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
        case(ABF_TdS_HP):
            Provide = "ABF -TdS(x) - cov(dH/dx,Epot)";
        break;
    // -------------------
        case(ABF_TdS_HR):
            Provide = "ABF -TdS(x) - cov(dH/dx,Erst)";
        break;
    // -------------------
        case(ABF_TdS_HK):
            Provide = "ABF -TdS(x) - cov(dH/dx,Ekin)";
        break;
    // -------------------
        case(ABF_TdS_HV):
            Provide = "ABF -TdS(x) - cov(dH/dx,pV)";
        break;

    // -------------------
        case(ABF_TdS_BP):
            Provide = "ABF -TdS(x) - cov(bias,Epot)";
        break;
    // -------------------
        case(ABF_TdS_BR):
            Provide = "ABF -TdS(x) - cov(bias,Erst)";
        break;
    // -------------------
        case(ABF_TdS_BK):
            Provide = "ABF -TdS(x) -cov(bias,Ekin)";
        break;
    // -------------------
        case(ABF_TdS_BV):
            Provide = "ABF -TdS(x) -cov(bias,pV)";
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

    double  nsamples = Accu->GetData("NTDS",ibin);
    return(nsamples);
}

//------------------------------------------------------------------------------

void CABFProxy_mTdS::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CABFProxy_mTdS::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = Accu->GetData("NTDS",ibin);
    double  ncorr    = Accu->GetNCorr();
    double  temp     = Accu->GetTemperature();
    double  value    = 0.0;


    double  c11     = 0.0;
    double  m2icf   = 0.0;
    double  m2ene   = 0.0;

    if( nsamples <= 0 ) return(value);

    switch(Type){
    // -------------------
        case(ABF_TdS_HH):{
            double m2pp = Accu->GetData("M2PP",ibin,icv);
            double m2pn = Accu->GetData("M2PN",ibin,icv);
            c11 = 0.25*(m2pp-m2pn)/nsamples;
            m2icf   = Accu->GetData("M2HICF",ibin,icv);
            m2ene   = Accu->GetData("M2ETOT",ibin);
        }
        break;

    // -------------------
        case(ABF_TdS_HP):
            c11     = Accu->GetData("C11HP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2HICF",ibin,icv);
            m2ene   = Accu->GetData("M2EPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_HR):
            c11     = Accu->GetData("C11HR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2HICF",ibin,icv);
            m2ene   = Accu->GetData("M2ERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_HK):
            c11     = Accu->GetData("C11HK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2HICF",ibin,icv);
            m2ene   = Accu->GetData("M2EKIN",ibin);
        break;
    // -------------------
        case(ABF_TdS_HV):
            c11     = Accu->GetData("C11HV",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2HICF",ibin,icv);
            m2ene   = Accu->GetData("M2EPV",ibin);
        break;
    // -------------------
        case(ABF_TdS_BP):
            c11     = Accu->GetData("C11BP",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2BICF",ibin,icv);
            m2ene   = Accu->GetData("M2EPOT",ibin);
        break;
    // -------------------
        case(ABF_TdS_BR):
            c11     = Accu->GetData("C11BR",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2BICF",ibin,icv);
            m2ene   = Accu->GetData("M2ERST",ibin);
        break;
    // -------------------
        case(ABF_TdS_BK):
            c11     = Accu->GetData("C11BK",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2BICF",ibin,icv);
            m2ene   = Accu->GetData("M2EKIN",ibin);
        break;
    // -------------------
        case(ABF_TdS_BV):
            c11     = Accu->GetData("C11BV",ibin,icv)/nsamples;
            m2icf   = Accu->GetData("M2BICF",ibin,icv);
            m2ene   = Accu->GetData("M2EPV",ibin);
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



