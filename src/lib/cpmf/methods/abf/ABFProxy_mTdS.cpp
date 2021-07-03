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
    SetType(ABF_C11HH);
    Requires.push_back("TABF");
    Requires.push_back("ABF");
}

//------------------------------------------------------------------------------

CABFProxy_mTdS::~CABFProxy_mTdS(void)
{
}

//------------------------------------------------------------------------------

void CABFProxy_mTdS::SetType(EABFTdSType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(ABF_C11HH):
            Provide = "ABF -TdS(x)";
        break;
    // -------------------
        case(ABF_C11PP):
            Provide = "ABF -TdS(x) cov(dH_p/dx,Epot)";
        break;
    // -------------------
        case(ABF_C11PK):
            Provide = "ABF -TdS(x) cov(dH_p/dx,Ekin)";
        break;
    // -------------------
        case(ABF_C11PR):
            Provide = "ABF -TdS(x) cov(dH_p/dx,Erst)";
        break;
    // -------------------
        case(ABF_C11KP):
            Provide = "ABF -TdS(x) cov(dH_k/dx,Epot)";
        break;
    // -------------------
        case(ABF_C11KK):
            Provide = "ABF -TdS(x) cov(dH_k/dx,Ekin)";
        break;
    // -------------------
        case(ABF_C11KR):
            Provide = "ABF -TdS(x) cov(dH_k/dx,Erst)";
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

    double  nsamples = Accu->GetData("NSAMPLES",ibin);
    double  ncorr    = Accu->GetNCorr();

    double  c11     = 0.0;
    double  m2icf   = 0.0;
    double  m2ene   = 0.0;

    switch(Type){
    // -------------------
        case(ABF_C11HH):
            c11     = Accu->GetData("C11HH",ibin,icv);
            m2icf   = Accu->GetData("M2ICF",ibin,icv);
            m2ene   = Accu->GetData("M2ETOT",ibin);
        break;
    // -------------------
        case(ABF_C11PP):
            c11     = Accu->GetData("C11PP",ibin,icv);
            m2icf   = Accu->GetData("M2ICF_POT",ibin,icv);
            m2ene   = Accu->GetData("M2EPOT",ibin);
        break;
    // -------------------
        case(ABF_C11PK):
            c11     = Accu->GetData("C11PK",ibin,icv);
            m2icf   = Accu->GetData("M2ICF_POT",ibin,icv);
            m2ene   = Accu->GetData("M2EKIN",ibin);
        break;
    // -------------------
        case(ABF_C11PR):
            c11     = Accu->GetData("C11PR",ibin,icv);
            m2icf   = Accu->GetData("M2ICF_POT",ibin,icv);
            m2ene   = Accu->GetData("M2ERST",ibin);
        break;
    // -------------------
        case(ABF_C11KP):
            c11     = Accu->GetData("C11KP",ibin,icv);
            m2icf   = Accu->GetData("M2ICF_KIN",ibin,icv);
            m2ene   = Accu->GetData("M2EPOT",ibin);
        break;
    // -------------------
        case(ABF_C11KK):
            c11     = Accu->GetData("C11KK",ibin,icv);
            m2icf   = Accu->GetData("M2ICF_KIN",ibin,icv);
            m2ene   = Accu->GetData("M2EKIN",ibin);
        break;
    // -------------------
        case(ABF_C11KR):
            c11     = Accu->GetData("C11KR",ibin,icv);
            m2icf   = Accu->GetData("M2ICF_KIN",ibin,icv);
            m2ene   = Accu->GetData("M2ERST",ibin);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    double  temp     = Accu->GetTemperature();

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
    // -------------------
        case(E_PROXY_VALUE): {
            return( (c11 / nsamples) / (temp * PMF_Rgas) );
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



