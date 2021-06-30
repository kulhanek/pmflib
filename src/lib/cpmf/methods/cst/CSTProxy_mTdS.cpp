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

#include <CSTProxy_mTdS.hpp>
#include <PMFConstants.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_mTdS::CCSTProxy_mTdS(void)
{
    Requires.push_back("CST");
    Provide = "CST -TdS(x)^{c}";    // entropy of the constrained system
    Type = CST_C11HH;
}

//------------------------------------------------------------------------------

CCSTProxy_mTdS::~CCSTProxy_mTdS(void)
{
}

//------------------------------------------------------------------------------

void CCSTProxy_mTdS::SetType(ECSTTdSType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(CST_C11HH):
            Provide = "CST -TdS(x)^{c}";    // entropy of the constrained system
        break;
    // -------------------
        case(CST_C11HP):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Epot)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        case(CST_C11HK):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Ekin)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        case(CST_C11HR):
            Provide = "CST -TdS(x)^{c} cov(dH/dx,Erst)";    // entropy of the constrained system  - contribution
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_mTdS::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_mTdS::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NSAMPLES",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_mTdS::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples = Accu->GetData("NSAMPLES",ibin);
    double  m2lam    = Accu->GetData("M2LAMBDA",ibin,icv);
    double  ncorr    = Accu->GetNCorr();

    double  c11     = 0.0;
    double  m2ene   = 0.0;

    switch(Type){
    // -------------------
        case(CST_C11HH):
            c11     = Accu->GetData("C11HH",ibin,icv);
            m2ene   = Accu->GetData("M2ETOT",ibin);
        break;
    // -------------------
        case(CST_C11HP):
            c11     = Accu->GetData("C11HP",ibin,icv);
            m2ene   = Accu->GetData("M2EPOT",ibin);
        break;
    // -------------------
        case(CST_C11HK):
            c11     = Accu->GetData("C11HK",ibin,icv);
            m2ene   = Accu->GetData("M2EKIN",ibin);
        break;
    // -------------------
        case(CST_C11HR):
            c11     = Accu->GetData("C11HR",ibin,icv);
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
            // negative value due to lambda vs dG/dx
            return( - (c11 / nsamples) / (temp * PMF_Rgas) );
        }
    // -------------------
        case(E_PROXY_SIGMA): {
            // approximation
            return( sqrt(m2lam / nsamples) * sqrt( m2ene / nsamples )  / (temp * PMF_Rgas) );
        }
    // -------------------
        case(E_PROXY_ERROR): {
            // approximation
            return( sqrt(ncorr) * sqrt(m2lam / nsamples) * sqrt( m2ene / nsamples ) / sqrt(nsamples) / (temp * PMF_Rgas) );
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



