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

#include <PMFProxy_dH.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFProxy_dH::CPMFProxy_dH(void)
{
    Requires.push_back("ABF");
    Requires.push_back("CST");
    SetType(PMF_EINT);
}

//------------------------------------------------------------------------------

CPMFProxy_dH::~CPMFProxy_dH(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFProxy_dH::SetType(EPMFdHType type)
{
    Type = type;

    switch(Type){
    // -------------------
        case(PMF_EINT):
            Provide = "dH(x)=<Eint>";
        break;
    // -------------------
        case(PMF_ETOT):
            Provide = "dH(x)=<Etot>";
        break;
    // -------------------
        case(PMF_EPOT):
            Provide = "dH(x)=<Epot>";
        break;
    // -------------------
        case(PMF_EKIN):
            Provide = "<Ekin>";
        break;
    // -------------------
        case(PMF_ERST):
            Provide = "<Erst>";
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFProxy_dH::GetVSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

double CPMFProxy_dH::GetValue( int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  nsamples    = Accu->GetData("NTDS",ibin);
    double  mene        = 0.0;
    double  m2ene       = 0.0;
    double  ncorr       = Accu->GetNCorr();

    switch(Type){
    // -------------------
        case(PMF_EINT):
            if( Accu->HasSectionData("MEINT") ){
                mene    = Accu->GetData("MEINT",ibin);
                m2ene   = Accu->GetData("M2EINT",ibin);
            } else {
                mene    = Accu->GetData("MEPOT",ibin) + Accu->GetData("MERST",ibin);
                m2ene    = Accu->GetData("M2ETOT",ibin) + Accu->GetData("M2ERST",ibin);
            }
        break;
    // -------------------
        case(PMF_ETOT):
            mene    = Accu->GetData("METOT",ibin);
            m2ene   = Accu->GetData("M2ETOT",ibin);
        break;
    // -------------------
        case(PMF_EPOT):
            mene    = Accu->GetData("MEPOT",ibin);
            m2ene   = Accu->GetData("M2EPOT",ibin);
        break;
    // -------------------
        case(PMF_EKIN):
            mene    = Accu->GetData("MEKIN",ibin);
            m2ene   = Accu->GetData("M2EKIN",ibin);
        break;
    // -------------------
        case(PMF_ERST):
            mene    = Accu->GetData("MERST",ibin);
            m2ene   = Accu->GetData("M2ERST",ibin);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
// mean force
        // -------------------
        case(E_PROXY_VALUE):
            return( mene );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(m2ene / nsamples) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(m2ene * ncorr) / nsamples );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



