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
    double  nsamples = 0.0;
    double  micf     = 0.0;
    double  m2icf    = 0.0;

    switch(Type){
    // -------------------
        case(CST_dG): {
            nsamples = Accu->GetData("NSAMPLES",ibin);
            micf     = Accu->GetData("MLAMBDA",ibin,icv);
            m2icf    = Accu->GetData("M2LAMBDA",ibin,icv);
            // this requires MTC correction
        }
        break;
    // ------------------
        case(CST_MICF): {
            nsamples = Accu->GetData("NTDS",ibin);
            micf     = Accu->GetData("MICF",ibin,icv);
            m2icf    = Accu->GetData("M2ICF",ibin,icv);
        }
        break;
    // -------------------
        case(CST_MICFFW): {
            double mfw  = Accu->GetData("MFWTDS",ibin);
            double m2fw = Accu->GetData("M2FWTDS",ibin);

            nsamples    = Accu->GetData("NTDS",ibin);
            double mup  = Accu->GetData("MICFFW",ibin,icv);
            double m2up = Accu->GetData("M2ICFFW",ibin,icv);

            micf  = mup / mfw;
            m2icf = micf*micf * ( m2up/(mup*mup) + m2fw/(mfw*mfw) );
        }
        break;
    // -------------------
        case(CST_MICFPFW):  {
            double mfw  = Accu->GetData("MFWTDS",ibin);
            double m2fw = Accu->GetData("M2FWTDS",ibin);

            nsamples    = Accu->GetData("NTDS",ibin);
            double mup  = Accu->GetData("MICFPFW",ibin,icv);
            double m2up = Accu->GetData("M2ICFPFW",ibin,icv);

            micf  = mup / mfw;
            m2icf = micf*micf * ( m2up/(mup*mup) + m2fw/(mfw*mfw) );
        }
        break;
    // -------------------
        case(CST_MICFKFW):  {
            double mfw  = Accu->GetData("MFWTDS",ibin);
            double m2fw = Accu->GetData("M2FWTDS",ibin);

            nsamples    = Accu->GetData("NTDS",ibin);
            double mup  = Accu->GetData("MICFKFW",ibin,icv);
            double m2up = Accu->GetData("M2ICFKFW",ibin,icv);

            micf  = mup / mfw;
            m2icf = micf*micf * ( m2up/(mup*mup) + m2fw/(mfw*mfw) );
        }
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
            return( micf );
        // -------------------
        case(E_PROXY_SIGMA):
            return( sqrt(m2icf / nsamples) );
        // -------------------
        case(E_PROXY_ERROR):
            return( sqrt(m2icf * ncorr) / nsamples );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



