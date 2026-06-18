// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <CSTProxy_dY.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTProxy_dY::CCSTProxy_dY(void)
{
    RegisterRealm(CST_ETOT,     "<Etot>",   "CST", "<Etot>");
    RegisterRealm(CST_ETOTFW,   "<Etot>FW", "CST", "<Etot>_FW");
    RegisterRealm(CST_EINT,     "<Eint>",   "CST", "<Eint>");
    RegisterRealm(CST_EINTFW,   "<Eint>FW", "CST", "<Eint>_FW");
    RegisterRealm(CST_EPOT,     "<Epot>",   "CST", "<Epot>");
    RegisterRealm(CST_EPOTFW,   "<Epot>FW", "CST", "<Epot>_FW");
    RegisterRealm(CST_ERST,     "<Erst>",   "CST", "<Erst>");
    RegisterRealm(CST_ERSTFW,   "<Erst>FW", "CST", "<Erst>_FW");
    RegisterRealm(CST_EKIN,     "<Ekin>",   "CST", "<Ekin>");
    RegisterRealm(CST_EKINFW,   "<Ekin>FW", "CST", "<Ekin>_FW");
}

//------------------------------------------------------------------------------

CCSTProxy_dY::~CCSTProxy_dY(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTProxy_dY::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NTDS",ibin));
}

//------------------------------------------------------------------------------

void CCSTProxy_dY::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NTDS",ibin,nsamples);
}

//------------------------------------------------------------------------------

double CCSTProxy_dY::GetValue( int ibin,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double  value   = 0.0; // sample mean
    double  sd      = 0.0; // unbiased sample standard deviation
    double  sem     = 0.0; // standard error of the sample mean

    switch(RealmID){
    // -------------------
        case(CST_ETOT):
            GetMeanValue("ETOT",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        case(CST_ETOTFW):
            GetWMeanValue("ETOTFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        case(CST_EINT):
            GetMeanValue("EINT",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        case(CST_EINTFW):
            GetWMeanValue("EINTFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        case(CST_EPOT):
            GetMeanValue("EPOT",value,sd,sem,realm==E_PROXY_MEAN,ibin);  
        break;
    // -------------------
        case(CST_EPOTFW):
            GetWMeanValue("EPOTFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        case(CST_ERST):
            GetMeanValue("ERST",value,sd,sem,realm==E_PROXY_MEAN,ibin);  
        break;
    // -------------------
        case(CST_ERSTFW):
            GetWMeanValue("ERSTFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        case(CST_EKIN):
            GetMeanValue("EKIN",value,sd,sem,realm==E_PROXY_MEAN,ibin); 
        break;
    // -------------------
        case(CST_EKINFW):
            GetWMeanValue("EKINFW","FW",value,sd,sem,realm==E_PROXY_MEAN,ibin);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
        break;
    }

// return result
    switch(realm){
        // -------------------
        case(E_PROXY_MEAN):
            return( value );
        // -------------------
        case(E_PROXY_SD):
            return( sd );
        // -------------------
        case(E_PROXY_SEM):
            return( sem );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



