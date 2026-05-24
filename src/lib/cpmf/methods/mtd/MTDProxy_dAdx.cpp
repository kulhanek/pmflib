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

#include <MTDProxy_dAdx.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDProxy_dAdx::CMTDProxy_dAdx(void)
{
    RegisterRealm(MTD_MICF, "dA/dx", "MTD", "dA(x)=|MICF dx|");
}

//------------------------------------------------------------------------------

CMTDProxy_dAdx::~CMTDProxy_dAdx(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDProxy_dAdx::IsWTMeta(void)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    return( Accu->HasSectionData("MTD-WT") );
}

//------------------------------------------------------------------------------

double CMTDProxy_dAdx::GetValue(int ibin,int icv,EProxyRealm realm) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    double micf = 0.0;
    double fact = 1.0;

    if( Accu->HasSectionData("MTD-WT") ){
        // well-tempered metadynamics
        double temp = Accu->GetTemperature();
        double wtem = Accu->GetData("MTD-WT",0);
        if( wtem > 0 ){
            fact = (temp + wtem) / wtem;
        } else {
            RUNTIME_ERROR("MTD-WT temerature is not greater than zero");
        }
    }

    switch(RealmID){
    // -------------------
        case(MTD_MICF):
            micf     = Accu->GetData("APLFORCE",ibin,icv);
        break;
    // -------------------
        default:
            RUNTIME_ERROR("unsupported type");
    }

    switch(realm){
// mean force
        // -------------------
        case(E_PROXY_VALUE):
            return( micf * fact );
        // -------------------
        case(E_PROXY_SIGMA):
            return( 0.0 );
        // -------------------
        case(E_PROXY_ERROR):
            return( 0.0 );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(0.0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



