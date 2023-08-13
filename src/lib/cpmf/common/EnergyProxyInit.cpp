// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <EnergyProxyInit.hpp>
#include <PMFProxy_dH.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyProxyPtr CEnergyProxyInit::InitProxy(const CSmallString& realm,CPMFAccumulatorPtr& accu)
{
    CEnergyProxyPtr lproxy;

    if( realm == "<Eint>" ){
        CPMFProxy_dH_Ptr proxy = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EINT);
        lproxy = proxy;
    } else if( realm == "<Etot>" ){
        CPMFProxy_dH_Ptr proxy = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_ETOT);
        lproxy = proxy;
// -----------------------------------------------
    } else if ( realm == "<Epot>" ) {
        CPMFProxy_dH_Ptr proxy = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EPOT);
        lproxy = proxy;
// -----------------------------------------------
    } else if ( realm == "<Ekin>" ) {
        CPMFProxy_dH_Ptr proxy = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EKIN);
        lproxy = proxy;
// -----------------------------------------------
    } else if ( realm == "<Erst>" ) {
        CPMFProxy_dH_Ptr proxy = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_ERST);
        lproxy = proxy;
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << realm ;
        RUNTIME_ERROR(error);
    }

    return(lproxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
