// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <ABFProxy_dU.hpp>
#include <CSTProxy_dU.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergyProxyInit::InitProxyList(std::list<CEnergyProxyPtr>& ene_proxies)
{
    CEnergyProxyPtr proxy;

// add supported proxies
    proxy = CEnergyProxyPtr(new CABFProxy_dU);
    ene_proxies.push_back(proxy);

    proxy = CEnergyProxyPtr(new CCSTProxy_dU);
    ene_proxies.push_back(proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyProxyPtr CEnergyProxyInit::InitProxy(const CSmallString& realm,CPMFAccumulatorPtr& accu)
{
    std::list<CEnergyProxyPtr> ened_proxies;
    InitProxyList(ened_proxies);

    CEnergyProxyPtr proxy;

// find suitable proxy
    std::list<CEnergyProxyPtr>::iterator it = ened_proxies.begin();
    std::list<CEnergyProxyPtr>::iterator ie = ened_proxies.end();

    while( it != ie ){
        proxy = *it;
        it++;
        if( proxy->IsCompatible(accu) == false ) continue;
        if( proxy->SetType(realm) ) return(proxy);
    }

    CSmallString error;
    error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
    RUNTIME_ERROR(error);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
