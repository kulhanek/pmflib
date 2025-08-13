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
#include <iomanip>
#include <ABFProxy_dU.hpp>
#include <CSTProxy_dU.hpp>
#include <CSTProxy_Ecorr.hpp>

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

    proxy = CEnergyProxyPtr(new CCSTProxy_Ecorr);
    ene_proxies.push_back(proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyProxyPtr CEnergyProxyInit::InitProxy(const CSmallString& realm,CPMFAccumulatorPtr& accu)
{
    std::list<CEnergyProxyPtr> ene_proxies;
    InitProxyList(ene_proxies);

    CEnergyProxyPtr proxy;

// find suitable proxy
    std::list<CEnergyProxyPtr>::iterator it = ene_proxies.begin();
    std::list<CEnergyProxyPtr>::iterator ie = ene_proxies.end();

    while( it != ie ){
        proxy = *it;
        it++;
        if( proxy->IsCompatible(accu) == false ) continue;
        if( proxy->SetRealm(realm) ) return(proxy);
    }

    CSmallString error;
    error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
    RUNTIME_ERROR(error);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergyProxyInit::EnumerateRealms(std::list<CProxyRealmDescr>& dlist)
{
    std::list<CEnergyProxyPtr> ene_proxies;
    InitProxyList(ene_proxies);

    CEnergyProxyPtr proxy;

// find suitable proxy
    std::list<CEnergyProxyPtr>::iterator it = ene_proxies.begin();
    std::list<CEnergyProxyPtr>::iterator ie = ene_proxies.end();

    while( it != ie ){
        proxy = *it;
        proxy->EnumerateRealms(dlist);
        it++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergyProxyInit::PrintRealms(std::ostream& fout)
{
    std::list<CProxyRealmDescr> dlist;
    EnumerateRealms(dlist);
    dlist.sort(CProxyRealmDescr::Compare);

    std::list<CProxyRealmDescr>::iterator it = dlist.begin();
    std::list<CProxyRealmDescr>::iterator ie = dlist.end();

    fout << std::endl;
    fout << "# Realm              Method Description                                           " << std::endl;
    fout << "# ------------------ ------ ------------------------------------------------------" << std::endl;

    while( it != ie ){
        fout << std::left << std::setw(20) << (*it).Realm << " " << std::setw(6) << (*it).Method << " " << (*it).Description << std::endl;
        it++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
