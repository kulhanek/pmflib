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

#include <EnergyDerProxyInit.hpp>
#include <ABFProxy_dGdx.hpp>
#include <ABFProxy_dUdx.hpp>
#include <ABFProxy_mTdSdx.hpp>
#include <CSTProxy_dGdx.hpp>
#include <CSTProxy_dUdx.hpp>
#include <CSTProxy_mTdSdx.hpp>
#include <iomanip>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergyDerProxyInit::InitProxyList(std::list<CEnergyDerProxyPtr>& eneder_proxies)
{
    CEnergyDerProxyPtr proxy;

// add supported proxies
    proxy = CEnergyDerProxyPtr(new CABFProxy_dGdx);
    eneder_proxies.push_back(proxy);

    proxy = CEnergyDerProxyPtr(new CABFProxy_dUdx);
    eneder_proxies.push_back(proxy);

    proxy = CEnergyDerProxyPtr(new CABFProxy_mTdSdx);
    eneder_proxies.push_back(proxy);

    proxy = CEnergyDerProxyPtr(new CCSTProxy_dGdx);
    eneder_proxies.push_back(proxy);

    proxy = CEnergyDerProxyPtr(new CCSTProxy_dUdx);
    eneder_proxies.push_back(proxy);

    proxy = CEnergyDerProxyPtr(new CCSTProxy_mTdSdx);
    eneder_proxies.push_back(proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyDerProxyPtr CEnergyDerProxyInit::InitProxy(const CSmallString& realm,CPMFAccumulatorPtr& accu)
{
    std::list<CEnergyDerProxyPtr> eneder_proxies;
    InitProxyList(eneder_proxies);

    CEnergyDerProxyPtr proxy;

// find suitable proxy
    std::list<CEnergyDerProxyPtr>::iterator it = eneder_proxies.begin();
    std::list<CEnergyDerProxyPtr>::iterator ie = eneder_proxies.end();

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

void CEnergyDerProxyInit::EnumerateRealms(std::list<CProxyRealmDescr>& dlist)
{
    std::list<CEnergyDerProxyPtr> eneder_proxies;
    InitProxyList(eneder_proxies);

    CEnergyDerProxyPtr proxy;

// find suitable proxy
    std::list<CEnergyDerProxyPtr>::iterator it = eneder_proxies.begin();
    std::list<CEnergyDerProxyPtr>::iterator ie = eneder_proxies.end();

    while( it != ie ){
        proxy = *it;
        proxy->EnumerateRealms(dlist);
        it++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergyDerProxyInit::PrintRealms(std::ostream& fout)
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
