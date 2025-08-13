#ifndef EnergyDerProxyInitH
#define EnergyDerProxyInitH
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

#include <PMFMainHeader.hpp>
#include <EnergyDerProxy.hpp>
#include <list>
#include <iostream>

//------------------------------------------------------------------------------

class PMF_PACKAGE CEnergyDerProxyInit {
public:
    /// init requested proxy
    static CEnergyDerProxyPtr InitProxy(const CSmallString& realm,CPMFAccumulatorPtr& accu);

    /// print supported realms
    static void PrintRealms(std::ostream& fout);

    /// enumerate supported realms
    static void EnumerateTypes(std::list<CProxyRealmDescr>& dlist);

private:
    /// create list of all supported proxies
    static void InitProxyList(std::list<CEnergyDerProxyPtr>& eneder_proxies);
};

//------------------------------------------------------------------------------

#endif
