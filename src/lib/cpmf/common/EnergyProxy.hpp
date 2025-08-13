#ifndef EnergyProxyH
#define EnergyProxyH
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

#include <PMFMainHeader.hpp>
#include <BaseProxy.hpp>

//------------------------------------------------------------------------------

class PMF_PACKAGE CEnergyProxy : public CBaseProxy {
public:
// constructor and destructor --------------------------------------------------
    CEnergyProxy(void);
    virtual ~CEnergyProxy(void);

// access methods --------------------------------------------------------------
    // get energy and its error
    virtual double GetValue( int ibin,EProxyRealm realm) const;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CEnergyProxy>    CEnergyProxyPtr;

//------------------------------------------------------------------------------

#endif
