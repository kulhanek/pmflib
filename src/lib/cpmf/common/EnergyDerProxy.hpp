#ifndef EnergyDerProxyH
#define EnergyDerProxyH
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

#include <PMFMainHeader.hpp>
#include <BaseProxy.hpp>
#include <EnergyProxy.hpp>

//------------------------------------------------------------------------------

class PMF_PACKAGE CEnergyDerProxy : public CBaseProxy {
public:
// constructor and destructor --------------------------------------------------
    CEnergyDerProxy(void);
    virtual ~CEnergyDerProxy(void);

// access methods -------------------------------------------------------------
    // get optional energy correction
    virtual CEnergyProxyPtr GetEnergyCorrection(void);

    // get derivative and its error
    double GetValue(int ibin,int cv,EProxyRealm realm) const;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CEnergyDerProxy>    CEnergyDerProxyPtr;

//------------------------------------------------------------------------------

#endif
