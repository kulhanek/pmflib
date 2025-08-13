#ifndef ABFProxy_dG_H
#define ABFProxy_dG_H
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
#include <EnergyDerProxy.hpp>

//------------------------------------------------------------------------------

enum EABFdGType {
    ABF_MICF,
};

//------------------------------------------------------------------------------

/** \brief ABF proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CABFProxy_dGdx : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CABFProxy_dGdx(void);
    ~CABFProxy_dGdx(void);

//------------------------------------------------------------------------------
    // get energy derivative and its error
    virtual double GetValue( int ibin,int icv,EProxyRealm realm) const;

};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CABFProxy_dGdx>    CABFProxy_dGdx_Ptr;

//------------------------------------------------------------------------------

#endif
