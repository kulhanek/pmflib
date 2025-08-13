#ifndef CSTProxy_dU_H
#define CSTProxy_dU_H
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
#include <EnergyProxy.hpp>

//------------------------------------------------------------------------------

enum ECSTdHType {
    CST_dU,         // ETOT - Fixman weighted
    CST_ETOT,       // ETOT
    CST_ETOTFW,     // ETOT - Fixman weighted
    CST_EINT,       // EPOT+ERST
    CST_EINTFW,     // EPOT+ERST - Fixman weighted
    CST_EPOT,       // EPOT
    CST_EPOTFW,     // EPOT - Fixman weighted
    CST_ERST,       // ERST
    CST_ERSTFW,     // ERST - Fixman weighted
    CST_EKIN,       // EKIN
    CST_EKINFW,     // EKIN - Fixman weighted
};

//------------------------------------------------------------------------------

/** \brief PMF proxy providing enthalpy
*/

class PMF_PACKAGE CCSTProxy_dU : public CEnergyProxy {
public:
// constructor and destructor --------------------------------------------------
    CCSTProxy_dU(void);
    ~CCSTProxy_dU(void);
//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy derivative and its error
    virtual double GetValue( int ibin,EProxyRealm realm) const;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CCSTProxy_dU>    CCSTProxy_dH_Ptr;

//------------------------------------------------------------------------------

#endif
