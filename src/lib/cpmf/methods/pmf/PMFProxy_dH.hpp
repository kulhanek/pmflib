#ifndef PMFProxy_dH_H
#define PMFProxy_dH_H
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
#include <EnergyProxy.hpp>

//------------------------------------------------------------------------------

enum EPMFdHType {
    PMF_ETOT,
    PMF_ETOT_MTC,
    PMF_EPOT,
    PMF_EKIN,
    PMF_ERST,
};

//------------------------------------------------------------------------------

/** \brief PMF proxy providing enthalpy
*/

class PMF_PACKAGE CPMFProxy_dH : public CEnergyProxy {
public:
// constructor and destructor --------------------------------------------------
    CPMFProxy_dH(void);
    ~CPMFProxy_dH(void);

//------------------------------------------------------------------------------
    // set type
    void SetType(EPMFdHType type);

//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy derivative and its error
    virtual double GetValue( int ibin,EProxyRealm realm) const;

// section of private data ----------------------------------------------------
private:
    EPMFdHType  Type;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CPMFProxy_dH>    CPMFProxy_dH_Ptr;

//------------------------------------------------------------------------------

#endif
