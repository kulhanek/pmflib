#ifndef EnergyProxyH
#define EnergyProxyH
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
#include <SmallString.hpp>
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

enum EProxyRealm {
    E_PROXY_VALUE        = 1,  // energy
    E_PROXY_SIGMA        = 2,  // fluctuation of energy
    E_PROXY_ERROR        = 3,  // error of energy
};

//------------------------------------------------------------------------------

class PMF_PACKAGE CEnergyProxy {
public:
// constructor and destructor --------------------------------------------------
    CEnergyProxy(void);
    virtual ~CEnergyProxy(void);

// setup methods ---------------------------------------------------------------
    // set accumulator and perform sanity checks
    virtual void Init(CPMFAccumulatorPtr accu);

// access methods --------------------------------------------------------------
    // get PMF accumulator
    CPMFAccumulatorPtr GetAccu(void);

    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy and its error
    virtual double GetValue( int ibin,EProxyRealm realm) const;

    // set NCorr
    virtual void SetNCorr(double ncorr);

// protected data --------------------------------------------------------------
protected:
    CSmallString        Require;
    CSmallString        Provide;
    CPMFAccumulatorPtr  Accu;
};

//------------------------------------------------------------------------------

typedef std::shared_ptr<CEnergyProxy>    CEnergyProxyPtr;

//------------------------------------------------------------------------------

#endif
