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
#include <SmallString.hpp>
#include <PMFAccumulator.hpp>
#include <vector>

//------------------------------------------------------------------------------

// if a given realm cannot provide a value, return zero.

enum EProxyRealm {
    E_PROXY_VALUE        = 1,  // sample mean
    E_PROXY_SIGMA        = 2,  // sample variance (fluctuation)
    E_PROXY_ERROR        = 3,  // variance of sample mean (error estimate)
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

    // get realm
    CSmallString GetRealm(void);

    /// return number of cvs
    int GetNumOfCVs(void) const;

    /// return number of bins
    int GetNumOfBins(void) const;

    // get number of samples (NSAMPLES)
    int GetNSamples(int ibin) const;

    // get energy and its error
    virtual double GetValue( int ibin,EProxyRealm realm) const;

// protected data --------------------------------------------------------------
protected:
    std::vector<std::string>    Requires;
    CSmallString                Provide;
    CPMFAccumulatorPtr          Accu;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CEnergyProxy>    CEnergyProxyPtr;

//------------------------------------------------------------------------------

#endif
