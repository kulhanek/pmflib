#ifndef BaseProxyH
#define BaseProxyH
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
#include <PMFAccumulator.hpp>
#include <SmallString.hpp>
#include <set>
#include <map>

//------------------------------------------------------------------------------

// if a given realm cannot provide a value, return zero.

enum EProxyRealm {
    E_PROXY_VALUE        = 1,  // sample mean
    E_PROXY_SIGMA        = 2,  // sample variance (fluctuation)
    E_PROXY_ERROR        = 3,  // variance of sample mean (error estimate)
};

//------------------------------------------------------------------------------

class CProxyRealmDescr {
// constructor and destructor --------------------------------------------------
public:
    CProxyRealmDescr(void);

// public data -----------------------------------------------------------------
    CSmallString    Method;
    CSmallString    Realm;
    int             RealmID;
    CSmallString    Description;

public:
    static bool Compare(const CProxyRealmDescr& left,const CProxyRealmDescr& right);
};

//------------------------------------------------------------------------------

class PMF_PACKAGE CBaseProxy {
public:
// constructor and destructor --------------------------------------------------
    CBaseProxy(void);
    virtual ~CBaseProxy(void);

// setup methods ---------------------------------------------------------------
    // register realm
    void RegisterRealm(int realmid,const CSmallString& realm,const CSmallString& method,const CSmallString& descr);

    // set accumulator and perform sanity checks
    virtual void Init(CPMFAccumulatorPtr accu);

    // is compatible with PMF Accumulator method
    virtual bool IsCompatible(CPMFAccumulatorPtr accu);

    // set type if it is supported
    virtual bool SetRealm(const CSmallString& realm);

    // set type - if not supported throw runtime_error
    virtual void SetRealm(int realmid);

    // enumerate types
    virtual void EnumerateRealms(std::list<CProxyRealmDescr>& dlist);

// access methods -------------------------------------------------------------
    // get PMF accumulator
    CPMFAccumulatorPtr GetAccu(void);

    // get realm
    CSmallString GetRealm(void);

    // get methods
    CSmallString GetMethods(void);

    // get realm description
    CSmallString GetDescription(void);

    /// return number of cvs
    int GetNumOfCVs(void) const;

    /// return number of bins
    int GetNumOfBins(void) const;

    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

// protected data --------------------------------------------------------------
protected:
    CPMFAccumulatorPtr                      Accu;
    std::set<CSmallString>                  Requires;
    std::map<CSmallString,CProxyRealmDescr> SupportedRealms;
    int                                     RealmID;
};

//------------------------------------------------------------------------------

#endif
