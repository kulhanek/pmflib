#ifndef CSTProxy_mTdS_H
#define CSTProxy_mTdS_H
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
#include <EnergyDerProxy.hpp>

//------------------------------------------------------------------------------

enum ECSTTdSType {
    CST_C11HH,
    CST_C11HP,
    CST_C11HK,
    CST_C11HR,
};

//------------------------------------------------------------------------------

/** \brief CST proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CCSTProxy_mTdS : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CCSTProxy_mTdS(void);
    ~CCSTProxy_mTdS(void);

//------------------------------------------------------------------------------
    // set type
    void SetType(ECSTTdSType type);

//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy derivative and its error
    virtual double GetValue( int ibin,int icv,EProxyRealm realm) const;

    // is compatible with PMFAccumulator method
    static bool IsCompatible(CPMFAccumulatorPtr accu);

// section of private data -----------------------------------------------------
private:
    ECSTTdSType    Type;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CCSTProxy_mTdS>    CCSTProxy_mTdS_Ptr;

//------------------------------------------------------------------------------

#endif
