#ifndef ABFProxy_mTdS_H
#define ABFProxy_mTdS_H
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

enum EABFTdSType {
    ABF_TdS_HH,

    ABF_TdS_HP,
    ABF_TdS_HR,
    ABF_TdS_HK,
    ABF_TdS_HV,

    ABF_TdS_BP,
    ABF_TdS_BR,
    ABF_TdS_BK,
    ABF_TdS_BV,
};

//------------------------------------------------------------------------------

/** \brief ABF proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CABFProxy_mTdS : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CABFProxy_mTdS(void);
    ~CABFProxy_mTdS(void);

//------------------------------------------------------------------------------
    // set type
    void SetType(EABFTdSType type);

//------------------------------------------------------------------------------
    // get energy derivative and its error
    virtual double GetValue( int ibin,int icv,EProxyRealm realm) const;

    // is compatible with PMF Accumulator method
    static bool IsCompatible(CPMFAccumulatorPtr accu);

// section of private data -----------------------------------------------------
private:
    EABFTdSType    Type;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CABFProxy_mTdS>    CABFProxy_mTdS_Ptr;

//------------------------------------------------------------------------------

#endif
