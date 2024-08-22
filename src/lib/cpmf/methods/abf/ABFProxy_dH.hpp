#ifndef ABFProxy_dH_H
#define ABFProxy_dH_H
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
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

enum EABFdHType {
    ABF_dH,
    ABF_MICFP,
};

//------------------------------------------------------------------------------

/** \brief ABF proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CABFProxy_dH : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CABFProxy_dH(void);
    ~CABFProxy_dH(void);

//------------------------------------------------------------------------------
    // set type
    void SetType(EABFdHType type);

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
    EABFdHType    Type;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CABFProxy_dH>    CABFProxy_dH_Ptr;

//------------------------------------------------------------------------------

#endif
