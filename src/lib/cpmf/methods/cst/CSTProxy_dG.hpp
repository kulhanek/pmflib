#ifndef CSTProxy_dG_H
#define CSTProxy_dG_H
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

/** \brief CST proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CCSTProxy_dG : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CCSTProxy_dG(void);
    ~CCSTProxy_dG(void);

//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy derivative and its error
    virtual double GetValue( int ibin,int icv,EProxyRealm realm) const;

    // is compatible with PMFAccumulator method
    static bool IsCompatible(CPMFAccumulatorPtr accu);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CCSTProxy_dG>    CCSTProxy_dG_Ptr;

//------------------------------------------------------------------------------

#endif
