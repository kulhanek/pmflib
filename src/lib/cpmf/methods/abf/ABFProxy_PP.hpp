#ifndef ABFProxy_PP_H
#define ABFProxy_PP_H
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

enum EABFPPNPType {
    ABF_PPNP_HH,
};

//------------------------------------------------------------------------------

/** \brief ABF proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CABFProxy_PP : public CEnergyProxy {
public:
// constructor and destructor --------------------------------------------------
    CABFProxy_PP(void);
    ~CABFProxy_PP(void);

//------------------------------------------------------------------------------
    // set type
    void SetType(EABFPPNPType type);

//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;


    // get energy derivative and its error
    virtual double GetValue( int ibin,EProxyRealm realm) const;

    // is compatible with PMF Accumulator method
    static bool IsCompatible(CPMFAccumulatorPtr accu);

// section of private data -----------------------------------------------------
private:
    EABFPPNPType    Type;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CABFProxy_PP>    CABFProxy_PP_Ptr;

//------------------------------------------------------------------------------

#endif
