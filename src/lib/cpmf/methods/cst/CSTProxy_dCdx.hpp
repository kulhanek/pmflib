#ifndef CSTProxy_dCdx_H
#define CSTProxy_dCdx_H
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

enum ECSTdCdxType {

    // TDS subsystem
    CST_TdS_LT,         // Cov(lambda,Etot)
    CST_TdS_LI,         // Cov(lambda,Eint)
    CST_TdS_LP,         // Cov(lambda,Epot)
    CST_TdS_LR,         // Cov(lambda,Erst)
    CST_TdS_LK,         // Cov(lambda,Ekin)

    CST_TdS_LTFW,       // Cov(lambda,Etot) - Fixman weighted

    // TDS+VF subsystem
    CST_TdS_II,         // Cov(ICF,Eint)
    CST_TdS_IIFW,       // Cov(ICF,Eint)    - Fixman weighted
    CST_TdS_PIFW,       // Cov(ICFP,Eint)   - Fixman weighted
    CST_TdS_KIFW,       // Cov(ICFK,Eint)   - Fixman weighted
};

//------------------------------------------------------------------------------

/** \brief ABF proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CCSTProxy_dCdx : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CCSTProxy_dCdx(void);
    ~CCSTProxy_dCdx(void);

//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy derivative and its error
    virtual double GetValue( int ibin,int icv,EProxyRealm realm) const;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CCSTProxy_dCdx>    CCSTProxy_dCdx_Ptr;

//------------------------------------------------------------------------------

#endif
