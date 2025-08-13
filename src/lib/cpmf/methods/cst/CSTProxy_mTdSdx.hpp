#ifndef CSTProxy_mTdSdx_H
#define CSTProxy_mTdSdx_H
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

enum ECSTTdSdxType {
    CST_mTdSdx,

    CST_TdS_LT,         // Cov(lambda,Etot)

    CST_TdS_LTFW,       // Cov(lambda,Etot) - Fixman weighted
    CST_TdS_LIFW,       // Cov(lambda,Eint) - Fixman weighted
    CST_TdS_LPFW,       // Cov(lambda,Epot) - Fixman weighted
    CST_TdS_LRFW,       // Cov(lambda,Erst) - Fixman weighted
    CST_TdS_LKFW,       // Cov(lambda,Ekin) - Fixman weighted

    CST_TdS_II,         // Cov(ICF,Eint)
    CST_TdS_IIFW,       // Cov(ICF,Eint) - Fixman weighted
};

//------------------------------------------------------------------------------

/** \brief CST proxy providing mean force for the free energy integration
*/

class PMF_PACKAGE CCSTProxy_mTdSdx : public CEnergyDerProxy {
public:
// constructor and destructor --------------------------------------------------
    CCSTProxy_mTdSdx(void);
    ~CCSTProxy_mTdSdx(void);

//------------------------------------------------------------------------------
    // get optional energy correction - MTC
    virtual CEnergyProxyPtr GetEnergyCorrection(void);

//------------------------------------------------------------------------------
    // get number of samples
    virtual int GetNumOfSamples(int ibin) const;

    // set number of samples
    virtual void SetNumOfSamples(int ibin,int nsamples);

    // get energy derivative and its error
    virtual double GetValue( int ibin,int icv,EProxyRealm realm) const;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr   <CCSTProxy_mTdSdx>    CCSTProxy_mTdSdx_Ptr;

//------------------------------------------------------------------------------

#endif
