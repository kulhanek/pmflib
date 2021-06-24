#ifndef ABFAccumulatorH
#define ABFAccumulatorH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//                       Martin Petrek, petrek@chemi.muni.cz
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

//------------------------------------------------------------------------------

enum EABFAccuValue {
    EABF_DG_VALUE       = 0,  // dG/dksi
    EABF_DG_SIGMA       = 1,  // fluctuation of ICF
    EABF_DG_ERROR       = 2,  // error of dG/dksi
    EABF_H_VALUE        = 3,  // enthalpy
    EABF_H_SIGMA        = 4,  // fluctuation of potential energy
    EABF_H_ERROR        = 5,  // error of enthalpy
    EABF_TDS_VALUE      = 6,  // -TdS/dksi
    EABF_TDS_SIGMA      = 7,  // fluctuation of -TdS/dksi
    EABF_TDS_ERROR      = 8,  // error of -TdS/dksi
};

//------------------------------------------------------------------------------

/** \brief ABF accumulator
*/

class PMF_PACKAGE CABFAccumulator : public CPMFAccumulator {
public:
// constructor and destructor --------------------------------------------------
    CABFAccumulator(void);
    ~CABFAccumulator(void);

//------------------------------------------------------------------------------
    /// return number of samples for a given bin position
    int GetNumOfSamples(const CSimpleVector<int>& position) const;

    /// return number of samples for a given bin index
    int GetNumOfSamples(int ibin) const;

    /// get value - dG/dksi
    double GetValue(int icv,int ibin,EABFAccuValue realm) const;

    /// get value - H=<PotEne>
    double GetValue(int ibin,EABFAccuValue realm) const;

//------------------------------------------------------------------------------
    /// set number of samples
    void SetNumOfSamples(int ibin,int samples);

    /// set mask weight for a given bin index
    void SetMaskWeight(int ibin,double weight);

    /// set number of statistically correlated samples (it influences calculated errors of mean forces)
    void SetNCorr(double ncorr);

// section of private data ----------------------------------------------------
private:
    double                      NCorr;          // number of correlated samples (for error evaluation)
};

typedef std::shared_ptr<CABFAccumulator>    CABFAccumulatorPtr;

//------------------------------------------------------------------------------

#endif
