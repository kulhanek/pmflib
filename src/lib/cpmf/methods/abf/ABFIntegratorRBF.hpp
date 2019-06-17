#ifndef ABFIntegratorRBFH
#define ABFIntegratorRBFH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>

//------------------------------------------------------------------------------

class CABFAccumulator;
class CEnergySurface;

// integration method

enum EARBFMethod {
    EARBF_SVD = 1,
    EARBF_QRLQ = 2,
};

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing radial basis functions
*/

class PMF_PACKAGE CABFIntegratorRBF {
public:
// constructor and destructor -------------------------------------------------
    CABFIntegratorRBF(void);
    virtual ~CABFIntegratorRBF(void);

// setup methods --------------------------------------------------------------
    /// set input ABF accumulator, only ABF forces are integrated
    void SetInputABFAccumulator(const CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

    /// set verbosity level
    void SetVerbosity(bool set);

    /// multiply of bin sizes
    void SetGaussianWidth(int order);

    /// set limit for number of samples
    void SetLimit(int limit);

    /// should we apply periodicity?
    void SetPeriodicity(bool set);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout);

// section of private data ----------------------------------------------------
private:
    const CABFAccumulator*  Accumulator;
    CEnergySurface*         FES;

    bool                    Verbose;
    int                     WidthOrder;
    bool                    Periodicity;
    int                     Limit;
    EARBFMethod             Method;

    // RBF data
    int                     NumOfRBFs;
    CSimpleVector<double>   Weights;
    CSimpleVector<int>      NumOfRBFBins;
    int                     NumOfCVs;
    CSimpleVector<double>   Sigmas;

    // SVD setup
    double                  RCond;

    bool IntegrateByLS(CVerboseStr& vout);

    void    GetRBFPosition(unsigned int index,CSimpleVector<double>& position);
    double  GetValue(const CSimpleVector<double>& position);
};

//------------------------------------------------------------------------------

#endif
