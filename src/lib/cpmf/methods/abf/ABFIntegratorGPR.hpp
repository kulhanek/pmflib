#ifndef ABFIntegratorGPRH
#define ABFIntegratorGPRH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <FortranMatrix.hpp>

//------------------------------------------------------------------------------

class CABFAccumulator;
class CEnergySurface;

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing radial basis functions
*/

class PMF_PACKAGE CABFIntegratorGPR {
public:
// constructor and destructor -------------------------------------------------
    CABFIntegratorGPR(void);
    virtual ~CABFIntegratorGPR(void);

// setup methods --------------------------------------------------------------
    /// set input ABF accumulator, only ABF forces are integrated
    void SetInputABFAccumulator(const CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

    /// multiply of bin sizes
    void SetWFac(double wfac);

    /// set sigmaf2
    void SetSigmaF2(double sigf2);

    /// set include error
    void SetIncludeError(bool set);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout);

    /// get root mean square residuals
    double GetRMSR(void);

// section of private data ----------------------------------------------------
private:
    const CABFAccumulator*  Accumulator;
    CEnergySurface*         FES;

    // GPR data
    int                     GPRSize;
    int                     NCVs;
    double                  WFac;
    CSimpleVector<double>   CVLengths2;
    double                  SigmaF2;
    CFortranMatrix          K;          // kernels
    CSimpleVector<double>   GPRModel;   // weights
    bool                    IncludeError;

    CSimpleVector<double>   ipos;
    CSimpleVector<double>   jpos;
    CSimpleVector<double>   gpos;
    CSimpleVector<double>   rk;
    CSimpleVector<double>   lk;
    CSimpleVector<double>   ik;

    bool TrainGP(CVerboseStr& vout);
    void CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout); // gpos - position of global minimum

    double GetValue(const CSimpleVector<double>& position);
    double GetMeanForce(const CSimpleVector<double>& position,int icoord);
    double GetCov(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos);
};

//------------------------------------------------------------------------------

#endif
