#ifndef ABFIntegratorRFD2H
#define ABFIntegratorRFD2H
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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
#include <SimpleVector.hpp>
#include <FortranMatrix.hpp>
#include <VerboseStr.hpp>
#include <ABFAccumulator.hpp>

//------------------------------------------------------------------------------

// how to find solution for least square problem

enum ERFDLLSMethod {
    ERFDLLS_SVD   = 1,
    ERFDLLS_QR    = 2,    // we have overdetermined system, thus only QR is applicable
};

//------------------------------------------------------------------------------

class CABFAccumulator;
class CEnergySurface;

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator based on finite differences
    ABF points with number of samples higher than zero are considered,
    only ABF data part is integrated
*/

class PMF_PACKAGE CABFIntegratorRFD2 {
public:
// constructor and destructor -------------------------------------------------
    CABFIntegratorRFD2(void);
    virtual ~CABFIntegratorRFD2(void);

// setup methods --------------------------------------------------------------
    /// set input ABF accumulator, only ABF forces are integrated
    void SetInputABFAccumulator(const CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

    /// set FD number of points (3 or 4)
    void SetFDPoints(int npts);

    /// should we apply periodicity?
    void SetPeriodicity(bool set);

    /// set algorithm for LLS
    void SetLLSMehod(ERFDLLSMethod set);

    /// set rcond for SVD
    void SetRCond(double rcond);

// execution method -----------------------------------------------------------
    /// integrate data, for errors the FES must be already allocated!!!
    bool Integrate(CVerboseStr& vout);

    /// get root mean square residuals
    double GetRMSR(int k);

// section of private data ----------------------------------------------------
private:
    const CABFAccumulator*  Accumulator;
    CEnergySurface*         FES;

    int                     FDLevel;
    bool                    Periodicity;

    int                     NumOfVariables;
    int                     NumOfEquations;
    ERFDLLSMethod           Method;

    CSimpleVector<int>      XMap;       // translation between global index and X index
    CSimpleVector<int>      IPoint;     // index points
    CSimpleVector<double>   Rhs;        // right hand side
    CSimpleVector<int>      RhsCv;      // mapping to CV for RMSR calculation
    CSimpleVector<double>   Lhs;        // left hand side
    CSimpleVector<double>   X;          // unknow free energy
    CFortranMatrix          A;          // main matrix with differentation schemes
    CFortranMatrix          BA;         // copy of A
    CSimpleVector<double>   BRhs;       // copy of Rhs

    // SVD setup
    double                  RCond;

    /// build system of equations describing differentation scheme
    bool BuildSystemOfEquations(CVerboseStr& vout);

    /// solve redundant system of linear equations
    bool SolveSystemOfEquations(CVerboseStr& vout);

    /// recursively build system of equations
    void BuildEquations(bool trial);

    /// get global index of point
    int GetFBinIndex(const CSimpleVector<int>& position,int ifcoord,int offset) const;

    /// get integrated value
    double GetIntegratedValue(int icoord,int ibin);
};

//------------------------------------------------------------------------------

#endif
