#ifndef IntegratorRFDH
#define IntegratorRFDH
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
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <EnergyDerProxy.hpp>
#include <EnergySurface.hpp>

extern "C" {
#include <cs.h>
}

//------------------------------------------------------------------------------

/** \brief integrator of PMF accumulator based on finite differences
    only PMF points with number of samples higher than zero are considered
*/

class PMF_PACKAGE CIntegratorRFD {
public:
// constructor and destructor -------------------------------------------------
    CIntegratorRFD(void);
    virtual ~CIntegratorRFD(void);

// setup methods --------------------------------------------------------------
    /// set input energy derivative proxy
    void SetInputEnergyDerProxy(CEnergyDerProxyPtr p_proxy);

    /// set output free energy surface
    void SetOutputES(CEnergySurfacePtr p_surf);

    /// set FD number of points (3 or 4)
    void SetFDPoints(int npts);

    /// should we apply periodicity?
    void SetPeriodicity(bool set);

    /// use old fit
    void SetUseOldRFDMode(bool set);

    /// set position of global minimum - spec in real units
    void SetGlobalMin(const CSmallString& spec);

    /// set position of global minimum - in internal units
    void SetGlobalMin(const CSimpleVector<double>& pos);

    /// get position of global minima - in internal units
    CSimpleVector<double> GetGlobalMin(void);

    /// get position of global minima - bin
    int GetGlobalMinBin(void);

// execution method -----------------------------------------------------------
    /// integrate data, for errors the FES must be already allocated!!!
    bool Integrate(CVerboseStr& vout);

    /// get root mean square residuals
    double GetRMSR(int k);

// section of private data ----------------------------------------------------
private:
    CPMFAccumulatorPtr      Accu;
    CEnergyDerProxyPtr      DerProxy;
    CEnergySurfacePtr       EneSurf;

    int                     FDLevel;
    bool                    Periodicity;
    bool                    UseOldRFDMode;

    int                     NumOfVariables;
    int                     NumOfEquations;
    int                     NumOfNonZeros;

    bool                    GlobalMinSet;   // true if gpos set by SetGlobalMin()
    CSimpleVector<double>   GPos;           // global position, either detected or use
    bool                    GPosSet;        // true is gpos set by any means, either SetGlobalMin() or from FES
    int                     GPosBin;

    CSimpleVector<int>      XMap;       // translation between global index and X index
    CSimpleVector<int>      IPoint;     // index points
    CSimpleVector<double>   Rhs;        // right hand side
    CSimpleVector<int>      RhsCv;      // mapping to CV for RMSR calculation
    CSimpleVector<double>   X;          // unknow free energy
    cs*                     A;          // main matrix with differentation schemes
    cs*                     cA;         // compressed A
    int                     LocIter;

    /// release all resources
    void ReleaseAllResources(void);

    /// build system of equations describing differentation scheme
    bool BuildSystemOfEquations(CVerboseStr& vout);

    /// solve redundant system of linear equations
    bool SolveSystemOfEquations(void);

    /// recursively build system of equations
    void BuildEquations(bool trial);

    /// get global index of point
    int GetFBinIndex(const CSimpleVector<int>& position,int ifcoord,int offset) const;

    /// get integrated value
    double GetIntegratedValue(int icoord,int ibin);
};

//------------------------------------------------------------------------------

#endif
