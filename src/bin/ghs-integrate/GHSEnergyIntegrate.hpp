#ifndef CABFIntegrateH
#define CABFIntegrateH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Martin Petrek, petrek@chemi.muni.cz
//                       Petr Kulhanek, kulhanek@enzim.hu
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

#include "GHSEnergyIntOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <SmallTimeAndDate.hpp>
#include <PMFAccumulator.hpp>
#include <EnergyDerProxy.hpp>
#include <EnergySurface.hpp>

//------------------------------------------------------------------------------

class CIntegratorGPR;
class CSmootherGPR;

//------------------------------------------------------------------------------

/// utility to integrate ABF accumulator

class CGHSEnergyIntegrate {
public:
    CGHSEnergyIntegrate(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CGHSEnergyIntOptions    Options;

    // input
    CPMFAccumulatorPtr      Accu;
    CEnergyDerProxyPtr      GDerProxy;
    CEnergyDerProxyPtr      HDerProxy;
    CEnergyProxyPtr         HEneProxy;
    CEnergyDerProxyPtr      SDerProxy;

    // output
    CEnergySurfacePtr       FES;
    CEnergySurfacePtr       HES;
    CEnergySurfacePtr       SES;

    // sampled data
    CSimpleVector<int>      FFSeeds;
    CSimpleVector<int>      IPos;
    CSimpleVector<int>      TPos;

    // helpers
    CTerminalStr            Console;
    CVerboseStr             vout;
    CSmallTimeAndDate       StartTime;
    int                     State;

    void PrintAccuStat(void);
    bool Integrate0C(void);
    void WriteES(CEnergySurfacePtr& surf,const CSmallString& name);
};

//------------------------------------------------------------------------------

#endif
