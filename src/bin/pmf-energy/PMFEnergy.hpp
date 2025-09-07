#ifndef CPMFEnergyH
#define CPMFEnergyH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "PMFEneOptions.hpp"
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <PMFAccumulator.hpp>
#include <EnergySurface.hpp>
#include <EnergyProxy.hpp>
#include <SmootherGPR.hpp>

//------------------------------------------------------------------------------

/// utility to extract enthalpy from  accumulator

class CPMFEnergy {
public:
    CPMFEnergy(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CPMFEneOptions          Options;
    CStdIOFile              OutputFile;
    CPMFAccumulatorPtr      Accu;
    CEnergyProxyPtr         EneProxy;
    CEnergySurfacePtr       ENE;
    int                     State;

    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    /// helper methods
    void GetRawEnthalpy(void);
    void LoadGPRHyprms(CSmootherGPR& gpr);
    bool PrintENE(void);
    void WriteHeader(void);
    void PrintSampledStat(void);
    void AdjustGlobalMin(void);
};

//------------------------------------------------------------------------------

#endif
