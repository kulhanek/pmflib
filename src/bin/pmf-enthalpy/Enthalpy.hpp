#ifndef CEnthalpyH
#define CEnthalpyH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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

#include "EnthalpyOptions.hpp"
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

class CEnthalpy {
public:
    CEnthalpy(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CEnthalpyOptions                Options;
    CSmallString                    HEOutputName;
    CStdIOFile                      OutputFile;
    std::vector<CEnergyProxyPtr>    EnergyProxies;
    std::vector<CPMFAccumulatorPtr> Accumulators;
    CEnergySurfacePtr               HES;
    int                             State;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// helper methods
    void GetRawEnthalpy(void);
    void LoadGPRHyprms(CSmootherGPR& gpr);
    bool PrintHES(void);
    void WriteHeader(void);
    void PrintSampledStat(void);
    void AdjustGlobalMin(void);
};

//------------------------------------------------------------------------------

#endif
