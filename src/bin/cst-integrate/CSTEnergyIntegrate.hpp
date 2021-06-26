#ifndef CCSTIntegrateH
#define CCSTIntegrateH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include "CSTEnergyIntOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <SmallTimeAndDate.hpp>
#include <PMFAccumulator.hpp>
#include <EnergyDerProxy.hpp>
#include <EnergySurface.hpp>

//------------------------------------------------------------------------------

class CIntegratorGPR;

//------------------------------------------------------------------------------

/// utility to integrate CST accumulator

class CCSTEnergyIntegrate {
public:
    CCSTEnergyIntegrate(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CCSTEnergyIntOptions    Options;
    CStdIOFile              InputFile;
    CStdIOFile              OutputFile;
    CPMFAccumulatorPtr      Accu;
    CEnergySurfacePtr       FES;
    CEnergyDerProxyPtr      DerProxy;
    CSmallTimeAndDate       StartTime;
    CSimpleVector<double>   GPos;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;
    CSmallString        CSTAccuName;
    CSmallString        FEOutputName;
    CSmallString        FullFEOutputName;
    std::vector<bool>   KeepCVs;

    void PrepareAccumulatorII(void);
    void PrintSampledStat(void);

    bool IntegrateForEcut(void);
    bool Integrate(void);
    void WriteHeader(void);

    void PrintGPRHyprms(FILE* p_fout);
    void LoadGPRHyprms(CIntegratorGPR& gpr);
    void SyncFESWithACCU(void);
    void DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname);
    bool ReduceFES(void);
};

//------------------------------------------------------------------------------

#endif
