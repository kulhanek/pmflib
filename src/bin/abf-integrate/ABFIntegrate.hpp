#ifndef CABFIntegrateH
#define CABFIntegrateH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include "ABFIntOptions.hpp"
#include <ABFAccumulator.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <SmallTimeAndDate.hpp>
#include <ABFIntegratorGPR.hpp>

//------------------------------------------------------------------------------

/// utility to integrate ABF accumulator

class CABFIntegrate {
public:
    CABFIntegrate(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CABFIntOptions      Options;
    CStdIOFile          InputFile;
    CStdIOFile          OutputFile;
    CABFAccumulator     Accumulator;
    CEnergySurface      FES;
    CSmallTimeAndDate   StartTime;
    CSimpleVector<int>  FFSeeds;
    CSimpleVector<int>  IPos;
    CSimpleVector<int>  TPos;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;
    CSmallString        ABFAccuName;
    CSmallString        FEOutputName;
    CSmallString        FullFEOutputName;

    void PrepareAccumulatorI(void);
    void PrepareAccumulatorII(void);
    void PrintSampledStat(void);
    void FloodFillTest(void);
    bool InstallNewSeed(int seedid,bool unsampled);
    int  FillSeed(int seedid,bool unsampled);
    void GetTPoint(CSimpleVector<int>& ipos,int d,CSimpleVector<int>& tpos);
    bool IntegrateForMFZScore(int pass);
    bool IntegrateForEcut(void);
    bool Integrate(void);
    void WriteHeader(void);
    void GlueingFES(int factor);
    void GlueHoles(void);
    void SeedSampled(int seedid);
    bool IsHole(int seedid);
    void MarkAsHole(int seedid);
    void PrintGPRHyprms(FILE* p_fout);
    void LoadGPRHyprms(CABFIntegratorGPR& gpr);
};

//------------------------------------------------------------------------------

#endif
