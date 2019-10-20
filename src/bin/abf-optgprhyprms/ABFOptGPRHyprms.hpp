#ifndef CABFOptGPRHyprmsH
#define CABFOptGPRHyprmsH
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

#include "ABFOptGPRHyprmsOptions.hpp"
#include <ABFAccumulator.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <SmallTimeAndDate.hpp>
#include <ABFIntegratorGPR.hpp>
#include <vector>

//------------------------------------------------------------------------------

/// utility to integrate ABF accumulator

class CABFOptGPRHyprms {
public:
    CABFOptGPRHyprms(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CABFOptGPRHyprmsOptions Options;
    CStdIOFile              InputFile;
    CStdIOFile              OutputFile;
    CABFAccumulator         Accumulator;
    CEnergySurface          FES;
    CSmallTimeAndDate       StartTime;

    double                  SigmaF2;
    CSimpleVector<double>   NCorr;
    CSimpleVector<double>   WFac;

// L-BFGS setup
    int                     NumOfCorrections;
    bool                    SigmaF2Enabled;
    bool                    SplitNCorr;
    std::vector<bool>       NCorrEnabled;
    std::vector<bool>       WFacEnabled;
    int                     NumOfPrms;
    int                     NumOfOptPrms;
    CSimpleVector<double>   Hyprms;
    CSimpleVector<double>   HyprmsGrd;
    CSimpleVector<double>   Work;
    CSimpleVector<double>   TmpXG;
    double                  logML;

// helper
    std::vector<bool>       HyprmsEnabled;
    int                     State;

    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    void PrintSampledStat(void);
    void InitOptimizer(void);
    void Test(void);
    void Optimize(void);
    void RunGPRAnalytical(void);
    void RunGPRNumerical(void);
    void WriteResults(int istep);
    bool WriteHyperPrms(FILE* p_fout);

    void    ScatterHyprms(CSimpleVector<double>& hyprsm);
    double  GetLogML(CABFIntegratorGPR& gpr);
    void    DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname);
    void    DecodeVList(const CSmallString& spec, CSimpleVector<double>& vlist,const CSmallString& optionname,double defv);

};

//------------------------------------------------------------------------------

#endif
