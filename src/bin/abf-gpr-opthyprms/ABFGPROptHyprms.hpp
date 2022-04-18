#ifndef ABFGPROptHyprmsH
#define ABFGPROptHyprmsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include "ABFGPROptHyprmsOptions.hpp"
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <PMFAccumulator.hpp>
#include <FortranMatrix.hpp>
#include <GPFilter.hpp>

//------------------------------------------------------------------------------

enum EGPROptTarget {
    EGOT_LOGML  = 1,    // log of marginal likelihood
    EGOT_LOGPL  = 2,     // log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
};

//------------------------------------------------------------------------------

enum EGPROptRealm {
    EGOR_CVS    = 1,    // CVS
    EGOR_ICF    = 2,    // ICF
    EGOR_KIN    = 3,    // kinetic energy
    EGOR_TOT    = 4,    // total energy
};

//------------------------------------------------------------------------------

/// utility to optimize ABF GPR hyperparameters

class CABFGPROptHyprms {
public:
    CABFGPROptHyprms(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CABFGPROptHyprmsOptions     Options;
    CStdIOFile                  InputFile;
    CStdIOFile                  OutputFile;
    CPMFAccumulatorPtr          InAccu;
    int                         State;

    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    size_t                  NSTLimit;
    size_t                  NumOfBins;
    size_t                  NumOfCVs;
    double                  TimeStep;

// other GPR setup
    int                     GPRLen;
    CVectorDataPtr          InData;
    int                     InCVS;

// hyperparameters
    double                  Width;
    double                  Noise;

// L-BFGS setup
    int                     NumOfCorrections;
    bool                    WidthEnabled;
    bool                    NoiseEnabled;
    int                     NumOfPrms;
    int                     NumOfOptPrms;
    CSimpleVector<double>   Hyprms;
    CSimpleVector<double>   HyprmsGrd;
    CSimpleVector<double>   Work;
    CSimpleVector<double>   TmpXG;
    EGPROptRealm            OptRealm;
    EGPROptTarget           OptTarget;
    double                  logTarget;
    int                     KRank;

// helper
    std::vector<bool>       HyprmsEnabled;

// hessian
    CFortranMatrix          Hessian;
    CSimpleVector<double>   EigenValues;


    void InitOptimizer(void);
    void Test(void);
    bool Optimize(void);
    void RunGPRAnalytical(void);
    void RunGPRNumerical(void);
    void WriteResults(int istep);
    bool WriteHyperPrms(FILE* p_fout);
    void LoadGPRHyprms(void);
    void PrintGradientSummary(void);
    void ShowGPRStat(void);
    void CalcHessian(void);
    double GetGNorm(void);
    double GetTargetFromGPFilter(CGPFilter& gpr,std::vector<bool>& hypen,CSimpleVector<double>& grd,int& krank);

    void ScatterHyprms(CSimpleVector<double>& hyprsm);
    bool ResetOpt(int& numofreset);

    void GetEtot(void);
};

//------------------------------------------------------------------------------

#endif
