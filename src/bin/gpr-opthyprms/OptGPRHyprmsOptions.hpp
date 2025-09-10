#ifndef OptGPRHyprmsOptionsH
#define OptGPRHyprmsOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class COptGPRHyprmsOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    COptGPRHyprmsOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "gpr-opthyprms"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program finds optimal GPR hyperparameters, which maximize logarithm of marginal likelihood. "
    "The optimization is performed by the L-BFGS optimizer employing either analytical or numerical gradients of logML w.r.t. hyperparameters."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,AccuFile)
    CSO_ARG(CSmallString,Realm)
    CSO_ARG(CSmallString,HyprmsFile)
    // options ------------------------------
    CSO_OPT(bool,ListRealms)
    CSO_OPT(CSmallString,Target)
    CSO_OPT(int,Limit)
    CSO_OPT(CSmallString,LAMethod)
    CSO_OPT(double,RCond)

    CSO_OPT(CSmallString,SigmaF2)
    CSO_OPT(double,MinSigmaF2)
    CSO_OPT(CSmallString,CoVar)
    CSO_OPT(double,MinCoVar)
    CSO_OPT(CSmallString,WFac)
    CSO_OPT(double,MinWFac)
    CSO_OPT(CSmallString,SigmaN2)
    CSO_OPT(double,MinSigmaN2)

    CSO_OPT(CSmallString,SigmaF2Enabled)
    CSO_OPT(CSmallString,CoVarEnabled)
    CSO_OPT(CSmallString,WFacEnabled)
    CSO_OPT(CSmallString,SigmaN2Enabled)

    CSO_OPT(bool,Numeric)
    CSO_OPT(int,NOptSteps)
    CSO_OPT(int,NumOfResets)
    CSO_OPT(double,TermEps)
    CSO_OPT(double,TermVal)
    CSO_OPT(int,NumOfLBFGSCorr)
    CSO_OPT(bool,Test)
    CSO_OPT(bool,PrintStat)
    CSO_OPT(bool,SPType)
    CSO_OPT(bool,CD5)
    CSO_OPT(CSmallString,LoadHyprms)
    CSO_OPT(CSmallString,GPRKernel)
    CSO_OPT(bool,UseFDKernel)
    CSO_OPT(bool,GPRUseInv)
    CSO_OPT(bool,GPRCalcLogPL)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, AccuFile, NULL, true, "ACCU",
            "Name of the file containing the input PMF accumulator.")
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, Realm, NULL, true, "REALM",
            "Realm for GPR hyperparameter optimization. The list of supported realms can be obtained by --listrealms.")
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, HyprmsFile, NULL, true, "HYPRMS",
            "Name of the file containing the optimized hyperparameters.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, ListRealms, false, false, '\0', "listrealms", NULL,
            "List supported realms.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, Target, "logml", false, 't', "target", "NAME",
            "Specify the optimized target, which can be either logml "
            "(log of marginal likelihood) or logpl "
            "(log of pseudo-likelihood from leave-one-out cross-validation, LOO-CV).")
    // -------------------------------------------
        CSO_MAP_OPT(int, Limit, 1000, false, 'l', "limit", "LIMIT",
            "Only bins containing more samples than NUMBER are considered "
            "properly sampled.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, LAMethod, "default", false, 'a', "lmethod", "NAME",
            "Linear algebra method for matrix inversion. Supported algorithms: "
            "default, svd (SVD – singular value decomposition, divide-and-conquer driver), "
            "svd2 (SVD – singular value decomposition, simple driver), "
            "lu (LU factorization), and ll (LL – Cholesky factorization).")
    // -------------------------------------------
        CSO_MAP_OPT(double, RCond, 1e-6, false, 'r', "rcond", "NUMBER",
            "Rank condition for SVD. The chosen value must be carefully tested. "
            "Calculation at computer precision is requested with -1 (not recommended).")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaF2, "15.0", false, 's', "sigmaf2", "SPEC",
            "Variance of the reconstructed free-energy surface (signal variance) "
            "in the form SigmaF2(1)[xSigmaF2(2)[x...]]. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MinSigmaF2, 0.1, false, 0, "minsigmaf2", "NUMBER",
            "Minimum value of SigmaF2.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, CoVar, "0.0", false, 'o', "covar", "SPEC",
            "Covariances between GPR tasks in the form CoVar1[xCoVar2[x...]]. "
            "The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MinCoVar, -1000, false, 0, "mincovar", "NUMBER",
            "Minimum value of CoVar.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, WFac, "3.0", false, 'w', "wfac", "SPEC",
            "Characteristic scale of collective variables in the form WFac1[xWFac2[x...]]. "
            "The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MinWFac, 0.1, false, 0, "minwfac", "NUMBER",
            "Minimum value of WFac.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaN2, "1e-5", false, 'n', "sigman2", "SPEC",
            "Noise variances (σ²) for each CV in the form SigmaN2(1)[xSigmaN2(2)[x...]]. "
            "The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MinSigmaN2, 0.0, false, 0, "minsigman2", "NUMBER",
            "Minimum value of SigmaN2.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaF2Enabled, "F", false, 0, "enablesigmaf2", "SPEC",
            "Enable optimization of SigmaF2 hyperparameters. Flags are specified "
            "in the form SigmaF2(1)Enabled[xSigmaF2(2)Enabled[x...]] with F and T "
            "for disabled and enabled, respectively. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, CoVarEnabled, "F", false, 0, "enablecovar", "SPEC",
            "Enable optimization of CoVar hyperparameters. Flags are specified "
            "in the form CoVar1Enabled[xCoVar2Enabled[x...]] with F and T "
            "for disabled and enabled, respectively. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, WFacEnabled, "F", false, 0, "enablewfac", "SPEC",
            "Enable optimization of WFac hyperparameters. Flags are specified "
            "in the form WFac1Enabled[xWFac2Enabled[x...]] with F and T "
            "for disabled and enabled, respectively. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaN2Enabled, "F", false, 0, "enablesigman2", "SPEC",
            "Enable optimization of SigmaN2 hyperparameters. Flags are specified "
            "in the form SigmaN2(1)Enabled[xSigmaN2(2)Enabled[x...]] with F and T "
            "for disabled and enabled, respectively. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Numeric, false, false, 0, "numeric", NULL,
            "Use numerical gradients.")
    // -------------------------------------------
        CSO_MAP_OPT(int, NOptSteps, 100, false, 0, "noptsteps", "NUMBER",
            "Maximum number of optimization steps.")
    // -------------------------------------------
        CSO_MAP_OPT(int, NumOfResets, 3, false, 0, "nresets", "NUMBER",
            "Maximum number of resets due to insufficient optimization progress.")
    // -------------------------------------------
        CSO_MAP_OPT(double, TermEps, 1e-5, false, 0, "termeps", "NUMBER",
            "Termination criterion for the L-BFGS optimizer (see L-BFGS code).")
    // -------------------------------------------
        CSO_MAP_OPT(double, TermVal, 1e-6, false, 0, "termval", "NUMBER",
            "Termination criterion for the L-BFGS optimizer. "
            "Minimum change of the optimized property.")
    // -------------------------------------------
        CSO_MAP_OPT(int, NumOfLBFGSCorr, 10, false, 0, "nlbfgscorr", "NUMBER",
            "Number of corrections used in an L-BFGS update.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Test, false, false, 0, "test", NULL,
            "Compare analytical and numerical gradients for input hyperparameters.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, CD5, false, false, 0, "cd5", NULL,
            "Use the 5-point stencil for numerical differentiation.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, LoadHyprms, NULL, false, 0, "loadhyprms", "NAME",
            "Name of the file containing the GPR hyperparameters.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, PrintStat, false, false, 0, "stat", NULL,
            "Calculate detailed GPR statistics using the optimized hyperparameters.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, SPType, false, false, 0, "sptype", NULL,
            "Determine the type of stationary point.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, GPRKernel, "default", false, '\0', "kernel", "NAME",
            "GPR kernel type. Supported types: ardse (ARD squared exponential), "
            "ardmc52 (ARD Matérn class 5/2), default (=ardse).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, UseFDKernel, false, false, 0, "fdkernel", NULL,
            "Use the first derivatives of the GPR kernel (only for SmootherGPR).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRUseInv, false, false, 0, "useinv", NULL,
            "Use the matrix inversion pathway (for testing only).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRCalcLogPL, false, false, 0, "calclogpl", NULL,
            "Calculate logPL for --stat if --target is not logpl.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, GlobalMin, NULL, false, '\0', "globalmin", "SPEC",
            "Position of the global minimum provided as a single string in the form "
            "CV1xCV2x...xCVn (relevant for error determination). "
            "If not set, the position is determined automatically.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Verbose, false, false, 'v', "verbose", NULL,
            "Increase output verbosity.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Version, false, false, '\0', "version", NULL,
            "Output version information and exit.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Help, false, false, 'h', "help", NULL,
            "Display this help and exit.")
    // -------------------------------------------
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
