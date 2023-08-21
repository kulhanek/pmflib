#ifndef OptGPRHyprmsOptionsH
#define OptGPRHyprmsOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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
    CSO_OPT(CSmallString,NCorr)
    CSO_OPT(double,MinNCorr)
    CSO_OPT(CSmallString,SigmaN2)
    CSO_OPT(double,MinSigmaN2)

    CSO_OPT(CSmallString,SigmaF2Enabled)
    CSO_OPT(CSmallString,CoVarEnabled)
    CSO_OPT(CSmallString,WFacEnabled)
    CSO_OPT(CSmallString,NCorrEnabled)
    CSO_OPT(CSmallString,SigmaN2Enabled)

    CSO_OPT(bool,EnableConstraints)

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
    CSO_OPT(bool,GPRUseInv)
    CSO_OPT(bool,GPRCalcLogPL)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                AccuFile,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "ACCU",                        /* parameter name */
                "Name of file containing the input ABF accumulator.\n")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Realm,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "REALM",                        /* parameter name */
                "Realm for GPR hyperparameter optimization. Supported realms: dG/dx, dH/dx, dH, and mTDS/dx (-TdS/dx).\n")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                HyprmsFile,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "HYPRMS",                        /* parameter name */
                "Name of file containing the optimized hyperparameters.\n")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Target,                        /* option name */
                "logml",                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "target",                      /* long option name */
                "NAME",                           /* parameter name */
                "Specify optimized targed, which can be either logml (log of marginal likelihood) "
                "or logpl (log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                1000,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "LIMIT",                           /* parameter name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                LAMethod,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                'a',                           /* short option name */
                "lmethod",                      /* long option name */
                "NAME",                           /* parameter name */
                "Linear algebra method for matrix inversion. Supported algorithms are: "
                "default, svd (SVD - singular value decomposition, divide and conquer driver), "
                "svd2 (SVD - singular value decomposition, simple driver), lu (LU factorization), "
                "and ll (LL - Cholesky factorization).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                RCond,                        /* option name */
                1e-6,                          /* default value */
                false,                          /* is option mandatory */
                'r',                           /* short option name */
                "rcond",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Rank condition for SVD. Used value must be carefully tested. Calculation at computer precision is requested with -1 (not recommended).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SigmaF2,                        /* option name */
                "15.0",                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "sigmaf2",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Variance of the reconstructed free energy surface (signal variance) in the form SigmaF2(1)[xSigmaF2(2)[x...]]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MinSigmaF2,                        /* option name */
                0.1,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "minsigmaf2",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Minimal value of SigmaF2.")   /* option description */
//----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                CoVar,                        /* option name */
                "0.0",                          /* default value */
                false,                          /* is option mandatory */
                'o',                           /* short option name */
                "covar",                      /* long option name */
                "SPEC",                           /* parameter name */
                "CoVariances between GPR tasks in the form CoVar1[xCoVar2[x...]]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MinCoVar,                        /* option name */
                -1000,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "mincovar",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Minimal value of CoVariances.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                WFac,                        /* option name */
                "3.0",                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "wfac",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Characteristic scale of collective variables in the form WFac1[xWFac2[x...]]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MinWFac,                        /* option name */
                0.1,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "minwfac",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Minimal value of WFac.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                NCorr,                        /* option name */
                "0.0",                          /* default value */
                false,                          /* is option mandatory */
                'c',                           /* short option name */
                "ncorr",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Number of statistically correlated samples in the form NCorr1[xNCorr2[x...]]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MinNCorr,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "minncorr",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Minimal value of NCorr.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SigmaN2,                        /* option name */
                "0.0",                          /* default value */
                false,                          /* is option mandatory */
                'n',                           /* short option name */
                "sigman2",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Values of noise sigma squared for each CV in the form SigmaN2(1)[xSigmaN2(2)[x...]]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MinSigmaN2,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "minsigman2",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Minimal value of SigmaN2.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SigmaF2Enabled,                        /* option name */
                "F",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablesigmaf2",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Enable optimization of SigmaF2 hyperparamters. Flags are specified in the form SigmaF2(1)Enabled[xSigmaF2(1)2Enabled[x...]] with F and T for disabled and enabled, respectively. "
                "The last value pads the rest.")   /* option description */
//----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                CoVarEnabled,                        /* option name */
                "F",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablecovar",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Enable optimization of CoVar hyperparamters. Flags are specified in the form CoVar1Enabled[xCoVar2Enabled[x...]] with F and T for disabled and enabled, respectively. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                WFacEnabled,                        /* option name */
                "F",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablewfac",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Enable optimization of Wfac hyperparamters. Flags are specified in the form WFac1Enabled[xWFac2Enabled[x...]] with F and T for disabled and enabled, respectively. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                NCorrEnabled,                        /* option name */
                "F",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablencorr",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Enable optimization of NCorr hyperparamters. Flags are specified in the form NCorr1Enabled[xNCorr2Enabled[x...]] with F and T for disabled and enabled, respectively. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SigmaN2Enabled,                        /* option name */
                "F",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablesigman2",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Enable optimization of SigmaN2 hyperparamters. Flags are specified in the form SigmaN2(1)Enabled[xSigmaN2(1)Enabledx...] with F and T for disabled and enabled, respectively. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                EnableConstraints,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "constraints",                      /* long option name */
                NULL,                           /* parameter name */
                "Enable constraints in multitask GPR, e.g., impose dG(x)/dx - dH(x)/dx - (-TdS(x)/dx) = 0 constraints.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Numeric,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "numeric",                      /* long option name */
                NULL,                           /* parameter name */
                "Use numeric gradients.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                NOptSteps,                        /* option name */
                100,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "noptsteps",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Maximum number of optimization steps.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                NumOfResets,                        /* option name */
                3,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "nresets",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Maximum number of resets due to insufficient optimization progress.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                TermEps,                        /* option name */
                1e-5,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "termeps",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Termination criteria for L-BFGS optimizer (see L-BFGS code).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                TermVal,                        /* option name */
                1e-7,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "termval",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Termination criteria for L-BFGS optimizer. Minimum change of optimized property.")   /* option description */
//----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                NumOfLBFGSCorr,                        /* option name */
                10,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "nlbfgscorr",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Number of corrections used in an L-BFGS update.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Test,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "test",                      /* long option name */
                NULL,                           /* parameter name */
                "Compare analytical and numerical gradient for input hyperparameters.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                CD5,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "cd5",                      /* long option name */
                NULL,                           /* parameter name */
                "Use 5-point stencil for numerical differentiation.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                LoadHyprms,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "loadhyprms",                      /* long option name */
                "NAME",                           /* parameter name */
                "Name of file with GPR hyperparameters.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                PrintStat,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "stat",                      /* long option name */
                NULL,                           /* parameter name */
                "Calculate detailed GPR status employing found hyperparameters.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                SPType,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "sptype",                      /* long option name */
                NULL,                           /* parameter name */
                "Determine type of stationary point.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                GPRKernel,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "kernel",                      /* long option name */
                "NAME",                           /* parameter name */
                "GPR: Kernel type. Supported types: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), default(=ardse)")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRUseInv,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "useinv",                      /* long option name */
                NULL,                           /* parameter name */
                "Use matrix inversion pathway (for testing only).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRCalcLogPL,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "calclogpl",                      /* long option name */
                NULL,                           /* parameter name */
                "Calculate logPL for --stat if --target is not logpl.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                GlobalMin,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "globalmin",                      /* long option name */
                "SPEC",                           /* parameter name */
                "position of global minimum provided as a single string in the form CV1xCV2x...xCVn (relevant for error determination), if not set the position is determined automatically.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Verbose,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'v',                           /* short option name */
                "verbose",                      /* long option name */
                NULL,                           /* parameter name */
                "Increase output verbosity.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Version,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "version",                      /* long option name */
                NULL,                           /* parameter name */
                "Output version information and exit.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Help,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'h',                           /* short option name */
                "help",                      /* long option name */
                NULL,                           /* parameter name */
                "Display this help and exit.")   /* option description */
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
