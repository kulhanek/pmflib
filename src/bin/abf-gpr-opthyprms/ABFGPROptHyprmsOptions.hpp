#ifndef ABFGPROptHyprmsOptionsH
#define ABFGPROptHyprmsOptionsH
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

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class CABFGPROptHyprmsOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABFGPROptHyprmsOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-gpr-opthyprms"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program calculates optimal hyperparameters for ABF GPR algorithm."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,InAccuName)
    CSO_ARG(CSmallString,HyprmsName)
    // options ------------------------------
    CSO_OPT(CSmallString,OptTarget)
    CSO_OPT(CSmallString,OptRealm)
    CSO_OPT(int,GPRLen)
    CSO_OPT(double,Width)
    CSO_OPT(double,Noise)
    CSO_OPT(bool,WidthEnabled)
    CSO_OPT(bool,NoiseEnabled)
    CSO_OPT(bool,Numeric)
    CSO_OPT(int,NOptSteps)
    CSO_OPT(int,NumOfResets)
    CSO_OPT(double,TermEps)
    CSO_OPT(double,TermVal)
    CSO_OPT(int,NumOfLBFGSCorr)
    CSO_OPT(bool,Test)
    CSO_OPT(bool,SPType)
    CSO_OPT(bool,CD5)
    CSO_OPT(CSmallString,GPRKernel)
    CSO_OPT(double,RCond)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                InAccuName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "inaccu",                        /* parameter name */
                "Name of input PMF accumulator. If the name is '-' then the PFM accumulator is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                HyprmsName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "outaccu",                        /* parameter name */
                "Name for file with optimized hyperparameters. If the name is '-' then the file will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OptTarget,                        /* option name */
                "logml",                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "target",                      /* long option name */
                "NAME",                           /* parameter name */
                "Specify optimized targed, which can be either logml (log of marginal likelihood) "
                "or logpl (log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OptRealm,                        /* option name */
                "CVS",                          /* default value */
                false,                          /* is option mandatory */
                'r',                           /* short option name */
                "realm",                      /* long option name */
                "NAME",                           /* parameter name */
                "Specify optimized realm: CVS, ICF, KIN, TOT.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                GPRLen,                        /* option name */
                1000,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "gprlen",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "The length of series for GPR.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                RCond,                        /* option name */
                1e-6,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "rcond",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Rank condition for SVD. Used value must be carefully tested. Calculation at computer precision is requested with -1 (not recommended).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Width,                        /* option name */
                15.0,                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "width",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Kernel width.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Noise,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                'n',                           /* short option name */
                "noise",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Noise.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                WidthEnabled,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablewidth",                      /* long option name */
                NULL,                           /* parameter name */
                "Enable optimization of Width hyperparameter.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoiseEnabled,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "enablenoise",                      /* long option name */
                NULL,                           /* parameter name */
                "Enable optimization of Noise hyperparameter.")   /* option description */
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
                1e-6,                          /* default value */
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
