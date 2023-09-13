#ifndef ABFIntOptionsH
#define ABFIntOptionsH
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

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class CGHSEnergyIntOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CGHSEnergyIntOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-integrate"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program numerically integrates data from the ABF calculation."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,AccuFile)
    CSO_ARG(CSmallString,FESFile)
    CSO_ARG(CSmallString,HESFile)
    CSO_ARG(CSmallString,SESFile)
    // options ------------------------------
    CSO_OPT(CSmallString,Realm)
    CSO_OPT(bool,EnableConstraints)
    CSO_OPT(CSmallString,LAMethod)
    CSO_OPT(double,RCond)
    CSO_OPT(int,Limit)
    CSO_OPT(CSmallString,SigmaF2)
    CSO_OPT(CSmallString,CoVar)
    CSO_OPT(CSmallString,WFac)
    CSO_OPT(CSmallString,SigmaN2)
    CSO_OPT(CSmallString,LoadHyprms)
    CSO_OPT(CSmallString,RFac)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(double,Offset)
    CSO_OPT(bool,WithError)
    CSO_OPT(bool,NoEnergy)
    CSO_OPT(bool,BalanceResiduals)
    CSO_OPT(CSmallString,OutputFormat)
    CSO_OPT(bool,UnsampledAsMaxE)
    CSO_OPT(double,MaxEnergy)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(bool,IncludeBinStat)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OEFormat)
    CSO_OPT(CSmallString,MFInfo)
    CSO_OPT(CSmallString,GPRKernel)
    CSO_OPT(bool,GPRNumDiff)
    CSO_OPT(bool,GPRUseInv)
    CSO_OPT(bool,GPRCalcLogPL)
    CSO_OPT(bool,GPRNoFastError)
    CSO_OPT(double,SLevel)
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
                FESFile,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "FES",                        /* parameter name */
                "Name of file containing the output free energy surface [dG(x)].\n")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                HESFile,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "HES",                        /* parameter name */
                "Name of file containing the output enthalpy surface [dH(x)].\n")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                SESFile,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "SES",                        /* parameter name */
                "Name of file containing the output entropic energy surface [-TdS(x)].\n")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                EnableConstraints,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "constraints",                      /* long option name */
                NULL,                           /* parameter name */
                "Impose dG(x)/dx - dH(x)/dx - (-TdS(x)/dx) = 0 constraints.")   /* option description */
// description of options ---------------------------------------------------
   CSO_MAP_OPT(CSmallString,                           /* option type */
                Realm,                        /* option name */
                "GHS_dH",                          /* default value */
                false,                          /* is option mandatory */
                'r',                           /* short option name */
                "realm",                      /* long option name */
                "NAME",                           /* parameter name */
                "Requested output from the integration:\n"
                "**  GHS_dH    - 0B\n"
                "**  GHS_dH/dx - 0A\n"
                )
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                LAMethod,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                'a',                           /* short option name */
                "lmethod",                      /* long option name */
                "NAME",                           /* parameter name */
                "Linear algebra method for LLS solution or matrix inversion. Supported algorithms are: "
                "default, svd (SVD - singular value decomposition, divide and conquer driver), "
                "svd2 (SVD - singular value decomposition, simple driver), qr (QR factorization), "
                "lu (LU factorization), ll (LL - Cholesky factorization). "
                "Possible combinations are: GPR(LU,SVD,SVD2,LL,default).")   /* option description */
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
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                1000,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SigmaF2,                        /* option name */
                "15.0",                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "sigmaf2",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Values of signal variances for each realm in the form SigmaF2(1)[xSigmaF2(2)x...]. "
                "The last value pads the rest. There are three realms.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                CoVar,                        /* option name */
                "0.0",                          /* default value */
                false,                          /* is option mandatory */
                'c',                           /* short option name */
                "covar",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Values of signal co-variances for each realm in the form CoVar(1)[CoVar(2)x...]. "
                "The last value pads the rest. There are three realms.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                WFac,                        /* option name */
                "3.0",                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "wfac",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Factors influencing widths of RBFs or square exponential kernels. The width is distance between "
                "the adjacent square exponential functions multiplied by this factors in the form WFac1[xWFac2x...]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SigmaN2,                        /* option name */
                "0.0",                          /* default value */
                false,                          /* is option mandatory */
                'n',                           /* short option name */
                "sigman2",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Values of noise sigma squared for each CV in the form SigmaN2(1)[xSigmaN2(2)x...]. "
                "The last value pads the rest.")   /* option description */
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
    CSO_MAP_OPT(CSmallString,                           /* option type */
                GlobalMin,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "globalmin",                      /* long option name */
                "SPEC",                           /* parameter name */
                "position of global minimum provided as a single string in the form CV1xCV2x...xCVn (relevant for error determination), if not set the position is determined automatically.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Offset,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                'o',                           /* short option name */
                "offset",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Specify an integration constant.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                WithError,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'e',                           /* short option name */
                "witherror",                      /* long option name */
                NULL,                           /* parameter name */
                "Estimate free energy errors.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoEnergy,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "noenergy",                      /* long option name */
                NULL,                           /* parameter name */
                "Skip calculation of energy and errors (it can save some time when only logML is required).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                BalanceResiduals,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "balres",                      /* long option name */
                NULL,                           /* parameter name */
                "Balance residual errors between dG and dH and -TdS.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                MFInfo,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "mfinfo",                      /* long option name */
                "NAME",                           /* parameter name */
                "Name of file with input and predicted mean forces.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OutputFormat,                        /* option name */
                "gnuplot",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "output",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT to print the free energy surface. Supported formats are: plain, and gnuplot.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                UnsampledAsMaxE,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "unsampledasmax",                      /* long option name */
                NULL,                           /* parameter name */
                "Set energy values in unsampled region to maximum energy from sampled region or to value provided by --maxenergy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MaxEnergy,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "maxenergy",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "If set, this is the energy used of unsampled regions.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoHeader,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "noheader",                      /* long option name */
                NULL,                           /* parameter name */
                "Do not print header to the output file.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                IncludeBinStat,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "includebinstat",                      /* long option name */
                NULL,                           /* parameter name */
                "Include bin statuses (1=sampled, 0=unsampled, -1=glued) into resulting FES.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IXFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fx",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT, which will be used to print values of collective variables.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OEFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT, which will be used to print values of free energy.")   /* option description */
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
                GPRNumDiff,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "numdiff",                      /* long option name */
                NULL,                           /* parameter name */
                "Use numerical differentiation of kernel function (for testing only).")   /* option description */
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
                "Calculate logPL.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRNoFastError,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "nofasterror",                      /* long option name */
                NULL,                           /* parameter name */
                "Do not use faster algorithm for error calculation.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                SLevel,                        /* option name */
                1.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "slevel",                      /* long option name */
                "VALUE",                           /* parameter name */
                "Sigma-level for confidence interval.")   /* option description */
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
    virtual int CheckArguments(void);
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
