#ifndef ABFEnthalpyOptionsH
#define ABFEnthalpyOptionsH
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

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class CABFEnthalpyOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABFEnthalpyOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-enthalpy"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program provides enthalpy from the ABF accumulator."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,ABFAccuName)
    CSO_ARG(CSmallString,OutputName)
    // options ------------------------------
    CSO_OPT(CSmallString,Method)
    CSO_OPT(int,Limit)
    CSO_OPT(bool,Absolute)
    CSO_OPT(bool,WithError)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(double,Offset)
    CSO_OPT(CSmallString,GPRKernel)
    CSO_OPT(double,SigmaF2)
    CSO_OPT(double,NCorr)
    CSO_OPT(CSmallString,WFac)
    CSO_OPT(CSmallString,LoadHyprms)
    CSO_OPT(CSmallString,OutputFormat)
    CSO_OPT(double,SLevel)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(bool,PrintAll)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OEFormat)
    CSO_OPT(CSmallString,LAMethod)
    CSO_OPT(double,RCond)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                ABFAccuName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "accuname",                        /* parameter name */
                "Name of file containing the ABF accumulator. If the name is '-' then the accumulator is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                OutputName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "enthalpy",                        /* parameter name */
                "Name of file where results will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Method,                        /* option name */
                "raw",                          /* default value */
                false,                          /* is option mandatory */
                'm',                           /* short option name */
                "method",                      /* long option name */
                "NAME",                           /* parametr name */
                "Supported methods are: raw (data taken directly from ABF accumulator) and gpr (gaussian process filtered data).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Absolute,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'a',                           /* short option name */
                "absolute",                      /* long option name */
                NULL,                           /* parametr name */
                "absolute enthalpy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "LIMIT",                           /* parameter name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                WithError,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'e',                           /* short option name */
                "witherror",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Estimate free energy errors. RAW: Print free energy errors from ABF accumulator.")   /* option description */
   //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                GPRKernel,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "kernel",                      /* long option name */
                "NAME",                           /* parametr name */
                "GPR: Kernel type. Supported types: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), default(=ardse)")   /* option description */
//----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                SigmaF2,                        /* option name */
                15.0,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "sigmaf2",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "GPR: Variance of the reconstructed free energy surface (signal variance).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                NCorr,                        /* option name */
                1.0,                          /* default value */
                false,                          /* is option mandatory */
                'c',                           /* short option name */
                "ncorr",                      /* long option name */
                "VALUE",                           /* parametr name */
                "Number of statistically correlated samples.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                WFac,                        /* option name */
                "3.0",                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "wfac",                      /* long option name */
                "SPEC",                           /* parametr name */
                "GPR: Factors influencing widths of square exponential kernels. The width is distance between "
                "the adjacent square exponential functions multiplied by this factors in the form WFac1[xWFac2x...]. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                LoadHyprms,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "loadhyprms",                      /* long option name */
                "NAME",                           /* parametr name */
                "GPR: Name of file with GPR hyperparameters.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                SLevel,                        /* option name */
                1.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "slevel",                      /* long option name */
                "VALUE",                           /* parametr name */
                "Sigma-level for confidence interval.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                GlobalMin,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "globalmin",                      /* long option name */
                "SPEC",                           /* parametr name */
                "GPR: position of global minimum provided as a single string in the form CV1xCV2x...xCVn (relevant for error determination), if not set the position is determined automatically.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Offset,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                'o',                           /* short option name */
                "offset",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Specify an integration constant.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OutputFormat,                        /* option name */
                "gnuplot",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "output",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT to print the free energy surface. Supported formats are: plain, gnuplot, fes.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoHeader,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "noheader",                      /* long option name */
                NULL,                           /* parametr name */
                "Do not print header to the output file.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                PrintAll,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "printall",                      /* long option name */
                NULL,                           /* parametr name */
                "Print results for all bins even if not properly sampled.")   /* option description */
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
                "Output FORMAT, which will be used to print results.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                LAMethod,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "lmethod",                      /* long option name */
                "NAME",                           /* parametr name */
                "Linear algebra method for LLS solution or matrix inversion. Supported algorithms are: "
                "default, svd (SVD - singular value decomposition, divide and conquer driver), "
                "svd2 (SVD - singular value decomposition, simple driver), "
                "lu (LU factorization), ll (LL - Cholesky factorization). "
                "Possible combinations are: GPR(LU,SVD,SVD2,LL,default).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                RCond,                        /* option name */
                1e-6,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "rcond",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "GPR: Rank condition for SVD. Used value must be carefully tested. Calculation at computer precision is requested with -1 (not recommended).")   /* option description */
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
