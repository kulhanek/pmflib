#ifndef ABFIntOptionsH
#define ABFIntOptionsH
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

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class CABFIntOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABFIntOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-integrate"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program numericaly integrates data from the ABF calculation. The integration is performed either by "
    "reverse finite difference (RFD or RFD2) method, employing radial basis functions (RBF), or gaussian process (GPR)."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

    CSO_PROG_ARGS_SHORT_DESC_BEGIN
    "accuname fename [fullfename]"
    CSO_PROG_ARGS_SHORT_DESC_END

    CSO_PROG_ARGS_LONG_DESC_BEGIN
    "<cyan><b>accuname</b></cyan>                   Name of file containing the ABF accumulator. If the name is '-' then the accumulator is read from the standard input.\n"
    "<cyan><b>fename</b></cyan>                     Name of file where the resulting free energy surface will be printed. If the name is '-' then the output will be written to the standard output.\n"
    "<cyan><b>fullfename</b></cyan>                 Optional name of file with free energy surface containing all rigions (sampled and unsampled)."
    CSO_PROG_ARGS_LONG_DESC_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN

    // options ------------------------------
    CSO_OPT(CSmallString,Method)
    CSO_OPT(CSmallString,EcutMethod)
    CSO_OPT(CSmallString,LAMethod)
    CSO_OPT(double,RCond)
    CSO_OPT(int,Limit)
    CSO_OPT(bool,SkipFFTest)
    CSO_OPT(double,EnergyLimit)
    CSO_OPT(bool,EraseNegativeEnergy)
    CSO_OPT(bool,SkipLastEnergyLimit)
    CSO_OPT(double,MFMaxZScore)
    CSO_OPT(int,MFZTestPasses)
    CSO_OPT(double,SigmaF2)
    CSO_OPT(CSmallString,NCorr)
    CSO_OPT(bool,SplitNCorr)
    CSO_OPT(CSmallString,WFac)
    CSO_OPT(CSmallString,LoadHyprms)
    CSO_OPT(CSmallString,RFac)
    CSO_OPT(int,Overhang)
    CSO_OPT(bool,IncludeGluedRegions)
    CSO_OPT(int,GlueingFactor)
    CSO_OPT(bool,GlueHoles)
    CSO_OPT(bool,Periodicity)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(bool,UseRealGlobalMin)
    CSO_OPT(double,Offset)
    CSO_OPT(bool,WithError)
    CSO_OPT(bool,NoEnergy)
    CSO_OPT(CSmallString,OutputFormat)
    CSO_OPT(bool,PrintAll)
    CSO_OPT(bool,UnsampledAsMaxE)
    CSO_OPT(double,MaxEnergy)   
    CSO_OPT(bool,NoHeader)
    CSO_OPT(bool,IncludeBinStat)
    CSO_OPT(bool,UseOldRFD)
    CSO_OPT(int,FDPoints)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OEFormat)
    CSO_OPT(CSmallString,MFInfo)
    CSO_OPT(CSmallString,SaveABF)
    CSO_OPT(CSmallString,GPRKernel)
    CSO_OPT(bool,GPRNumDiff)
    CSO_OPT(bool,GPRUseInv)
    CSO_OPT(bool,GPRCalcLogPL)
    CSO_OPT(bool,GPRIncludeZPE)
    CSO_OPT(bool,GPRFastError)
    CSO_OPT(CSmallString,KeepCVs)
    CSO_OPT(CSmallString,ReducedFES)
    CSO_OPT(double,Temperature)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of options ---------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Method,                        /* option name */
                "rfd",                          /* default value */
                false,                          /* is option mandatory */
                'm',                           /* short option name */
                "method",                      /* long option name */
                "NAME",                           /* parametr name */
                "Integration method. Supported methods are: rfd (reverse finite differences via csparse), rfd2 (reverse finite differences via lapack), rbf (radial basis functions), and gpr (gaussian process).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                EcutMethod,                        /* option name */
                "rfd",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "emethod",                      /* long option name */
                "NAME",                           /* parametr name */
                "Integration method for energy cut-off. Supported methods are: rfd (reverse finite differences via csparse), rfd2 (reverse finite differences via lapack), rbf (radial basis functions), and gpr (gaussian process).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                LAMethod,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                'a',                           /* short option name */
                "lmethod",                      /* long option name */
                "NAME",                           /* parametr name */
                "Linear algebra method for LLS solution or matrix inversion. Supported algorithms are: "
                "default, svd (SVD - singular value decomposition, divide and conquer driver), "
                "svd2 (SVD - singular value decomposition, simple driver), qr (QR factorization), "
                "lu (LU factorization), ll (LL - Cholesky factorization). "
                "Possible combinations are: RFD(LU,default), RFD2(QR,SVD,default), "
                "RBF(QR,SVD,default), and GPR(LU,SVD,SVD2,LL,default).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                RCond,                        /* option name */
                1e-6,                          /* default value */
                false,                          /* is option mandatory */
                'r',                           /* short option name */
                "rcond",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RFD2+RBF+GPR: Rank condition for SVD. Used value must be carefully tested. Calculation at computer precision is requested with -1 (not recommended).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                100,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                SkipFFTest,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "skipfftest",                      /* long option name */
                NULL,                           /* parametr name */
                "Skip flood fill test for discountinous regions.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                EnergyLimit,                        /* option name */
                -1.0,                          /* default value */
                false,                          /* is option mandatory */
                'q',                           /* short option name */
                "energylimit",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Integrate data only if the free energy is below NUMBER. "
                "This limit enforces two integration runs. Negative value disables the limit.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                EraseNegativeEnergy,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "rmnegene",                      /* long option name */
                NULL,                           /* parametr name */
                "Remove regions with negative free energy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                SkipLastEnergyLimit,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "skiplastelimit",                      /* long option name */
                NULL,                           /* parametr name */
                "Skip energy limit filter after final integration.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MFMaxZScore,                        /* option name */
                -1.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "mfmaxzscore",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBF+GPR: Reject mean forces whose errors in prediction have z-score above NUMBER. In the test, it is assumend that mean force errors have zero mean and they follow normal distribution. "
                "This limit is aplied in the each pass. Negative value disables the limit.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                MFZTestPasses,                        /* option name */
                1,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "mfnumofztests",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBF+GPR: Repeat z-score test for mean force errrors NUMBER times.")   /* option description */
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
    CSO_MAP_OPT(CSmallString,                           /* option type */
                NCorr,                        /* option name */
                "1.0",                          /* default value */
                false,                          /* is option mandatory */
                'c',                           /* short option name */
                "ncorr",                      /* long option name */
                "SPEC",                           /* parametr name */
                "GPR: Number of statistically correlated samples in the form NCorr1[xNCorr2x...]. "
                "The last value pads the rest. Split ncorr mode is enabled by the --splitncorr option.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                SplitNCorr,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "splitncorr",                      /* long option name */
                NULL,                           /* parametr name */
                "Use indepenedent ncorr for each CV.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                WFac,                        /* option name */
                "3.0",                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "wfac",                      /* long option name */
                "SPEC",                           /* parametr name */
                "RBF+GPR: Factors influencing widths of RBFs or square exponential kernels. The width is distance between "
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
    CSO_MAP_OPT(CSmallString,                           /* option type */
                RFac,                        /* option name */
                "1.0",                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "rfac",                      /* long option name */
                "SPEC",                           /* parametr name */
                "RBF: Reduction factor for number of RBFs. Number of RBFs in given direction is number of bins in that "
                "direction divided by this factor in the form RFac1[xRFac2x...]. The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Overhang,                        /* option name */
                2,                          /* default value */
                false,                          /* is option mandatory */
                'g',                           /* short option name */
                "overhang",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBFs overhang to properly integrate areas near sampled edges. Ignored for periodic CVs.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                IncludeGluedRegions,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "includeglued",                      /* long option name */
                NULL,                           /* parametr name */
                "RBF+GPR: Explicitly include glued regions. This options is set ON when --glueing > 0.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                GlueingFactor,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "glueing",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBF+GPR: Calculate energy also for unsampled bins in close vicinity to sampled ones.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GlueHoles,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "glueholes",                      /* long option name */
                NULL,                           /* parametr name */
                "RBF+GPR: Calculate energy also for unsampled regions inside the FES.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Periodicity,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "periodic",                      /* long option name */
                NULL,                           /* parametr name */
                "RFD: Switch on periodicity for collective variables that are periodic.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                GlobalMin,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "globalmin",                      /* long option name */
                "SPEC",                           /* parametr name */
                "RBF+GPR: position of global minimum provided as a single string in the form CV1xCV2x...xCVn (relevant for error determination), if not set the position is determined automatically.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                UseRealGlobalMin,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "userealglbmin",                      /* long option name */
                NULL,                           /* parametr name */
                "RBF+GPR: ignore --globalmin in the final integration step")   /* option description */
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
    CSO_MAP_OPT(bool,                           /* option type */
                WithError,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'e',                           /* short option name */
                "witherror",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Estimate free energy errors.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoEnergy,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "noenergy",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Skip calculation of energy and errors (it can save some time when only logML is required).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SaveABF,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "saveabf",                      /* long option name */
                "NAME",                           /* parametr name */
                "Save the final ABF accumulator into the file with NAME.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                MFInfo,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "mfinfo",                      /* long option name */
                "NAME",                           /* parametr name */
                "RBF+GPR: name of file with input and predicted mean forces.")   /* option description */
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
                PrintAll,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "printall",                      /* long option name */
                NULL,                           /* parametr name */
                "Print results for all bins even if not properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                UnsampledAsMaxE,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "unsampledasmax",                      /* long option name */
                NULL,                           /* parametr name */
                "Set energy values in unsampled region to maximum energy from sampled region or to value provided by --maxenergy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MaxEnergy,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "maxenergy",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "If set, this is the energy used of unsampled regions.")   /* option description */
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
                IncludeBinStat,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "includebinstat",                      /* long option name */
                NULL,                           /* parametr name */
                "Include bin statuses (1=sampled, 0=unsampled, -1=glued) into resulting FES.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                UseOldRFD,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "oldrfd",                      /* long option name */
                NULL,                           /* parametr name */
                "RFD: Use old RFD implementation.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                FDPoints,                        /* option name */
                3,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "fdpoints",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RFD: Determine number of points employed in differenciation scheme (three or four is upported) in RFD method.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IXFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fx",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print values of collective variables.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OEFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print values of free energy.")   /* option description */
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
    CSO_MAP_OPT(bool,                           /* option type */
                GPRNumDiff,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "numdiff",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Use numerical differentiation of kernel function (for testing only).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRUseInv,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "useinv",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Use matrix inversion pathway (for testing only).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRCalcLogPL,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "calclogpl",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Calculate logPL.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRIncludeZPE,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "inczpe",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Include zero-point energy at position specified by --globalmin.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                GPRFastError,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "fasterror",                      /* long option name */
                NULL,                           /* parametr name */
                "GPR: Use faster algorithm for error calculation.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                KeepCVs,                        /* option name */
                "T",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "keepcvs",                      /* long option name */
                "SPEC",                           /* parametr name */
                "Which CVs should be kept during statistical reweighting of FES. Flags are specified in the form CV1[xCV2x...] with F and T for skip and kept CV, respectively. "
                "The last value pads the rest.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                ReducedFES,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "reducedfes",                      /* long option name */
                "NAME",                           /* parametr name */
                "Name of file for FES reduced by statistical reweighting containing only kept CVs.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Temperature,                        /* option name */
                300,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "temperature",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Absolute temperature for reduction of FES in [K].")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Verbose,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'v',                           /* short option name */
                "verbose",                      /* long option name */
                NULL,                           /* parametr name */
                "Increase output verbosity.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Version,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "version",                      /* long option name */
                NULL,                           /* parametr name */
                "Output version information and exit.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Help,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'h',                           /* short option name */
                "help",                      /* long option name */
                NULL,                           /* parametr name */
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
