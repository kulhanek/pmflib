#ifndef PMFEnergyOptionsH
#define PMFEnergyOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CPMFEneOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CPMFEneOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "pmf-energy"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program provides various types of energies from the PMF accumulator."
    CSO_PROG_DESC_END

    CSO_PROG_ARGS_SHORT_DESC_BEGIN
    "accuname1 [accuname2 ...] energy"
    CSO_PROG_ARGS_SHORT_DESC_END

    CSO_PROG_ARGS_LONG_DESC_BEGIN
    "<cyan><b>accuname1</b></cyan>                  Name of file containing the PMF accumulator.\n"
    "<cyan><b>energy</b></cyan>                     Resulting energy.\n"
    CSO_PROG_ARGS_LONG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // options ------------------------------
    CSO_OPT(CSmallString,Realm)
    CSO_OPT(bool,ListRealms)
    CSO_OPT(CSmallString,Method)
    CSO_OPT(int,Limit)
    CSO_OPT(bool,UnsampledAsMaxE)
    CSO_OPT(double,MaxEnergy)
    CSO_OPT(bool,IncludeBinStat)
    CSO_OPT(bool,Absolute)
    CSO_OPT(bool,WithError)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(double,Offset)
    CSO_OPT(CSmallString,GPRKernel)
    CSO_OPT(bool,GPRCalcLogPL)
    CSO_OPT(CSmallString,SigmaF2)
    CSO_OPT(CSmallString,WFac)
    CSO_OPT(CSmallString,NCorr)
    CSO_OPT(CSmallString,SigmaN2)
    CSO_OPT(CSmallString,LoadHyprms)
    CSO_OPT(double,SLevel)
    CSO_OPT(CSmallString,MFInfo)
    CSO_OPT(CSmallString,OutputFormat)
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
        CSO_MAP_OPT(CSmallString, Realm, "<Eint>", false, 'r', "realm", "NAME",
                "Intended output. The list of supported realms can be obtained by --listrealms.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, ListRealms, false, false, '\0', "listrealms", NULL,
                "List supported realms.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, Method, "raw", false, 'm', "method", "NAME",
                "Supported methods: raw (data taken directly from the accumulator) and gpr (Gaussian process filtered data).")
        // -------------------------------------------
            CSO_MAP_OPT(bool, Absolute, false, false, 'a', "absolute", NULL,
                "Absolute energy.")
        // -------------------------------------------
            CSO_MAP_OPT(int, Limit, 0, false, 'l', "limit", "LIMIT",
                "Only bins containing more samples than NUMBER are considered properly sampled.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, WithError, false, false, 'e', "witherror", NULL,
                "GPR: Estimate energy errors from the GPR model. RAW: Print energy errors from the PMF accumulator.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, GPRKernel, "default", false, '\0', "kernel", "NAME",
                "GPR: Kernel type. Supported types: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), "
                "ardmc32 (ARD Matern class 3/2), ardmc12 (ARD Matern class 1/2), default (=ardse).")
        // -------------------------------------------
            CSO_MAP_OPT(bool, GPRCalcLogPL, false, false, 0, "calclogpl", NULL,
                "GPR: Calculate logPL.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, SigmaF2, "15.0", false, 's', "sigmaf2", "NUMBER",
                "GPR: Variance of the reconstructed energy surface (signal variance).")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, WFac, "3.0", false, 'w', "wfac", "SPEC",
                "GPR: Factors influencing widths of squared exponential kernels. The width is the distance between "
                "adjacent squared exponential functions multiplied by these factors in the form WFac1[xWFac2x...]. "
                "The last value pads the rest.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, NCorr, "0.0", false, 'c', "ncorr", "VALUE",
                "Number of statistically correlated samples in the form NCorr1[NCorr2x...]. "
                "The last value pads the rest.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, SigmaN2, "0.0", false, 'n', "sigman2", "SPEC",
                "Values of noise sigma squared for each CV in the form SigmaN2(1)[xSigmaN2(2)x...]. "
                "The last value pads the rest.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, LoadHyprms, NULL, false, 0, "loadhyprms", "NAME",
                "GPR: Name of file containing the GPR hyperparameters.")
        // -------------------------------------------
            CSO_MAP_OPT(double, SLevel, 1.0, false, 0, "slevel", "VALUE",
                "Sigma-level for the confidence interval.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, MFInfo, NULL, false, '\0', "mfinfo", "NAME",
                "GPR: Name of file with input and predicted energy.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, GlobalMin, NULL, false, '\0', "globalmin", "SPEC",
                "GPR: Position of the global minimum provided as a single string in the form CV1xCV2x...xCVn "
                "(relevant for error determination). If not set, the position is determined automatically.")
        // -------------------------------------------
            CSO_MAP_OPT(double, Offset, 0.0, false, 'o', "offset", "NUMBER",
                "Specify the integration constant.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, UnsampledAsMaxE, false, false, 0, "unsampledasmax", NULL,
                "Set energy values in unsampled regions to the maximum energy from sampled regions "
                "or to the value provided by --maxenergy.")
        // -------------------------------------------
            CSO_MAP_OPT(double, MaxEnergy, 0.0, false, 0, "maxenergy", "NUMBER",
                "If set, this is the energy used for unsampled regions.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, IncludeBinStat, false, false, 0, "includebinstat", NULL,
                "Include bin statuses (1 = sampled, 0 = unsampled, -1 = glued) in the resulting FES.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, OutputFormat, "gnuplot", false, 0, "output", "FORMAT",
                "Output FORMAT for printing the energy surface. Supported formats: plain and gnuplot.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, NoHeader, false, false, 0, "noheader", NULL,
                "Do not print a header to the output file.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, PrintAll, false, false, 0, "printall", NULL,
                "Print results for all bins, even if not properly sampled.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, IXFormat, "%15.7e", false, '\0', "fx", "FORMAT",
                "Output FORMAT for printing values of collective variables.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, OEFormat, "%15.7e", false, '\0', "fe", "FORMAT",
                "Output FORMAT for printing results.")
        // -------------------------------------------
            CSO_MAP_OPT(CSmallString, LAMethod, "default", false, 0, "lmethod", "NAME",
                "GPR: Linear algebra method for LLS solution or matrix inversion. Supported algorithms: "
                "svd (SVD - singular value decomposition, divide and conquer driver), "
                "svd2 (SVD - singular value decomposition, simple driver), "
                "lu (LU factorization), ll (LL - Cholesky factorization). "
                "Possible combinations: GPR(LU, SVD, SVD2, LL, default).")
        // -------------------------------------------
            CSO_MAP_OPT(double, RCond, 1e-6, false, 0, "rcond", "NUMBER",
                "GPR: Rank condition for SVD. The chosen value must be carefully tested. "
                "Calculation at computer precision is requested with -1 (not recommended).")
        // -------------------------------------------
            CSO_MAP_OPT(bool, Verbose, false, false, 'v', "verbose", NULL,
                "Increase output verbosity.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, Version, false, false, '\0', "version", NULL,
                "Output version information and exit.")
        // -------------------------------------------
            CSO_MAP_OPT(bool, Help, false, false, 'h', "help", NULL,
                "Display this help and exit.")
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
