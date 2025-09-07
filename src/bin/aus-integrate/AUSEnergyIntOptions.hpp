#ifndef AUSEnergyIntOptionsH
#define AUSEnergyIntOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CAUSEnergyIntOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CAUSEnergyIntOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "aus-integrate"
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
    CSO_ARG(CSmallString,FENFile)
    CSO_ARG(CSmallString,INTFile)
    CSO_ARG(CSmallString,TDSFile)
    // options ------------------------------
    CSO_OPT(CSmallString,Realm)
    CSO_OPT(bool,ListRealms)
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
        CSO_MAP_ARG(CSmallString, AccuFile, NULL, true, "ACCU",
            "Name of the file containing the input ABF accumulator.")
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, FENFile, NULL, true, "FEN",
            "Name of the file containing the output free energy surface [dA(x)].")
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, INTFile, NULL, true, "INT",
            "Name of the file containing the output internal energy surface [dU(x)].")
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, TDSFile, NULL, true, "TDS",
            "Name of the file containing the output entropic energy surface [-TdS(x)].")
    // -------------------------------------------
        CSO_MAP_OPT(bool, EnableConstraints, false, false, 0, "constraints", NULL,
            "Impose constraints: dA(x)/dx - dU(x)/dx - (-TdS(x)/dx) = 0.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, Realm, "AUS", false, 'r', "realm", "NAME",
            "Requested realm for data processing. The list of supported realms can be obtained by --listrealms.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, ListRealms, false, false, '\0', "listrealms", NULL,
            "List supported realms.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, LAMethod, "default", false, 'a', "lmethod", "NAME",
            "Linear algebra method for LLS solution or matrix inversion. Supported algorithms: "
            "default, svd (SVD - divide and conquer), svd2 (SVD - simple), qr (QR factorization), "
            "lu (LU factorization), ll (Cholesky factorization). "
            "Possible combinations: GPR(LU, SVD, SVD2, LL, default).")
    // -------------------------------------------
        CSO_MAP_OPT(double, RCond, 1e-6, false, 0, "rcond", "NUMBER",
            "Rank condition for SVD. The chosen value must be carefully tested. "
            "Calculation at machine precision is requested with -1 (not recommended).")
    // -------------------------------------------
        CSO_MAP_OPT(int, Limit, 1000, false, 'l', "limit", "NUMBER",
            "Only bins containing more samples than NUMBER are considered properly sampled.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaF2, "15.0", false, 's', "sigmaf2", "SPEC",
            "Signal variances for each realm in the form SigmaF2(1)[xSigmaF2(2)x...]. "
            "The last value pads the rest. There are three realms.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, CoVar, "0.0", false, 'c', "covar", "SPEC",
            "Signal covariances for each realm in the form CoVar(1)[xCoVar(2)x...]. "
            "The last value pads the rest. There are three realms.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, WFac, "3.0", false, 'w', "wfac", "SPEC",
            "Factors influencing widths of RBFs or squared exponential kernels. "
            "The width is the distance between adjacent functions multiplied by these factors, "
            "in the form WFac1[xWFac2x...]. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaN2, "0.0", false, 'n', "sigman2", "SPEC",
            "Noise variances (sigma squared) for each CV in the form SigmaN2(1)[xSigmaN2(2)x...]. "
            "The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, LoadHyprms, NULL, false, 0, "loadhyprms", "NAME",
            "Name of the file containing GPR hyperparameters.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, GlobalMin, NULL, false, '\0', "globalmin", "SPEC",
            "Position of the global minimum provided as a single string in the form CV1xCV2x...xCVn "
            "(relevant for error determination). If not set, the position is determined automatically.")
    // -------------------------------------------
        CSO_MAP_OPT(double, Offset, 0.0, false, 'o', "offset", "NUMBER",
            "Specify the integration constant.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, WithError, false, false, 'e', "witherror", NULL,
            "Estimate free energy errors.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, NoEnergy, false, false, 0, "noenergy", NULL,
            "Skip calculation of energy and errors (saves time when only logML is required).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, BalanceResiduals, false, false, 0, "balres", NULL,
            "Balance residual errors between dG, dH, and -TdS.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, MFInfo, NULL, false, '\0', "mfinfo", "NAME",
            "Name of the file with input and predicted mean forces.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, OutputFormat, "gnuplot", false, 0, "output", "FORMAT",
            "Output format for printing the free energy surface. Supported formats: plain, gnuplot.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, UnsampledAsMaxE, false, false, 0, "unsampledasmax", NULL,
            "Set energy values in unsampled regions to the maximum energy from sampled regions "
            "or to the value provided by --maxenergy.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MaxEnergy, 0.0, false, 0, "maxenergy", "NUMBER",
            "Energy value used for unsampled regions if set.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, NoHeader, false, false, 0, "noheader", NULL,
            "Do not print a header in the output file.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, IncludeBinStat, false, false, 0, "includebinstat", NULL,
            "Include bin statuses (1=sampled, 0=unsampled, -1=glued) in the resulting FES.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, IXFormat, "%15.7e", false, '\0', "fx", "FORMAT",
            "Output format used to print values of collective variables.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, OEFormat, "%15.7e", false, '\0', "fe", "FORMAT",
            "Output format used to print values of free energy.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, GPRKernel, "default", false, '\0', "kernel", "NAME",
            "GPR kernel type. Supported types: ardse (ARD squared exponential), "
            "ardmc52 (ARD Matern class 5/2), default (=ardse).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRNumDiff, false, false, 0, "numdiff", NULL,
            "Use numerical differentiation of the kernel function (for testing only).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRUseInv, false, false, 0, "useinv", NULL,
            "Use matrix inversion pathway (for testing only).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRCalcLogPL, false, false, 0, "calclogpl", NULL,
            "Calculate logPL.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRNoFastError, false, false, 0, "nofasterror", NULL,
            "Do not use the faster algorithm for error calculation.")
    // -------------------------------------------
        CSO_MAP_OPT(double, SLevel, 1.0, false, 0, "slevel", "VALUE",
            "Sigma level for confidence intervals.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Verbose, false, false, 'v', "verbose", NULL,
            "Increase output verbosity.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Version, false, false, '\0', "version", NULL,
            "Output version information and exit.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Help, false, false, 'h', "help", NULL,
            "Display this help message and exit.")
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckArguments(void);
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
