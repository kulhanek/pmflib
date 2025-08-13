#ifndef ABFIntOptionsH
#define ABFIntOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

class CPMFEnergyIntOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CPMFEnergyIntOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "pmf-integrate"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program numerically integrates data from the PMF calculations. The integration is performed either by "
    "the reverse finite difference (RFD) method, radial basis functions (RBF), or Gaussian process (GPR)."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

    CSO_PROG_ARGS_SHORT_DESC_BEGIN
    "accuname1 [accuname2 [...]] fename"
    CSO_PROG_ARGS_SHORT_DESC_END

    CSO_PROG_ARGS_LONG_DESC_BEGIN
    "<cyan><b>accuname1</b></cyan>                  Name of file containing the ABF accumulator.\n"
    "<cyan><b>fename</b></cyan>                     Name of file where the resulting free energy surface will be printed. If the name is '-' then the output will be written to the standard output.\n"
    CSO_PROG_ARGS_LONG_DESC_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // options ------------------------------
    CSO_OPT(CSmallString,Realm)
    CSO_OPT(bool,ListRealms)
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
    CSO_OPT(CSmallString,SigmaF2)
    CSO_OPT(CSmallString,WFac)
    CSO_OPT(CSmallString,NCorr)
    CSO_OPT(CSmallString,SigmaN2)
    CSO_OPT(CSmallString,LoadHyprms)
    CSO_OPT(CSmallString,RFac)
    CSO_OPT(int,Overhang)
    CSO_OPT(bool,IncludeGluedRegions)
    CSO_OPT(int,GlueingFactor)
    CSO_OPT(bool,GlueHoles)
    CSO_OPT(bool,Periodicity)
    CSO_OPT(CSmallString,GlobalMin)
    CSO_OPT(double,Offset)
    CSO_OPT(bool,WithError)
    CSO_OPT(bool,NoEnergy)
    CSO_OPT(CSmallString,OutputFormat)
    CSO_OPT(CSmallString,PrintAll)
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
    CSO_OPT(bool,GPRNoFastError)
    CSO_OPT(double,SLevel)
    CSO_OPT(CSmallString,KeepCVs)
    CSO_OPT(CSmallString,ReducedFES)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
        CSO_MAP_OPT(CSmallString, Realm, "dG/dx", false, 'r', "realm", "NAME",
            "Requested realm for the integration. The list of supported realms can be obtained by --listrealms.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, ListRealms, false, false, '\0', "listrealms", NULL,
            "List supported realms.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, Method, "rfd", false, 'm', "method", "NAME",
            "Integration method. Supported methods are: rfd (reverse finite differences via csparse), "
            "rbf (radial basis functions), and gpr (Gaussian process).")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, EcutMethod, "rfd", false, '\0', "emethod", "NAME",
            "Integration method for the energy cut-off. Supported methods are: rfd (reverse finite differences via csparse), "
            "rbf (radial basis functions), and gpr (Gaussian process).")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, LAMethod, "default", false, 'a', "lmethod", "NAME",
            "Linear algebra method for LLS solution or matrix inversion. Supported algorithms are: "
            "default, svd (SVD - singular value decomposition, divide-and-conquer driver), "
            "svd2 (SVD - singular value decomposition, simple driver), qr (QR factorization), "
            "lu (LU factorization), ll (LL - Cholesky factorization). Possible combinations are: "
            "RFD(LU, default), RBF(QR, SVD, default), and GPR(LU, SVD, SVD2, LL, default).")
    // -------------------------------------------
        CSO_MAP_OPT(double, RCond, 1e-6, false, 0, "rcond", "NUMBER",
            "RBF+GPR: Rank condition for SVD. The selected value must be carefully tested. "
            "Calculation at computer precision is requested with -1 (not recommended).")
    // -------------------------------------------
        CSO_MAP_OPT(int, Limit, 1000, false, 'l', "limit", "NUMBER",
            "Only bins containing more samples than NUMBER are considered properly sampled.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, SkipFFTest, false, false, 0, "skipfftest", NULL,
            "Skip flood-fill test for discontinuous regions.")
    // -------------------------------------------
        CSO_MAP_OPT(double, EnergyLimit, -1.0, false, 'q', "energylimit", "NUMBER",
            "Integrate data only if the free energy is below NUMBER. "
            "This limit enforces two integration runs. A negative value disables the limit.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, EraseNegativeEnergy, false, false, '\0', "rmnegene", NULL,
            "Remove regions with negative free energy.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, SkipLastEnergyLimit, false, false, '\0', "skiplastelimit", NULL,
            "Skip energy limit filtering after the final integration.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MFMaxZScore, -1.0, false, 0, "mfmaxzscore", "NUMBER",
            "RBF+GPR: Reject mean forces whose prediction errors have a z-score above NUMBER. "
            "It is assumed that mean force errors have zero mean and follow a normal distribution. "
            "This limit is applied in each pass. A negative value disables the limit.")
    // -------------------------------------------
        CSO_MAP_OPT(int, MFZTestPasses, 1, false, 0, "mfnumofztests", "NUMBER",
            "RBF+GPR: Repeat the z-score test for mean force errors NUMBER times.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaF2, "15.0", false, 's', "sigmaf2", "NUMBER",
            "GPR: Variance of the reconstructed free energy surface (signal variance).")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, NCorr, "0.0", false, 'c', "ncorr", "NUMBER",
            "GPR: Number of statistically correlated samples.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, WFac, "3.0", false, 'w', "wfac", "SPEC",
            "RBF+GPR: Factors influencing widths of RBFs or squared exponential kernels. "
            "The width is the distance between adjacent squared exponential functions, multiplied by these factors, "
            "in the form WFac1[xWFac2x...]. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SigmaN2, "0.0", false, 'n', "sigman2", "SPEC",
            "GPR: Values of noise sigma squared for each CV, in the form SigmaN2(1)[xSigmaN2(2)x...]. "
            "The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, LoadHyprms, NULL, false, 0, "loadhyprms", "NAME",
            "GPR: Name of file containing GPR hyperparameters.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, RFac, "1.0", false, 't', "rfac", "SPEC",
            "RBF: Reduction factor for the number of RBFs. "
            "The number of RBFs in a given direction is the number of bins in that direction divided by this factor, "
            "in the form RFac1[xRFac2x...]. The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(int, Overhang, 2, false, 'g', "overhang", "NUMBER",
            "RBFs overhang to properly integrate areas near sampled edges. Ignored for periodic CVs.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, IncludeGluedRegions, 0, false, '\0', "includeglued", NULL,
            "RBF+GPR: Explicitly include glued regions. This option is set ON when --glueing > 0.")
    // -------------------------------------------
        CSO_MAP_OPT(int, GlueingFactor, 0, false, 0, "glueing", "NUMBER",
            "RBF+GPR: Calculate energy also for unsampled bins in close vicinity to sampled ones.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GlueHoles, 0, false, 0, "glueholes", NULL,
            "RBF+GPR: Calculate energy also for unsampled regions inside the FES.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Periodicity, false, false, 'p', "periodic", NULL,
            "RFD: Enable periodicity for collective variables that are periodic.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, GlobalMin, NULL, false, '\0', "globalmin", "SPEC",
            "RFD+RBF+GPR: Position of the global minimum provided as a single string in the form CV1xCV2x...xCVn "
            "(relevant for error determination). If not set, the position is determined automatically.")
    // -------------------------------------------
        CSO_MAP_OPT(double, Offset, 0.0, false, 'o', "offset", "NUMBER",
            "Specify an integration constant.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, WithError, false, false, 'e', "witherror", NULL,
            "GPR: Estimate free energy errors.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, NoEnergy, false, false, 0, "noenergy", NULL,
            "GPR: Skip calculation of energy and errors (can save time when only logML is required).")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, SaveABF, NULL, false, '\0', "saveabf", "NAME",
            "Save the final ABF accumulator to the file NAME.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, MFInfo, NULL, false, 0, "mfinfo", "NAME",
            "RBF+GPR: Name of file containing input and predicted mean forces.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, OutputFormat, "gnuplot", false, '\0', "output", "FORMAT",
            "Output FORMAT to print the free energy surface. Supported formats are: plain, gnuplot.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, PrintAll, false, NULL, '\0', "printall", "NAME",
            "Print results for all bins, even if they are not properly sampled.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, UnsampledAsMaxE, false, false, '\0', "unsampledasmax", NULL,
            "Set energy values in unsampled regions to the maximum energy from the sampled region, "
            "or to the value provided by --maxenergy.")
    // -------------------------------------------
        CSO_MAP_OPT(double, MaxEnergy, 0.0, false, '\0', "maxenergy", "NUMBER",
            "If set, this is the energy used for unsampled regions.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, NoHeader, false, false, '\0', "noheader", NULL,
            "Do not print a header in the output file.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, IncludeBinStat, false, false, '\0', "includebinstat", NULL,
            "Include bin statuses (1 = sampled, 0 = unsampled, -1 = glued) in the resulting FES.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, UseOldRFD, false, false, '\0', "oldrfd", NULL,
            "RFD: Use the old RFD implementation.")
    // -------------------------------------------
        CSO_MAP_OPT(int, FDPoints, 3, false, '\0', "fdpoints", "NUMBER",
            "RFD: Number of points used in the differentiation scheme "
            "(three or four points are supported) in the RFD method.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, IXFormat, "%15.7e", false, '\0', "fx", "FORMAT",
            "Output FORMAT to print values of collective variables.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, OEFormat, "%15.7e", false, '\0', "fe", "FORMAT",
            "Output FORMAT to print values of free energy.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, GPRKernel, "default", false, '\0', "kernel", "NAME",
            "GPR: Kernel type. Supported types: ardse (ARD squared exponential), "
            "ardmc52 (ARD Matern class 5/2), default (= ardse).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRNumDiff, false, false, '\0', "numdiff", NULL,
            "GPR: Use numerical differentiation of the kernel function (for testing only).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRUseInv, false, false, '\0', "useinv", NULL,
            "GPR: Use the matrix inversion pathway (for testing only).")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRCalcLogPL, false, false, '\0', "calclogpl", NULL,
            "GPR: Calculate logPL.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRIncludeZPE, false, false, '\0', "inczpe", NULL,
            "GPR: Include zero-point energy at the position specified by --globalmin.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, GPRNoFastError, false, false, '\0', "nofasterror", NULL,
            "GPR: Do not use the faster algorithm for error calculation.")
    // -------------------------------------------
        CSO_MAP_OPT(double, SLevel, 1.0, false, '\0', "slevel", "VALUE",
            "Sigma-level for the confidence interval.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, KeepCVs, "T", false, '\0', "keepcvs", "SPEC",
            "Which CVs should be kept during statistical reweighting of the FES. "
            "Flags are specified in the form CV1[xCV2x...] with F and T for skipped and kept CVs, respectively. "
            "The last value pads the rest.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, ReducedFES, NULL, false, '\0', "reducedfes", "NAME",
            "Name of file for FES reduced by statistical reweighting, containing only kept CVs.")
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
    virtual int CheckArguments(void);
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
