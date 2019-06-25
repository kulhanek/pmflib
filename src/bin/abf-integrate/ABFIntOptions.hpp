#ifndef ABFIntOptionsH
#define ABFIntOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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
    "It numericaly integrates data from the ABF calculation. The integration is performed either by "
    "reverse finite difference (RFD) method, employing radial basis functions (RBF), or gaussian process (GPR).\n"
    "<b><red>Warning:</red></b> "
    "While RFD is extremly fast it can fail on data with irregular sampling. This problem can be partially overcome "
    "by increasing sampling limit. "
    "Both RFD and RBF can fail on data with quick changes in mean forces. "
    "This situation can be partially overcome by increasing sampling limit and "
    "at the same time by energy limit, which, however, requires two integration passes. "
    "GPR is probably ultimate way how to obtain correctly integrated data but "
    "at cost of significantly higher computational demands."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,ABFAccuName)
    CSO_ARG(CSmallString,FEOutputName)
    // options ------------------------------
    CSO_OPT(int,Limit)
    CSO_OPT(double,EnergyLimit)
    CSO_OPT(double,Offset)
    CSO_OPT(CSmallString,Method)
    CSO_OPT(int,FDPoints)
    CSO_OPT(double,RCond)
    CSO_OPT(double,RFac)
    CSO_OPT(double,WFac)
    CSO_OPT(double,SigmaF2)
    CSO_OPT(int,Overhang)
    CSO_OPT(bool,Periodicity)
    CSO_OPT(bool,WithError)
    CSO_OPT(CSmallString,OutputFormat)
    CSO_OPT(bool,PrintWithLimit)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OEFormat)
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
                "accuname",                        /* parametr name */
                "Name of file containing the ABF accumulator. If the name is '-' then the accumulator is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                FEOutputName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "fename",                        /* parametr name */
                "Name of file where the resulting free energy surface will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                100,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
// description of options ---------------------------------------------------
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
                Method,                        /* option name */
                "rfd",                          /* default value */
                false,                          /* is option mandatory */
                'm',                           /* short option name */
                "method",                      /* long option name */
                "NAME",                           /* parametr name */
                "Integration method. Supported methods are: rfd (reverse finite differences), rbf (radial basis functions), and gpr (gaussian process).")   /* option description */
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
    CSO_MAP_OPT(double,                           /* option type */
                RCond,                        /* option name */
                -1.0,                          /* default value */
                false,                          /* is option mandatory */
                'r',                           /* short option name */
                "rcond",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBF: Rank condition for SVD.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                RFac,                        /* option name */
                3.0,                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "rfac",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBF: Reduction factor for number of RBFs. Number of RBFs in given direction is number of bins in that direction divided by this factor.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                WFac,                        /* option name */
                2.0,                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "wfac",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "RBF+GPR: Factor influencing widths of RBFs or square exponential kernels. The width is distance between "
                "the adjacent square exponential functions multiplied by this factor.")   /* option description */
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
                Periodicity,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "periodic",                      /* long option name */
                NULL,                           /* parametr name */
                "RFD: Switch on periodicity for collective variables that are periodic.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                WithError,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'e',                           /* short option name */
                "witherror",                      /* long option name */
                NULL,                           /* parametr name */
                "Integrate free energy and try to estimate its error.")   /* option description */
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
                PrintWithLimit,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "printwithlimit",                      /* long option name */
                NULL,                           /* parametr name */
                "Print results from bins that were adequately sampled.")   /* option description */
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
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
};

//------------------------------------------------------------------------------

#endif
