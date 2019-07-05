#ifndef ABFDeriOptionsH
#define ABFDeriOptionsH
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

class CABFDeriOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABFDeriOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-derivatives"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "<b>abf-derivatives</b> prints selected derivatives of the free energy with respect to "
    "given collective variables (mean forces) from the ABF accumulator."
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
    CSO_OPT(int,Limit)
    CSO_OPT(double,NCorr)
    CSO_OPT(int,Item)
    CSO_OPT(bool,Sigma)
    CSO_OPT(bool,Error)
    CSO_OPT(bool,NoGNUPlot)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OSFormat)
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
                OutputName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "dname",                        /* parametr name */
                "Name of file where results will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "LIMIT",                           /* parametr name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                NCorr,                        /* option name */
                1.0,                          /* default value */
                false,                          /* is option mandatory */
                'c',                           /* short option name */
                "ncorr",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Number of statistically correlated samples.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Item,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'i',                           /* short option name */
                "item",                      /* long option name */
                "CV",                           /* parametr name */
                "Forces for all CVs are printed (item == 0) unless item is set for given CV index (counted from 1).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Sigma,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "sigma",                      /* long option name */
                NULL,                           /* parametr name */
                "Print the standard deviation of instanteous forces (fluctuations).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Error,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'e',                           /* short option name */
                "error",                      /* long option name */
                NULL,                           /* parametr name */
                "Print the standard deviation of mean forces (standard error).")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoGNUPlot,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "nognuplot",                      /* long option name */
                NULL,                           /* parametr name */
                "Do not print delimiters between records.")   /* option description */
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
                OSFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print results.")   /* option description */
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
