#ifndef PMFAccuInfoOptionsH
#define PMFAccuInfoOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CPMFAccuInfoOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CPMFAccuInfoOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "accu-info"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "<b>accu-info</b> prints info about the PMF accumulator file and extracts data from it."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

    CSO_PROG_ARGS_SHORT_DESC_BEGIN
    "accuname action [section]"
    CSO_PROG_ARGS_SHORT_DESC_END

    CSO_PROG_ARGS_LONG_DESC_BEGIN
    "<cyan><b>accuname</b></cyan>                   Input name of PMF accumulator or '-' to read it from the standard input.\n"
    "<cyan><b>action</b></cyan>                     Requested action:\n"
    "                                               ** <b>info</b>          - print summary about the accumulator\n"
    "                                               ** <b>list-sections</b> - print available sections\n"
    "                                               ** <b>get-section</b>   - get data from the section\n"
    "                                               ** <b>nsamples</b>      - print number of samples in bins\n"
    "                                               ** <b>micf</b>          - print mean ICF forces (only the ABF simulations)"
    CSO_PROG_ARGS_LONG_DESC_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // options ------------------------------
    CSO_OPT(int,Limit)
    CSO_OPT(int,CV)
    CSO_OPT(bool,Sigma)
    CSO_OPT(bool,Errors)
    CSO_OPT(bool,NoGNUPlot)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OSFormat)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of options ---------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "LIMIT",                           /* parameter name */
                "Only bins containing more samples than NUMBER will be reported.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                CV,                        /* option name */
                1,                          /* default value */
                false,                          /* is option mandatory */
                'i',                           /* short option name */
                "item",                      /* long option name */
                "CV",                           /* parameter name */
                "Select CV for data printing.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoGNUPlot,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "nognuplot",                      /* long option name */
                NULL,                           /* parameter name */
                "Do not print delimiters between records.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoHeader,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "noheader",                      /* long option name */
                NULL,                           /* parameter name */
                "Do not print header to the output.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IXFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fx",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "FORMAT for printing values of collective variables.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OSFormat,                        /* option name */
                "%15.7e",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "FORMAT for printing values of section data.")   /* option description */
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
    virtual int CheckArguments(void);
};

//------------------------------------------------------------------------------

#endif
