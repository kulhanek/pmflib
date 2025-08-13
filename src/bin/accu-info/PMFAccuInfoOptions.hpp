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
    "                                               ** <b>info</b>              - print summary about the accumulator\n"
    "                                               ** <b>list-sections</b>     - print available sections\n"
    "                                               ** <b>get-section</b>       - get data from the section\n"
    "                                               ** <b>get-derivative</b>    - get derivative for given realm\n"
    "                                               ** <b>get-energy</b>        - get energy for given realm\n"
    "                                               ** <b>get-mean</b>          - get sample mean, variance, and mean error\n"
    "                                               ** <b>get-tseries</b>       - get time series\n"
    "                                               ** <b>nsamples</b>          - print number of samples in bins\n"
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
        CSO_MAP_OPT(int, Limit, 0, false, 'l', "limit", "LIMIT",
            "Report only bins containing more samples than the specified NUMBER."
            )
        // -------------------------------------------
        CSO_MAP_OPT(int, CV, 1, false, 'i', "item", "CV",
            "Select the collective variable (CV) for data output."
            )
        // -------------------------------------------
        CSO_MAP_OPT(bool, NoGNUPlot, false, false, 0, "nognuplot", NULL,
            "Disable printing of delimiters between records."
            )
        // -------------------------------------------
        CSO_MAP_OPT(bool, NoHeader, false, false, 0, "noheader", NULL,
            "Do not print the header in the output."
            )
        // -------------------------------------------
        CSO_MAP_OPT(CSmallString, IXFormat, "%15.7e", false, '\0', "fx", "FORMAT",
            "Format string for printing values of collective variables."
            )
        // -------------------------------------------
        CSO_MAP_OPT(CSmallString, OSFormat, "%15.7e", false, '\0', "fe", "FORMAT",
            "Format string for printing values of section data."
            )
        // -------------------------------------------
        CSO_MAP_OPT(bool, Verbose, false, false, 'v', "verbose", NULL,
            "Increase the verbosity of the output."
            )
        // -------------------------------------------
        CSO_MAP_OPT(bool, Version, false, false, '\0', "version", NULL,
            "Display version information and exit."
            )
        // -------------------------------------------
        CSO_MAP_OPT(bool, Help, false, false, 'h', "help", NULL,
            "Display this help message and exit."
            )
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
    virtual int CheckArguments(void);
};

//------------------------------------------------------------------------------

#endif
