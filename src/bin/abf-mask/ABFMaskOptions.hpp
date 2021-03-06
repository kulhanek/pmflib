#ifndef ABFMaskOptionsH
#define ABFMaskOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CABFMaskOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABFMaskOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-mask"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "It create the ABF accumulator mask according to user specified criteria."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,ABFAccuName)
    CSO_ARG(CSmallString,MaskOutputName)
    // options ------------------------------
    CSO_OPT(double,MaxEnergy)
    CSO_OPT(int,Limit)
    CSO_OPT(int,FDPoints)
    CSO_OPT(double,Offset)
    CSO_OPT(bool,Periodicity)
    CSO_OPT(CSmallString,OutputFormat)
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
                MaskOutputName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "mask",                        /* parametr name */
                "Name of file where the info about samples will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MaxEnergy,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                'm',                           /* short option name */
                "maxene",                      /* long option name */
                "ENERGY",                           /* parametr name */
                "Generate mask with zero weights for bins with the free energy greater than ENERGY.")   /* option description */
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
    CSO_MAP_OPT(int,                           /* option type */
                FDPoints,                        /* option name */
                3,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "fdpoints",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Determines differenciation scheme. Currently implemented schemes use either three or four points.")   /* option description */
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
                Periodicity,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "periodic",                      /* long option name */
                NULL,                           /* parametr name */
                "Switch on periodicity for collective variables that are periodic.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OutputFormat,                        /* option name */
                "mask",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "output",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print free energy surface. Supported formats are: mask, gnuplot.")   /* option description */
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
