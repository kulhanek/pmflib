#ifndef EntropyDerOptionsH
#define EntropyDerOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

class CEntropyDerOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CEntropyDerOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "pmf-entropy-der"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The program provides entropy from the PMF accumulator containing time progress of system properties."
    CSO_PROG_DESC_END

    CSO_PROG_ARGS_SHORT_DESC_BEGIN
    "inaccu outaccu"
    CSO_PROG_ARGS_SHORT_DESC_END

    CSO_PROG_ARGS_LONG_DESC_BEGIN
    "<cyan><b>accuname1</b></cyan>                  Name of file containing the PMF accumulator.\n"
    "<cyan><b>outaccu</b></cyan>                    Resulting PMF accumulator.\n"
    CSO_PROG_ARGS_LONG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,InAccuName)
    CSO_ARG(CSmallString,OutAccuName)
    // options ------------------------------
    CSO_OPT(CSmallString,KSKernel)
    CSO_OPT(CSmallString,KSWFac)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                InAccuName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "inaccu",                        /* parameter name */
                "Name of input PMF accumulator. If the name is '-' then the PFM accumulator is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                OutAccuName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "outaccu",                        /* parameter name */
                "Name of resulting PMF accumulator. If the name is '-' then the PFM accumulator will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                KSKernel,                        /* option name */
                "default",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "kskernel",                      /* long option name */
                "NAME",                           /* parameter name */
                "Kernel for data smoothing. Supported types: epanechnikov, triweight, gaussian, and default(=triweight)")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                KSWFac,                        /* option name */
                "1.0",                          /* default value */
                false,                          /* is option mandatory */
                'w',                           /* short option name */
                "kswfac",                      /* long option name */
                "SPEC",                           /* parameter name */
                "Factors influencing kernel smoothing for data histogram construction. The CV width is distance between "
                "the adjacent bins multiplied by this factors in the form WFac1[xWFac2x...]. "
                "The last value pads the rest.")   /* option description */
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
