#ifndef ACCUCombOptionsH
#define ACCUCombOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CACCUCombOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CACCUCombOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "accu-combine"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "<b>accu-combine</b> combines more PMF accumulators into one. All PMF accumulators have to be of the same dimensions and contain data from the same PMF method."
    CSO_PROG_DESC_END

    CSO_PROG_ARGS_SHORT_DESC_BEGIN
    "accuname1 [accuname2 [..]] outaccuname"
    CSO_PROG_ARGS_SHORT_DESC_END

    CSO_PROG_ARGS_LONG_DESC_BEGIN
    "<cyan><b>accuname1</b></cyan>                  File name of the PMF accumulator.\n"
    "<cyan><b>outaccuname</b></cyan>                File name of the resulting PMF accumulator."
    CSO_PROG_ARGS_LONG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // options ------------------------------
    CSO_OPT(bool,Linear)
    CSO_OPT(bool,CommonOnly)
    CSO_OPT(int,Limit)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of options ---------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Linear,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "linear",                      /* long option name */
                NULL,                           /* parameter name */
                "Combine PMF accumulators by step-by-step algorithm.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                CommonOnly,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "common",                      /* long option name */
                NULL,                           /* parameter name */
                "Combine only common sections.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "LIMIT",                           /* parameter name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
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
