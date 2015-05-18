#ifndef STMTrajOptionsH
#define STMTrajOptionsH
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

class CSTMTrajOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CSTMTrajOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "stm-trajectory"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "It extracts a single snapshot from the STM trajectory."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    PMFLIB_VERSION
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,STMTrajName)
    CSO_ARG(CSmallString,STMResultName)
    // options ------------------------------
    CSO_OPT(int,Snapshot)
    CSO_OPT(bool,Verbose)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Help)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                STMTrajName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "trjname",                        /* parametr name */
                "Name of file containing the STM trajectory. If the name is '-' then the accumulator is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                STMResultName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "accuname",                        /* parametr name */
                "Name of file where a snapsphot (STM results) will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Snapshot,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "snapshot",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Index of extracted snapshot (counted from 0).")   /* option description */
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
