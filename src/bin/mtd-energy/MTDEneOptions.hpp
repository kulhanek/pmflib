#ifndef MTDEneOptionsH
#define MTDEneOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

class CMTDEneOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CMTDEneOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "mtd-energy"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "It calculate the free energy surface estimate from the metadynamics restart file."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    PMFLIB_VERSION
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,Input)
    CSO_ARG(CSmallString,Output)
    // options ------------------------------
    CSO_OPT(int,Time)
    CSO_OPT(double,Offset)
    CSO_OPT(int,Smooth)
    CSO_OPT(bool,PrintSD)
    CSO_OPT(CSmallString,OutputFormat)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,OEFormat)
    CSO_OPT(CSmallString,OSFormat)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Input,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "input",                        /* parametr name */
                "Name of file containing the metadynamics restart file. If the name is '-' then the file is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Output,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "output",                        /* parametr name */
                "Name of file where the resulting free energy surface will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Time,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "time",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Specifies metadynamics time that will be used for the energy calculation.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Offset,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "offset",                      /* long option name */
                "REAL",                           /* parametr name */
                "Determine position of global minima.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Smooth,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "smooth",                      /* long option name */
                "NUMBER",                           /* parametr name */
                "Calculate the energy surface average from time interval SMOOTH to TIME.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                PrintSD,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "sd",                      /* long option name */
                NULL,                           /* parametr name */
                "Print the standard deviation of the free energy when the smoothing is enabled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OutputFormat,                        /* option name */
                "gnuplot",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "output",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print free energy surface. Supported formats are: plain, gnuplot, fes.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoHeader,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "noheader",                      /* long option name */
                NULL,                           /* parametr name */
                "Do not print header to output file.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IXFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fx",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print values of collective variables.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OEFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print values of free energy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OSFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fs",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "Output FORMAT, which will be used to print values of standard deviations.")   /* option description */

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
