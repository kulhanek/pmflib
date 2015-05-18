#ifndef BMIntOptionsH
#define BMIntOptionsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
//                       Martin Petrek, petrek@chemi.muni.cz
//    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CBMIntOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CBMIntOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "con-integrate"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "Numerically integrate given data by Simpson method. Program reads three numbers from input file (x coordinate, derivative, and standard deviation of derivative). As output program prints copy of input data plus integrated values and its standard error."
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
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_OPT(int,SkipLines)
    CSO_OPT(int,AnalLines)
    CSO_OPT(int,PadLines)
    CSO_OPT(double,Offset)
    CSO_OPT(bool,NoSigma)
    CSO_OPT(bool,NoHeader)
    CSO_OPT(CSmallString,IXFormat)
    CSO_OPT(CSmallString,IYFormat)
    CSO_OPT(CSmallString,ISFormat)
    CSO_OPT(CSmallString,OIFormat)
    CSO_OPT(CSmallString,OEFormat)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Input,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "input",                        /* parametr name */
                "input data file or - for data taken from standard input")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Output,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "output",                        /* parametr name */
                "output data file or - for data written to standard output")   /* argument description */
// description of options ---------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Offset,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "offset",                      /* long option name */
                "REAL",                           /* parametr name */
                "specifies integration constant")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                SkipLines,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "skip",                      /* long option name */
                "LINES",                           /* parametr name */
                "number of lines skipped from the beginning of input file")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                AnalLines,                        /* option name */
                -1,                          /* default value */
                false,                          /* is option mandatory */
                'a',                           /* short option name */
                "anal",                      /* long option name */
                "LINES",                           /* parametr name */
                "number of LINES from input file to be analyzed")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                PadLines,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "pad",                      /* long option name */
                "LINES",                           /* parametr name */
                "number of padding LINES between used records from input file")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoSigma,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "nosigma",                      /* long option name */
                NULL,                           /* parametr name */
                "standard deviation of derivative is not present in input data")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoHeader,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "noheader",                      /* long option name */
                NULL,                           /* parametr name */
                "do not print header to output")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IXFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fx",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "output FORMAT of echoed x values")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IYFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fy",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "output FORMAT of echoed derivative values")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                ISFormat,                        /* option name */
                "%14.6le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fs",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "output FORMAT of echoed sigma values")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OIFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fi",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "output FORMAT of integrated values")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OEFormat,                        /* option name */
                "%14.6le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parametr name */
                "output FORMAT of integrated error values")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Verbose,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'v',                           /* short option name */
                "verbose",                      /* long option name */
                NULL,                           /* parametr name */
                "increase output verbosity")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Version,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "version",                      /* long option name */
                NULL,                           /* parametr name */
                "output version information and exit")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Help,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'h',                           /* short option name */
                "help",                      /* long option name */
                NULL,                           /* parametr name */
                "display this help and exit")   /* option description */
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int FinalizeOptions(void);
    virtual int CheckOptions(void);
};

//------------------------------------------------------------------------------

#endif
