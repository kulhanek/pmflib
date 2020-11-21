#ifndef ABPEneOptionsH
#define ABPEneOptionsH
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

class CABPEneOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABPEneOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abp-energy"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "It calculate the free energy surface estimate from the metadynamics restart file."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,Input)
    CSO_ARG(CSmallString,Output)
    // options ------------------------------
    CSO_OPT(CSmallString,Mode)
    CSO_OPT(int,Limit)
    CSO_OPT(bool,PrintAll)
    CSO_OPT(bool,UnsampledAsMaxE)
    CSO_OPT(double,MaxEnergy)
    CSO_OPT(bool,IncludeBinStat)
    CSO_OPT(double,Offset)
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
                "input",                        /* parameter name */
                "Name of file containing the ABP restart file. If the name is '-' then the file is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Output,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "output",                        /* parameter name */
                "Name of file where the resulting free energy surface will be printed. If the name is '-' then the output will be written to the standard output.")   /* argument description */
// description of options ---------------------------------------------------
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Mode,                        /* option name */
                "rl",                          /* default value */
                false,                          /* is option mandatory */
                'm',                           /* short option name */
                "mode",                      /* long option name */
                "NAME",                           /* parameter name */
                "Energy calculation mode: mollified (use raw mollified FES), rl (Lucy–Richardson deconvolution)")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                Limit,                        /* option name */
                100,                          /* default value */
                false,                          /* is option mandatory */
                'l',                           /* short option name */
                "limit",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Only bins containing more samples than NUMBER are considered as properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                Offset,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                'o',                           /* short option name */
                "offset",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "Specify an integration constant.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                PrintAll,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "printall",                      /* long option name */
                NULL,                           /* parameter name */
                "Print results for all bins even if not properly sampled.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                UnsampledAsMaxE,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "unsampledasmax",                      /* long option name */
                NULL,                           /* parameter name */
                "Set energy values in unsampled region to maximum energy from sampled region or to value provided by --maxenergy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(double,                           /* option type */
                MaxEnergy,                        /* option name */
                0.0,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "maxenergy",                      /* long option name */
                "NUMBER",                           /* parameter name */
                "If set, this is the energy used of unsampled regions.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                IncludeBinStat,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "includebinstat",                      /* long option name */
                NULL,                           /* parameter name */
                "Include bin statuses (1=sampled, 0=unsampled, -1=glued) into resulting FES.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OutputFormat,                        /* option name */
                "gnuplot",                          /* default value */
                false,                          /* is option mandatory */
                0,                           /* short option name */
                "output",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT, which will be used to print free energy surface. Supported formats are: plain, gnuplot, fes.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                NoHeader,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "noheader",                      /* long option name */
                NULL,                           /* parameter name */
                "Do not print header to output file.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                IXFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fx",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT, which will be used to print values of collective variables.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OEFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fe",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT, which will be used to print values of free energy.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                OSFormat,                        /* option name */
                "%15.7le",                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "fs",                      /* long option name */
                "FORMAT",                           /* parameter name */
                "Output FORMAT, which will be used to print values of standard deviations.")   /* option description */

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
