#ifndef ABFAdmOptionsH
#define ABFAdmOptionsH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
// ===============================================================================

#include <SimpleOptions.hpp>
#include <PMFMainHeader.hpp>

//------------------------------------------------------------------------------

class CABFAdmOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CABFAdmOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "abf-admin"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "The utility that controls the behaviour and state of the server implementing the multiple walkers extension of the ABF method."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    LibBuildVersion_PMF
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,Command)
    // options ------------------------------
    CSO_OPT(CSmallString,ServerKey)
    CSO_OPT(CSmallString,Password)
    CSO_OPT(bool,Force)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Command,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "abf://server[:port]/command",                        /* parametr name */
                "It provides the action specification. The server is either the DNS name or IP address of the server or word 'serverkey'. In the later case, the information about the server is read from the server key file. The port number, on which the server is listen, may be optionally provided. Finally, the command is administration task, which can be one of the following:\n"
                "info           = print information about registered clients\n"
                "flush          = save the accumulated ABF data on the server side\n"
                "get?file=NAME  = get the accumulated ABF data and saves them locally to the NAME file (the default name is _abfserver.rst)\n"
                "shutdown       = stop the server execution (use --force to skip protection by passphrase)\n"
                "errors         = print errors from the server stack\n")   /* argument description */
// description of options -----------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                ServerKey,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "serverkey",                      /* long option name */
                "FILE",                           /* parametr name */
                "Name of file containing the server key. The server key contains the server name, port, and password. This option is mutually exclusive with 'password' option.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Password,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "password",                        /* long option name */
                "FILE",                           /* parametr name */
                "Name of file containing the server magic word. If the pasword is not provided via this option or via the server key then it is read interactively from the keyboard. This option is mutually exclusive with 'serverkey' option.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Force,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'f',                           /* short option name */
                "force",                      /* long option name */
                NULL,                           /* parametr name */
                "Skip protection of server shutdown by passphrase.")   /* option description */
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
    virtual int CheckArguments(void);
};

//------------------------------------------------------------------------------

#endif
