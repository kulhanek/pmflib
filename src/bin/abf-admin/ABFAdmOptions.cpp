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

#include "ABFAdmOptions.hpp"
#include <ActionRequest.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFAdmOptions::CABFAdmOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

int CABFAdmOptions::CheckOptions(void)
{
    if(IsOptPasswordSet() && IsOptServerKeySet()) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --password and --serverkey options are mutually exclusive\n", (char*)GetProgramName());
        IsError = true;
        return(SO_OPTS_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CABFAdmOptions::FinalizeOptions(void)
{
    bool ret_opt = false;

    if(GetOptHelp() == true) {
        PrintUsage();
        ret_opt = true;
    }

    if(GetOptVersion() == true) {
        PrintVersion();
        ret_opt = true;
    }

    if(ret_opt == true) {
        printf("\n");
        return(SO_EXIT);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CABFAdmOptions::CheckArguments(void)
{
    CActionRequest server_desc;
    server_desc.SetProtocolName("abf");

    try {
        server_desc.SetQualifiedName(GetArgCommand());
    } catch(...) {
        if(IsVerbose()) {
            if(IsError == false) fprintf(stderr,"\n");
            fprintf(stderr,"%s: %s\n", (const char*)GetProgramName(), (const char*)ErrorSystem.GetLastError());
            IsError = true;
        }
        return(SO_OPTS_ERROR);
    }

    if(server_desc.GetAction() == NULL) {
        if(IsVerbose()) {
            if(IsError == false) fprintf(stderr,"\n");
            fprintf(stderr,"%s: no action is specified\n", (const char*)GetProgramName());
            IsError = true;
        }
        return(SO_OPTS_ERROR);
    }

    if((server_desc.GetAction() != "info") &&
            (server_desc.GetAction() != "shutdown") &&
            (server_desc.GetAction() != "errors") &&
            (server_desc.GetAction() != "get") &&
            (server_desc.GetAction() != "flush")) {
        if(IsVerbose()) {
            if(IsError == false) fprintf(stderr,"\n");
            fprintf(stderr,"%s: specified action '%s' is not supported\n", (const char*)GetProgramName(), (const char*)server_desc.GetAction());
            IsError = true;
        }
        return(SO_OPTS_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
