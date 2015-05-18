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

#include "STMTrajOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSTMTrajOptions::CSTMTrajOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

/*
 check validity of specified options
*/

int CSTMTrajOptions::CheckOptions(void)
{
// limit has to be grater than 0
    if(GetOptSnapshot() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: snapshot has to be grater than zero, but %d specified\n", (const char*)GetProgramName(),GetOptSnapshot());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

/*
 process special options (Help, Version) before arguments will be processed
*/

int CSTMTrajOptions::FinalizeOptions(void)
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
