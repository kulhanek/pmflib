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

#include "ABFCombOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFCombOptions::CABFCombOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

/*
 check validity of specified options
*/

int CABFCombOptions::CheckOptions(void)
{
    if((GetOptOperation() != "add") && (GetOptOperation() != "sub")) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: 'add' and 'sub' operation are supported only, but '%s' is provided\n", (const char*)GetProgramName(),(const char*)GetOptOperation());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

/*
 process special options (Help, Version) before arguments will be processed
*/

int CABFCombOptions::FinalizeOptions(void)
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

int CABFCombOptions::CheckArguments(void)
{
    if((GetArgABFAccuName1() == "-") && (GetArgABFAccuName1() == GetArgABFAccuName2())) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: standard input can be used only for one input file\n", (const char*)GetProgramName());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
