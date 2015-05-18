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

#include "BMIntOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CBMIntOptions::CBMIntOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

int CBMIntOptions::CheckOptions(void)
{
// check if options are in the proper ranges
// CSO_OPT(int,SkipLines)
// CSO_OPT(int,AnalLines)
// CSO_OPT(int,PadLines)

    if(GetOptSkipLines() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: number of skipped lines has to be grater or equal to zero, but %d specified\n", (const char*)GetProgramName(),GetOptSkipLines());
        IsError = true;
    }
    if(GetOptAnalLines() < -1) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: number of analysed lines has to be grater or equal to minus one, but %d specified\n", (const char*)GetProgramName(),GetOptAnalLines());
        IsError = true;
    }
    if(GetOptPadLines() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: number of padding lines has to be grater or equal to zero, but %d specified\n", (const char*)GetProgramName(),GetOptPadLines());
        IsError = true;
    }

// CSO_OPT(CSmallString,IXFormat)
// CSO_OPT(CSmallString,IYFormat)
// CSO_OPT(CSmallString,ISFormat)
// CSO_OPT(CSmallString,OIFormat)
// CSO_OPT(CSmallString,OEFormat)

    if(IsError == true) return(SO_OPTS_ERROR);

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CBMIntOptions::FinalizeOptions(void)
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
