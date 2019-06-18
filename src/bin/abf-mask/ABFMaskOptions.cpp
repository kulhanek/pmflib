// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "ABFMaskOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFMaskOptions::CABFMaskOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

/*
 check validity of specified options
*/

int CABFMaskOptions::CheckOptions(void)
{
// limit has to be grater than 0
    if(GetOptLimit() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: limit option has to be grater or equal to zero, but %d is specified\n", (const char*)GetProgramName(),GetOptLimit());
        IsError = true;
    }

    if((GetOptFDPoints() != 3) && (GetOptFDPoints() != 4)) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: FD number of points has to be either three or four, but %d is specified\n", (const char*)GetProgramName(),GetOptFDPoints());
        IsError = true;
    }

    if( (GetOptOutputFormat() != "mask") &&
        (GetOptOutputFormat() != "gnuplot")    ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: output FES format has to be either mask or gnuplot, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptOutputFormat());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

/*
 process special options (Help, Version) before arguments will be processed
*/

int CABFMaskOptions::FinalizeOptions(void)
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
