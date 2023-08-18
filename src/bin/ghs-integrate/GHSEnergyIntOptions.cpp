// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Martin Petrek, petrek@chemi.muni.cz
//                       Petr Kulhanek, kulhanek@enzim.hu
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

#include "GHSEnergyIntOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGHSEnergyIntOptions::CGHSEnergyIntOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

int CGHSEnergyIntOptions::CheckOptions(void)
{
    if(GetOptLimit() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: sampling limit has to be grater or equal to zero, but %d is specified\n", (const char*)GetProgramName(),GetOptLimit());
        IsError = true;
    }

    if( (GetOptLAMethod() != "svd") &&
        (GetOptLAMethod() != "svd2") &&
        (GetOptLAMethod() != "qr") &&
        (GetOptLAMethod() != "lu") &&
        (GetOptLAMethod() != "ll") &&
        (GetOptLAMethod() != "default") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: linear algebra method must be either svd, svd2, qr, lu, ll, or default, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptLAMethod());
        IsError = true;
    }

    if( (GetOptRCond() != -1) && (GetOptRCond() < 0) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: rcond has to be either -1 or greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptRCond());
        IsError = true;
    }

    if((GetOptOutputFormat() != "plain") &&
       (GetOptOutputFormat() != "gnuplot")) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: output FES format has to be either plain or gnuplot, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptOutputFormat());
        IsError = true;
    }


    if( IsOptGPRNoFastErrorSet() && (IsOptWithErrorSet() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --nofasterror can be set only alongside with --witherror\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptMaxEnergySet() && (GetOptUnsampledAsMaxE() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --maxenergy without --unsampledasmax does not have any effect\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptLoadHyprmsSet() && IsOptSigmaF2Set() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms is mutually exclusive with --sigmaf2\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptLoadHyprmsSet() && IsOptWFacSet() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms is mutually exclusive with --wfac\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptLoadHyprmsSet() && IsOptSigmaN2Set() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms is mutually exclusive with --sigman2\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CGHSEnergyIntOptions::CheckArguments(void)
{
    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

/*
 process special options (Help, Version) before arguments will be processed
*/

int CGHSEnergyIntOptions::FinalizeOptions(void)
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
