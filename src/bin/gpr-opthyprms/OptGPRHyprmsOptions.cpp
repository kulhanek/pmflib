// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "OptGPRHyprmsOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

COptGPRHyprmsOptions::COptGPRHyprmsOptions(void)
{
    SetShowMiniUsage(true);
    SetAllowProgArgs(true);
}

//------------------------------------------------------------------------------

int COptGPRHyprmsOptions::CheckOptions(void)
{
    if( (GetOptTarget() != "logml") &&
        (GetOptTarget() != "logpl") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: optimized target must be either logml or logpl, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptTarget());
        IsError = true;
    }

    if( (GetOptLAMethod() != "svd") &&
        (GetOptLAMethod() != "svd2") &&
        (GetOptLAMethod() != "ll") &&
        (GetOptLAMethod() != "lu") &&
        (GetOptLAMethod() != "default") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: linear algebra method must be either svd, svd2, qr, lu, or default, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptLAMethod());
        IsError = true;
    }

    if( GetOptNOptSteps() < 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --noptsteps: negative number of steps (%d) is not allowed\n", (const char*)GetProgramName(),GetOptNOptSteps());
        IsError = true;
    }

    if( GetOptSigmaF2() <= 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --sigmaf2 has to be greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptSigmaF2());
        IsError = true;
    }

    if( (GetOptRCond() != -1) && (GetOptRCond() < 0) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --rcond has to be either -1 or greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptRCond());
        IsError = true;
    }

    if( IsOptLoadHyprmsSet() && IsOptSigmaF2Set() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms is mutually exclusive with --sigmaf2\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptLoadHyprmsSet() && IsOptNCorrSet() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms is mutually exclusive with --ncorr\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptLoadHyprmsSet() && IsOptWFacSet() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms is mutually exclusive with --wfac\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptTestSet() && IsOptSPTypeSet() ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --sptype is mutually exclusive with --test\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRIncludeZPESet() && (IsOptGlobalMinSet() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --inczpe can be set only alongside with --globalmin\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int COptGPRHyprmsOptions::CheckArguments(void)
{
    if( GetNumberOfProgArgs() < 2 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: at least two arguments are expected, but %d is provided\n",
                (const char*)GetProgramName(),GetNumberOfProgArgs());
        IsError = true;
    }
    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

/*
 process special options (Help, Version) before arguments will be processed
*/

int COptGPRHyprmsOptions::FinalizeOptions(void)
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
