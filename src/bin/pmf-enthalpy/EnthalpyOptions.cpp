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

#include "EnthalpyOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnthalpyOptions::CEnthalpyOptions(void)
{
    SetShowMiniUsage(true);
    SetAllowProgArgs(true);
}

//------------------------------------------------------------------------------

int CEnthalpyOptions::CheckOptions(void)
{
    if( (GetOptMethod() != "raw") &&
        (GetOptMethod() != "gpr") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: method must be either raw or gpr, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptMethod());
        IsError = true;
    }

    if( IsOptGPRKernelSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --kernel can be set only for GPR method\n",(const char*)GetProgramName());
        IsError = true;
    }

    if( (GetOptLAMethod() != "svd") &&
        (GetOptLAMethod() != "svd2") &&
        (GetOptLAMethod() != "lu") &&
        (GetOptLAMethod() != "ll") &&
        (GetOptLAMethod() != "default") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: linear algebra method must be either svd, svd2, lu, ll, or default, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptLAMethod());
        IsError = true;
    }

    if( IsOptGPRCalcLogPLSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --calclogpl can be set only for GPR method\n",(const char*)GetProgramName());
        IsError = true;
    }

// limit has to be grater than 0
    if(GetOptLimit() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: limit has to be grater or equal to 0, but %d specified\n", (const char*)GetProgramName(),GetOptLimit());
        IsError = true;
    }

    if( IsOptSigmaF2Set() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --sigmaf2 can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( GetOptSigmaF2() <= 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: sigmaf2 has to be greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptSigmaF2());
        IsError = true;
    }

    if( IsOptWFacSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --wfac can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptNCorrSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --ncorr can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptMaxEnergySet() && (GetOptUnsampledAsMaxE() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --maxenergy without --unsampledasmax does not have any effect\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptWithErrorSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --witherror can be combined only with GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptLoadHyprmsSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --loadhyprms can be combined only with GPR\n",
                (const char*)GetProgramName());
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

    if( IsOptMFInfoSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --mfinfo can be combined only with GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if((GetOptOutputFormat() != "plain") &&
       (GetOptOutputFormat() != "gnuplot")) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: output FES format has to be either plain or gnuplot, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptOutputFormat());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CEnthalpyOptions::FinalizeOptions(void)
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
