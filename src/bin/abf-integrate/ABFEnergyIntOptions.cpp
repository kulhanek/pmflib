// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include "ABFEnergyIntOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFEnergyIntOptions::CABFEnergyIntOptions(void)
{
    SetShowMiniUsage(true);
    SetAllowProgArgs(true);
}

//------------------------------------------------------------------------------

int CABFEnergyIntOptions::CheckOptions(void)
{
    if(GetOptLimit() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: sampling limit has to be grater or equal to zero, but %d is specified\n", (const char*)GetProgramName(),GetOptLimit());
        IsError = true;
    }

    if((GetOptEnergyLimit() <= 0) && (GetOptEnergyLimit() != -1)) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: energy limit has to be grater than zero or equal to -1, but %f is specified\n", (const char*)GetProgramName(),GetOptEnergyLimit());
        IsError = true;
    }

    if( (GetOptMethod() != "rfd") &&
        (GetOptMethod() != "rbf") &&
        (GetOptMethod() != "gpr") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: method must be either rfd, rbf, or gpr, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptMethod());
        IsError = true;
    }

    if( (GetOptEcutMethod() != "rfd") &&
        (GetOptEcutMethod() != "rbf") &&
        (GetOptEcutMethod() != "gpr") ) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: emethod must be either rfd, rbf, or gpr, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptEcutMethod());
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

    if( GetOptMethod() == "rfd" ){
        if((GetOptFDPoints() != 3) && (GetOptFDPoints() != 4)) {
            if(IsError == false) fprintf(stderr,"\n");
            fprintf(stderr,"%s: RFD number of points has to be either three or four, but %d is specified\n", (const char*)GetProgramName(),GetOptFDPoints());
            IsError = true;
        }
    }

    if( GetOptSigmaF2() <= 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: sigmaf2 has to be greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptSigmaF2());
        IsError = true;
    }

    if( (GetOptRCond() != -1) && (GetOptRCond() < 0) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: rcond has to be either -1 or greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptRCond());
        IsError = true;
    }

    if( GetOptOverhang() < 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: overhang has to be greater than or equal zero, but %d is provided\n", (const char*)GetProgramName(),GetOptOverhang());
        IsError = true;
    }

    if((GetOptOutputFormat() != "plain") &&
       (GetOptOutputFormat() != "gnuplot")) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: output FES format has to be either plain or gnuplot, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptOutputFormat());
        IsError = true;
    }

    if( IsOptSigmaF2Set() && ( (GetOptMethod() != "gpr") && (GetOptEcutMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --sigmaf2 can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptNoEnergySet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --noenergy can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRNumDiffSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --numdiff can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRUseInvSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --useinv can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRCalcLogPLSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --calclogpl can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRIncludeZPESet() && ((GetOptMethod() != "gpr") || (IsOptGlobalMinSet() == false)) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --inczpe can be set only for GPR method and alongside with --globalmin\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRKernelSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --kernel can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRNoFastErrorSet() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --nofasterror can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptGPRNoFastErrorSet() && (IsOptWithErrorSet() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --nofasterror can be set only alongside with --witherror\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptKeepCVsSet() && (IsOptReducedFESSet() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --keepcvs can be set only alongside with --reducedfes\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptReducedFESSet() && (IsOptKeepCVsSet() == false) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --reducedfes can be set only alongside with --keepcvs\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptRFacSet() && (GetOptMethod() != "rbf") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --rfac can be set only for RBF method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptOverhangSet() && (GetOptMethod() != "rbf") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --overhang can be set only for RBF method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptWFacSet() && ( (GetOptMethod() != "rbf") && (GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --wfac can be set only for RBF or GPR method\n",
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

    if( IsOptMFInfoSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --mfinfo can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptMFMaxZScoreSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --maxzscore can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptMFZTestPassesSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --mfzpasses can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptGlueingFactorSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --glueing can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptGlueHolesSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --glueholes can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptGlobalMinSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --globalmin can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptUseRealGlobalMinSet() && ((GetOptMethod() != "rbf")&&(GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --userealglbmin can be combined only with RBF or GPR\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if( IsOptUseOldRFDSet() && (GetOptMethod() != "rfd") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: --oldrfd can be combined only with RFD\n",
                (const char*)GetProgramName());
        IsError = true;
    }
    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CABFEnergyIntOptions::CheckArguments(void)
{
    if( (GetNumberOfProgArgs() != 2) && (GetNumberOfProgArgs() != 3) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: two or three arguments are expected, but %d is provided\n",
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

int CABFEnergyIntOptions::FinalizeOptions(void)
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
