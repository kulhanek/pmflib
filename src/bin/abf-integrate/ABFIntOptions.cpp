// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include "ABFIntOptions.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntOptions::CABFIntOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

/*
 check validity of specified options
*/

int CABFIntOptions::CheckOptions(void)
{
    if(GetOptLimit() < 0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: sampling limit has to be grater or equal to zero, but %d is specified\n", (const char*)GetProgramName(),GetOptLimit());
        IsError = true;
    }

    if(GetOptNCorr() < 1.0) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: ncorr has to be grater or equal to one, but %f is specified\n", (const char*)GetProgramName(),GetOptNCorr());
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
        fprintf(stderr,"%s: method must be either rfd or rbf, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptMethod());
        IsError = true;
    }

    if( GetOptMethod() == "rfd" ){
        if((GetOptFDPoints() != 3) && (GetOptFDPoints() != 4)) {
            if(IsError == false) fprintf(stderr,"\n");
            fprintf(stderr,"%s: RFD number of points has to be either three or four, but %d is specified\n", (const char*)GetProgramName(),GetOptFDPoints());
            IsError = true;
        }
    }

    if( GetOptRFac() <= 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: rfac has to be greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptRFac());
        IsError = true;
    }

    if( GetOptWFac() <= 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: wfac has to be greater than zero, but %f is provided\n", (const char*)GetProgramName(),GetOptWFac());
        IsError = true;
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

    if( GetOptOverhang() <= 0 ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: overhang has to be greater than or equal zero, but %d is provided\n", (const char*)GetProgramName(),GetOptOverhang());
        IsError = true;
    }

    if((GetOptOutputFormat() != "plain") &&
            (GetOptOutputFormat() != "gnuplot") &&
            (GetOptOutputFormat() != "fes")) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: output FES format has to be either plain, gnuplot, or fes, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptOutputFormat());
        IsError = true;
    }

    if((GetOptOutputFormat() != "plain") &&
            (GetOptOutputFormat() != "gnuplot") &&
            (GetOptOutputFormat() != "fes")) {
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: output FES format has to be either plain, gnuplot, or fes, but %s is specified\n",
                (const char*)GetProgramName(),(const char*)GetOptOutputFormat());
        IsError = true;
    }

    if( IsOptSigmaF2Set() && (GetOptMethod() != "gpr") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: sigmaf2 can be set only for GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptRCondSet() && (GetOptMethod() != "rbf") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: rcond can be set only for RBF method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptRFacSet() && (GetOptMethod() != "rbf") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: rfac can be set only for RBF method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptOverhangSet() && (GetOptMethod() != "rbf") ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: overhang can be set only for RBF method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if( IsOptWFacSet() && ( (GetOptMethod() != "rbf") && (GetOptMethod() != "gpr")) ){
        if(IsError == false) fprintf(stderr,"\n");
        fprintf(stderr,"%s: wfac can be set only for RBF or GPR method\n",
                (const char*)GetProgramName());
        IsError = true;
    }

    if(IsError == true) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

/*
 process special options (Help, Version) before arguments will be processed
*/

int CABFIntOptions::FinalizeOptions(void)
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
