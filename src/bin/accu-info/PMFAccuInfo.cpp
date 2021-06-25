// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include "PMFAccuInfo.hpp"
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CPMFAccuInfo)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccuInfo::CPMFAccuInfo(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFAccuInfo::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << debug;
    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-info (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetProgArg(0) != "-") {
        vout << "# PMF accumulator (in)  : " << Options.GetProgArg(0) << endl;
    } else {
        vout << "# PMF accumulator (in)  : - (standard input)" << endl;
    }
    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetProgArg(0),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CPMFAccuInfo::Run(void)
{
// load accumulator
    vout << "#" << endl;
    vout << "# 1) Loading the PMF accumulator: " << Options.GetProgArg(0) << endl;
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input PMF accumulator file");
        return(false);
    }
    vout << "#   Done" << endl;

// perform action
    if( Options.GetNumberOfProgArgs() == 1 ) {
        PrintInfo();
// -----------------------------------------------
    } else  if( Options.GetProgArg(1)  == "info" ){
        PrintInfo();
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "list-sections" ){
        ListSections();
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-section" ){
        GetSection(Options.GetProgArg(2));
// -----------------------------------------------
    } else if( (Options.GetProgArg(1) == "NSAMPLES") || (Options.GetProgArg(1) == "nsamples") ){
        GetSection("NSAMPLES");
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unrecognized action: " << Options.GetProgArg(1);
        ES_ERROR(error);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuInfo::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << "#" << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-info terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        vout << high;
        ErrorSystem.PrintErrors(vout);
        vout << endl;
        vout << debug;
    } else {
        vout << "#" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuInfo::PrintInfo(void)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# 2) Info about the PMF accumulator."<< endl;

    vout << high;
    // print CVS info
    Accumulator.PrintAccuInfo(vout);
    Accumulator.PrintCVSInfo(vout);
    vout << endl;

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::ListSections(void)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# 2) List of data sections stored in the PMF accumulator."<< endl;

    vout << high;
    Accumulator.ListSections(vout);
    vout << endl;

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetSection(const CSmallString& name)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# Data section: " << name << endl;

    vout << high;


    vout << debug;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

