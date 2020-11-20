// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <errno.h>
#include "ABPTrajectory.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABPTrajectory)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABPTrajectory::CABPTrajectory(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPTrajectory::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABPIntOpts
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

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abp-trajectory (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABPTrajName() != "-") {
        vout << "# ABP accumulator file (in) : " << Options.GetArgABPTrajName() << endl;
    } else {
        vout << "# ABP accumulator file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgABPAccuName() != "-") {
        vout << "# Snapshot file (out)       : " << Options.GetArgABPAccuName() << endl;
    } else {
        vout << "# Snapshot file (out)       : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# Snapshot                  : " << Options.GetOptSnapshot() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABPTrajName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgABPAccuName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABPTrajectory::Run(void)
{
// read first line
    CSmallString   title;
    if(title.ReadLineFromFile(InputFile) == false) {
        ES_ERROR("unable to read first line from input stream");
        return(false);
    }

    if(title != "# ABPTRAJ\n") {
        ES_ERROR("input file is not ABP trajectory file");
        return(false);
    }

// read time cards
    CSmallString   time_cards;
    if(time_cards.ReadLineFromFile(InputFile) == false) {
        CSmallString error;
        error << "trajectory contains less snapshots than " << Options.GetOptSnapshot();
        ES_ERROR(error);
        return(false);
    }

// move to requested position
    int position = 0;
    while(feof(InputFile) == 0) {
        if( time_cards.FindSubString("# ABPSNAP") != 0 ) {
            CSmallString error;
            error << "trajectory contains less snapshots than " << Options.GetOptSnapshot();
            ES_ERROR(error);
            return(false);
        }
        position++;
        if(position == Options.GetOptSnapshot()) break;

        // find next time card
        while(feof(InputFile) == 0) {
            time_cards = NULL;
            if(time_cards.ReadLineFromFile(InputFile) == false) {
                CSmallString error;
                error << "trajectory contains less snapshots than " << Options.GetOptSnapshot();
                ES_ERROR(error);
                return(false);
            }
            if( time_cards.FindSubString("# ABPSNAP") == 0 ) break;
        }
    }

// read accumulator
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to read snapshot from input");
        return(false);
    }

// now write snapshot
    try {
        Accumulator.Save(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save selected snapshot");
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPTrajectory::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abp-trajectory terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

