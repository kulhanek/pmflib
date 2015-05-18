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

#include <errno.h>
#include "ABFTrajectory.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFTrajectory)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFTrajectory::CABFTrajectory(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFTrajectory::Init(int argc,char* argv[])
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

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-trajectory (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFTrajName() != "-") {
        vout << "# ABF accumulator file (in) : " << Options.GetArgABFTrajName() << endl;
    } else {
        vout << "# ABF accumulator file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgABFAccuName() != "-") {
        vout << "# Snapshot file (out)       : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# Snapshot file (out)       : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# Snapshot                  : " << Options.GetOptSnapshot() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFTrajName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgABFAccuName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFTrajectory::Run(void)
{
// read first line
    CSmallString   title;
    if(title.ReadLineFromFile(InputFile) == false) {
        ES_ERROR("unable to read first line from input stream");
        return(false);
    }

    if(title != "# ABFTRAJ\n") {
        ES_ERROR("input file is not ABF trajectory file");
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
        if( time_cards.FindSubString("# ABFSNAP") != 0 ) {
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
            if( time_cards.FindSubString("# ABFSNAP") == 0 ) break;
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

void CABFTrajectory::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-trajectory terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

