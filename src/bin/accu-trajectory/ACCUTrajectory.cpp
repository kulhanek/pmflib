// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "ACCUTrajectory.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

MAIN_ENTRY(CACCUTrajectory)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CACCUTrajectory::CACCUTrajectory(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CACCUTrajectory::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CACCUIntOpts
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
    vout << "# accu-trajectory (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgTrjName() != "-") {
        vout << "# PMF accumulator trajectory file (in) : " << Options.GetArgTrjName() << endl;
    } else {
        vout << "# PMF accumulator trajectory file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgOutName() != "-") {
        vout << "# Output PMF accumulator file (out)    : " << Options.GetArgOutName() << endl;
    } else {
        vout << "# Output PMF accumulator file (out)    : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
        vout << "# Snapshot                             : " << Options.GetOptSnapshot() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgTrjName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgOutName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CACCUTrajectory::Run(void)
{
    int state = 1;
    vout << endl;
    vout << format("%02d:Loading a snapshot from the trajectory of the PMF accumulators ...")%state << endl;

// read accumulator
    try {
        Accu.LoadSnapshot(InputFile,Options.GetOptSnapshot());
    } catch(...) {
        ES_ERROR("unable to read snapshot from input");
        return(false);
    }

    state++;
    vout << endl;
    vout << format("%02d:Saving the snapshot as the PMF accumulator ...")%state << endl;

// now write snapshot
    try {
        Accu.Save(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save selected snapshot");
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CACCUTrajectory::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-trajectory terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

