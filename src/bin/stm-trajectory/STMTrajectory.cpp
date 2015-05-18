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
#include "STMTrajectory.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CSTMTrajectory)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSTMTrajectory::CSTMTrajectory(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CSTMTrajectory::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CSTMIntOpts
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
    vout << "# stm-trajectory (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgSTMTrajName() != "-") {
        vout << "# STM accumulator file (in) : " << Options.GetArgSTMTrajName() << endl;
    } else {
        vout << "# STM accumulator file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgSTMResultName() != "-") {
        vout << "# Snapshot file (out)       : " << Options.GetArgSTMResultName() << endl;
    } else {
        vout << "# Snapshot file (out)       : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# Snapshot                  : " << Options.GetOptSnapshot() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgSTMTrajName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgSTMResultName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSTMTrajectory::Run(void)
{
// read first line
    CSmallString   title;
    if(title.ReadLineFromFile(InputFile) == false) {
        ES_ERROR("unable to read first line from input stream");
        return(false);
    }

    if( title.FindSubString("# STMTRAJ") != 0 ) {
        ES_ERROR("input file is not STM trajectory file");
        return(false);
    }

// move to requested position
    int  position = -1;
    bool found = false;
    CSmallString   line;
    while( (feof(InputFile) == 0) && (line.ReadLineFromFile(InputFile,true,true) == true) ) {
        // is it segment?
        if( line.FindSubString("# STMSNAP") == 0 ){
            position++;
            if( position == Options.GetOptSnapshot() ){
                found = true;
                break;
            }
            continue;
        }
        // is it comment?
        if( (line.GetLength() > 0) && (line[0] == '#') ){
            // add it to header
            Header << line << endl;
        }
    }

    if( found == false ) {
        CSmallString error;
        error << "trajectory contains less snapshots than " << Options.GetOptSnapshot();
        ES_ERROR(error);
        return(false);
    }

    // write header to output
    fputs(Header.str().c_str(),OutputFile);
    fprintf(OutputFile,"# STMSNAP %d\n",position);

    // copy snapshot to output
    while( (feof(InputFile) == 0) && (line.ReadLineFromFile(InputFile,true,true) == true) ) {
        if( line.GetLength() == 0 ) break; // end of segment
        fprintf(OutputFile,"%s\n",line.GetBuffer());
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMTrajectory::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# stm-trajectory terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

