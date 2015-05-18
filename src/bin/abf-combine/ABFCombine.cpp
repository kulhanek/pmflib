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
#include "ABFCombine.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFCombine)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFCombine::CABFCombine(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFCombine::Init(int argc,char* argv[])
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
    vout << "# abf-combine (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;
    if(Options.GetArgABFAccuName1() != "-") {
        vout << "# ABF accumulator file 1 (in) : " << Options.GetArgABFAccuName1() << endl;
    } else {
        vout << "# ABF accumulator file 1 (in) : - (standard input)" << endl;
    }

    if(Options.GetArgABFAccuName2() != "-") {
        vout << "# ABF accumulator file 2 (in) : " << Options.GetArgABFAccuName2() << endl;
    } else {
        vout << "# ABF accumulator file 2 (in) : - (standard input)" << endl;
    }

    if(Options.GetArgOutputName() != "-") {
        vout << "# ABF accumulator file (out)  : " << Options.GetArgOutputName() << endl;
    } else {
        vout << "# ABF accumulator file (out)  : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# Operation                   : " << Options.GetOptOperation() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile1.Open(Options.GetArgABFAccuName1(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( InputFile2.Open(Options.GetArgABFAccuName2(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgOutputName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFCombine::Run(void)
{
// load accumulators
    try {
        Accu1.Load(InputFile1);
    } catch(...) {
        ES_ERROR("unable to load the first ABF accumulator file");
        return(false);
    }

    try {
        Accu2.Load(InputFile2);
    } catch(...) {
        ES_ERROR("unable to load the second ABF accumulator file");
        return(false);
    }

    bool op_ok = false;

    if(Options.GetOptOperation() == "add") {
        Accu1.AddABFAccumulator(&Accu2);
        op_ok = true;
    }

    if(Options.GetOptOperation() == "sub") {
        Accu1.SubABFAccumulator(&Accu2);
        op_ok = true;
    }

    if(op_ok == false) {
        CSmallString error;
        error << "operation " << Options.GetOptOperation() << " is not implemented";
        ES_ERROR("unable to save final ABF accumulator file");
        return(false);
    }

// save final accumulator
    try {
        Accu1.Save(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save final ABF accumulator file");
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFCombine::Finalize(void)
{
    // close files if they are own by program
    InputFile1.Close();
    InputFile2.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-combine terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

