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
#include "ACCUCombine.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

MAIN_ENTRY(CACCUCombine)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CACCUCombine::CACCUCombine(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CACCUCombine::Init(int argc,char* argv[])
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
    vout << "# accu-combine (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    for(int i=0; i < Options.GetNumberOfProgArgs() - 1; i++){
        vout << format("# Input PMF accumulator #%03 : %s")%(i+1)%Options.GetProgArg(i) << endl;
    }
    vout << format("# Output PMF accumulator         : %s")%Options.GetProgArg(Options.GetNumberOfProgArgs()-1) << endl;

    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CACCUCombine::Run(void)
{
// load accumulators
    int state = 1;
    vout << endl;
    vout << format("%02d:Loading input ABF accumulators ...")%state << endl;
    try {
        for(int i=0; i < Options.GetNumberOfProgArgs() - 1; i++){
            vout << format("  >> %03d %s")%state%string(Options.GetProgArg(i)) << endl;
            CPMFAccumulatorPtr inaccu = CPMFAccumulatorPtr(new CPMFAccumulator);
            inaccu->Load(Options.GetProgArg(i));
            InAccus.push_back(inaccu);
        }
    } catch(...) {
        ES_ERROR("unable to load the input PMF accumulator");
        return(false);
    }

    state++;
    vout << endl;
    vout << format("%02d:Combining PMF accumulators ...")%state << endl;
    OutAccu = CPMFAccumulatorPtr(new CPMFAccumulator);


    // FIXME

    state++;
    vout << endl;
    vout << format("%02d:Saving the resulting PMF accumulator: %s")%state%Options.GetProgArg(Options.GetNumberOfProgArgs()-1) << endl;

// save final accumulator
    try {
        OutAccu->Save(Options.GetProgArg(Options.GetNumberOfProgArgs()-1));
    } catch(...) {
        ES_ERROR("unable to save final PMF accumulator file");
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CACCUCombine::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-combine terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

