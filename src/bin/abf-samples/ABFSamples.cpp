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
#include "ABFSamples.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFSamples)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFSamples::CABFSamples(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFSamples::Init(int argc,char* argv[])
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
    vout << "# abf-samples (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accumulator file (in) : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# ABF accumulator file (in) : - (standard input)" <<  endl;
    }
    if(Options.GetArgOutputName() != "-") {
        vout << "# Sample file (out)         : " << Options.GetArgOutputName() << endl;
    } else {
        vout << "# Sample file (out)         : - (standard output)" <<  endl;
    }
    vout << "# ------------------------------------------------" <<  endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Limit                     : all bins will be printed" <<  endl;
    } else {
        vout << "# Limit                     : " << Options.GetOptLimit() << endl;
    }
    vout << "# ------------------------------------------------" <<  endl;
    vout << "# No GNUPlot delimiters     : " << bool_to_str(Options.GetOptNoGNUPlot()) << endl;
    vout << "# X format                  : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format                  : " << Options.GetOptOSFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" <<  endl;
    vout << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFAccuName(),"r") == false ){
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

bool CABFSamples::Run(void)
{
// load accumulator
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }

    Point.CreateVector(Accumulator.GetNumOfCVs());

// print samples
    if(PrintSamples() == false) {
        ES_ERROR("unable to print samples");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFSamples::PrintSamples(void)
{
    CSimpleVector<double>   pos;
    CSimpleVector<int>      ipos;

    pos.CreateVector(Accumulator.GetNumOfCVs());
    ipos.CreateVector(Accumulator.GetNumOfCVs());

    int last_cv = -1;

    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++){

        // do we have enough samples?
        int nsamples = Accumulator.GetNumOfSamples(ibin);
        if( nsamples < Options.GetOptLimit() ) continue;

    // write block delimiter - required by GNUPlot
        if(Options.GetOptNoGNUPlot() == false) {
            int ncvs = Accumulator.GetNumOfCVs();
            Accumulator.GetIPoint(ibin,ipos);

            if( (last_cv >= 0) && (ipos[ncvs-1] != last_cv + 1) ){
                if(fprintf(OutputFile,"\n") <= 0) {
                    CSmallString error;
                    error << "unable to write to output (" << strerror(errno) << ")";
                    ES_ERROR(error);
                    return(false);
                }
            }
            last_cv = ipos[ncvs-1];
        }

        Accumulator.GetPoint(ibin,pos);

        CSmallString xformat,sformat;
        xformat = Options.GetOptIXFormat() + " ";

        // print point position
        for(int i=0; i < Accumulator.GetNumOfCVs(); i++) {
            double xvalue = pos[i];
            if(fprintf(OutputFile,xformat,xvalue) <= 0) {
                CSmallString error;
                error << "unable to write to output (" << strerror(errno) << ")";
                ES_ERROR(error);
                return(false);
            }
        }

        sformat = Options.GetOptOSFormat() + "\n";
        // and value
        if(fprintf(OutputFile,sformat,nsamples) <= 0) {
            CSmallString error;
            error << "unable to write to output (" << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFSamples::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-samples terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

