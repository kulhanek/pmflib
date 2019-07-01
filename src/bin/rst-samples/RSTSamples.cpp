// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "RSTSamples.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CRSTSamples)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CRSTSamples::CRSTSamples(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CRSTSamples::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CUMIntOpts
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
    vout << "# rst-samples (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgRSTHistName() != "-") {
        vout << "# RST accumulator file (in) : " << Options.GetArgRSTHistName() << endl;
    } else {
        vout << "# RST accumulator file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgOutputName() != "-") {
        vout << "# Sample file (out)         : " << Options.GetArgOutputName() << endl;
    } else {
        vout << "# Sample file (out)         : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Limit                     : all bins will be printed" << endl;
    } else {
        vout << "# Limit                     : " << Options.GetOptLimit() << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# No GNUPlot delimiters     : " << bool_to_str(Options.GetOptNoGNUPlot()) << endl;
    vout << "# X format                  : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format                  : " << Options.GetOptOSFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

// open files -----------------------------------
    if(Options.GetArgRSTHistName() == "-") {
        InputFile = stdin;
        OwnInputFile = false;
    } else {
        InputFile = fopen(Options.GetArgRSTHistName(),"r");
        if(InputFile == NULL) {
            CSmallString error;
            error << "unable to open input file (" << strerror(errno) << ")";
            ES_ERROR(error);
            return(SO_USER_ERROR);
        }
        OwnInputFile = true;
    }

    if(Options.GetArgOutputName() == "-") {
        OutputFile = stdout;
        OwnOutputFile = false;
    } else {
        OutputFile = fopen(Options.GetArgOutputName(),"w");
        if(OutputFile == NULL) {
            CSmallString error;
            error << "unable to open output file (" << strerror(errno) << ")";
            ES_ERROR(error);
            if((OwnInputFile == true) && (InputFile != NULL)) {
                fclose(InputFile);
                InputFile = NULL;
                OwnInputFile = false;
            }
            return(SO_USER_ERROR);
        }
        OwnOutputFile = true;
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CRSTSamples::Run(void)
{
// load history list
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to read RST accumulator file");
        return(false);
    }

    Point.CreateVector(Accumulator.GetNumberOfCoords());

// print samples
    if(PrintSamples(0) == false) {
        ES_ERROR("unable to print samples");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CRSTSamples::PrintSamples(int cv)
{
    if(cv >= Accumulator.GetNumberOfCoords()) {
        // get number of samples
        int value = Accumulator.GetNumberOfSamples(Point);

        if(value < Options.GetOptLimit()) return(true);      // not enough samples

        CSmallString xformat,sformat;
        xformat = Options.GetOptIXFormat() + " ";
        sformat = Options.GetOptOSFormat() + "\n";

        // print point position
        for(int i=0; i < Accumulator.GetNumberOfCoords(); i++) {
            const CColVariable* p_coord = Accumulator.GetCoordinate(i);
            double xvalue = p_coord->GetValue(Point[i]);
            if(fprintf(OutputFile,xformat,xvalue) <= 0) {
                CSmallString error;
                error << "unable to write to output (" << strerror(errno) << ")";
                ES_ERROR(error);
                return(false);
            }
        }
        // and value
        if(fprintf(OutputFile,sformat,value) <= 0) {
            CSmallString error;
            error << "unable to write to output (" << strerror(errno) << ")";
            ES_ERROR(error);
            return(false);
        }
        return(true);
    }

    const CColVariable* p_coord = Accumulator.GetCoordinate(cv);

// cycle through variable
    for(unsigned int i = 0; i < p_coord->GetNumberOfBins(); i++) {
        Point[cv] = i;
        if(PrintSamples(cv+1) == false) return(false);
    }

// write block delimiter - required by GNUPlot
    if(Options.GetOptNoGNUPlot() == false) {
        if(fprintf(OutputFile,"\n") <= 0) {
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

void CRSTSamples::Finalize(void)
{
// close files if they are own by program
    if((OwnInputFile == true) && (InputFile != NULL)) {
        fclose(InputFile);
        InputFile = NULL;
        OwnInputFile = false;
    }
    if((OwnOutputFile == true) && (OutputFile != NULL)) {
        fclose(OutputFile);
        OutputFile = NULL;
        OwnOutputFile = false;
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# rst-samples terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

