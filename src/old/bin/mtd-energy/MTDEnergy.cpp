// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

#include <stdio.h>
#include <math.h>
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <ESPrinter.hpp>
#include "MTDEnergy.hpp"

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CMTDEnergy)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDEnergy::CMTDEnergy(void)
{
    InputFile = NULL;
    OwnInputFile = false;
    OutputFile = NULL;
    OwnOutputFile = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDEnergy::Init(int argc,char* argv[])
{
    if((InputFile != NULL) || (OutputFile != NULL)) return(false);     // files already opened

// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to terminal stream and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

// print header if requested --------------------
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# mtd-energy (PMFLib utility) started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    vout << "#" << endl;
    if(Options.GetArgInput() != "-") {
        vout << "# MTD history file (in) : " << Options.GetArgInput() << endl;
    } else {
        vout << "# MTD history file (in) : " << Options.GetArgInput() << " (standard input)" << endl;
    }
    if(Options.GetArgOutput() != "-") {
        vout << "# Free energy file (out): " << Options.GetArgOutput() << endl;
    } else {
        vout << "# Free energy file (out): " << Options.GetArgOutput() << " (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptTime() == 0) {
        vout << "# Used MTD time         : infinity" << endl;
    } else {
        vout << "# Used MTD time         : " << Options.GetOptTime() << endl;
    }

    if(Options.IsOptOffsetSet() == false) {
        vout << "# Global minimum        : floating" << endl;
    } else {
        vout << "# Global minimum        : " << Options.GetOptOffset() << endl;
    }

    if(Options.GetOptSmooth() == 0) {
        vout << "# Smoothing from        : no smoothing" << endl;
    } else {
        vout << "# Smoothing from        : " << Options.GetOptSmooth() << endl;
    }
    vout << "# Print standard dev    : " << bool_to_str(Options.GetOptPrintSD());

    vout << "# ------------------------------------------------" << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;

// open files -----------------------------------
    if(Options.GetArgInput() == "-") {
        InputFile = stdin;
        OwnInputFile = false;
    } else {
        InputFile = fopen(Options.GetArgInput(),"r");
        if(InputFile == NULL) {
            CSmallString error;
            error << "unable to open input file '" << Options.GetArgInput() << "'";
            ES_ERROR(error);
            return(SO_USER_ERROR);
        }
        OwnInputFile = true;
    }

    if(Options.GetArgOutput() == "-") {
        OutputFile = stdout;
        OwnOutputFile = false;
    } else {
        OutputFile = fopen(Options.GetArgOutput(),"w");
        if(OutputFile == NULL) {
            CSmallString error;
            error << "unable to open input file '" << Options.GetArgOutput() << "'";
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

bool CMTDEnergy::Run(void)
{
// load history list
    try {
        MTDHistory.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input MTD history file");
        return(false);
    }

// print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        fprintf(OutputFile,"# Number of coordinates : %d\n",MTDHistory.GetNumOfCVs());
        fprintf(OutputFile,"# Total number of hills : %d\n",MTDHistory.GetNumOfHills());

        if(Options.GetOptTime() == 0) {
            fprintf(OutputFile,"# Used MTD time         : %d\n",MTDHistory.GetNumOfHills());
        } else {
            fprintf(OutputFile,"# Used MTD time         : %d\n",Options.GetOptTime());
        }
        if(Options.GetOptSmooth() > 0) {
            fprintf(OutputFile,"# Smooothed from        : %d\n",Options.GetOptSmooth());
        } else {
            fprintf(OutputFile,"# Smooothed from        : - no smoothing-\n");
        }
    }

// allocate surface
    EnergySurface.Allocate(&MTDHistory);

// calculate energy -----------------------------
    if(Options.GetOptSmooth() <= 0) {
        if(CalculateFES() == false) return(false);
    } else {
        if(CalculateSmoothedFES() == false) return(false);
    }

// should we apply offset?
    if(Options.IsOptOffsetSet() == true) {
        // get value of global minumum
        double offset = Options.GetOptOffset() - EnergySurface.GetGlobalMinimumValue();
        EnergySurface.ApplyOffset(offset);
    }

// print energy surface
    CESPrinter printer;

    printer.SetPrintedES(&EnergySurface);
    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    }
    if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    }
    if(Options.GetOptOutputFormat() == "fes") {
        printer.SetOutputFormat(EESPF_PMF_FES);
    }


    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to print energy surface to output file");
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDEnergy::Finalize(void)
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
    vout << "# mtd-energy terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDEnergy::CalculateFES(void)
{
// calculate surface
    // FIXME
    // EnergySurface.CalculateFES(MTDHistory,Options.GetOptTime());
    return(true);
}

//------------------------------------------------------------------------------

bool CMTDEnergy::CalculateSmoothedFES(void)
{
    CEnergySurface intermediate;

    intermediate.Allocate(&MTDHistory);

    int stop;
    if(Options.GetOptTime() == 0) {
        stop = MTDHistory.GetNumOfHills();
    } else {
        stop = Options.GetOptTime();
    }

// cumulate surfaces
    for(int time=Options.GetOptSmooth(); time <= stop; time++) {
        vout << " Processing time : " << time << endl;
        // FIXME
        // intermediate.CalculateFES(MTDHistory,time);
        EnergySurface += intermediate;
    }

    int length = stop - Options.GetOptSmooth() + 1;
    if(length <= 0) {
        ES_ERROR("total number of snapshots for smoothing is zero");
        return(false);
    }

// calculate average
    EnergySurface /= (double)length;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




