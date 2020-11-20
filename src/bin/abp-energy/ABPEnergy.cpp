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
#include "ABPEnergy.hpp"
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABPEnergy)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABPEnergy::CABPEnergy(void)
{
    InputFile = NULL;
    OwnInputFile = false;
    OutputFile = NULL;
    OwnOutputFile = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPEnergy::Init(int argc,char* argv[])
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
    vout << "# abp-energy (PMFLib utility) started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    vout << "#" << endl;
    if(Options.GetArgInput() != "-") {
    vout << "# ABP accu file (in)    : " << Options.GetArgInput() << endl;
    } else {
    vout << "# ABP accu file (in)    : " << Options.GetArgInput() << " (standard input)" << endl;
    }
    if(Options.GetArgOutput() != "-") {
    vout << "# Free energy file (out): " << Options.GetArgOutput() << endl;
    } else {
    vout << "# Free energy file (out): " << Options.GetArgOutput() << " (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# Mode                  : " << Options.GetOptMode() << endl;
    vout << "# Temperature           : " << Options.GetOptTemperature() << endl;
    if(Options.GetOptLimit() == 0) {
    vout << "# Sampling limit        : all bins will be taken into account" << endl;
    } else {
    vout << "# Sampling limit        : " << Options.GetOptLimit() << endl;
    }

    vout << "# ------------------------------------------------" << endl;
    vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# Include bin statuses  : " << bool_to_str(Options.GetOptIncludeBinStat()) << endl;
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

bool CABPEnergy::Run(void)
{
// load history list
    try {
        ABPAccumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABP accumulator file");
        return(false);
    }

// print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        // FIXME - FILL info
    }

// allocate surface
    FES.Allocate(&ABPAccumulator);

// calculate mollified FES
    vout << endl;
    CalculateFESMollified();
    vout << " Mollified FES SigmaF2 = " << setprecision(5) << FES.GetSigmaF2() << endl;

    if( Options.GetOptMode() == "rl" ){
        //deconvolution
        // FIXME
    } else if( Options.GetOptMode() == "mollified" ) {
        // nothing to be here
    } else {
        CSmallString error;
        error << "unsupported mode: " << Options.GetOptMode();
        ES_ERROR(error);
        return(false);
    }

    vout << " Final FES SigmaF2     = " << setprecision(5) << FES.GetSigmaF2() << endl;

// post-processing
    FES.ApplyOffset(Options.GetOptOffset());

    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            FES.AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            FES.AdaptUnsampledToMaxEnergy();
        }
    }

// print energy surface
    CESPrinter printer;

    if(Options.GetOptPrintAll()) {
        printer.SetSampleLimit(0);
    } else {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());

    printer.SetPrintedES(&FES);
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

void CABPEnergy::Finalize(void)
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
    vout << "# abp-energy terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPEnergy::CalculateFESMollified(void)
{
//! gas constant 8.314 472(15) J mol-1 K-1
//real(PMFDP), parameter  :: PMF_Rgas     = 0.0019872065d0     ! kcal mol-1 K-1 = 8.314 472 / 4184

    const double Rgas = 0.0019872065;

    double m = 1.0;
    for(int i=0; i < ABPAccumulator.GetNumberOfBins(); i++){
        double pop = ABPAccumulator.GetPop(i);
        if( pop > m ) m = pop;
    }
    for(int i=0; i < ABPAccumulator.GetNumberOfBins(); i++){
        double pop = ABPAccumulator.GetPop(i);
        double ene = - Options.GetOptTemperature()*Rgas*log(pop/m);
        FES.SetEnergy(i,ene);
        FES.SetNumOfSamples(i,ABPAccumulator.GetNumberOfABPSamples(i));
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




