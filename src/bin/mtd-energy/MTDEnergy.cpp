// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <iomanip>
#include "MTDEnergy.hpp"
#include <boost/format.hpp>
#include <MTDProxy_dG.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

MAIN_ENTRY(CMTDEnergy)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDEnergy::CMTDEnergy(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDEnergy::Init(int argc,char* argv[])
{

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
        vout << "# MTD accu file (in)    : " << Options.GetArgInput() << endl;
    } else {
        vout << "# MTD accu file (in)    : " << Options.GetArgInput() << " (standard input)" << endl;
    }
    if(Options.GetArgOutput() != "-") {
        vout << "# Free energy file (out): " << Options.GetArgOutput() << endl;
    } else {
        vout << "# Free energy file (out): " << Options.GetArgOutput() << " (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;

// open files -----------------------------------
    if( InputFile.Open(Options.GetArgInput(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
   if( OutputFile.Open(Options.GetArgOutput(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDEnergy::Run(void)
{
    int State = 1;

// -----------------------------------------------------------------------------
// setup accu, energy proxy, and output FES
    Accu        = CPMFAccumulatorPtr(new CPMFAccumulator);
    FES         = CEnergySurfacePtr(new CEnergySurface);
    EneProxy    = CMTDProxy_dG_Ptr(new CMTDProxy_dG);

// load accumulator
    vout << endl;
    vout << format("%02d:Loading MTD accumulator: %s")%State%string(Options.GetArgInput()) << endl;
    State++;

    try {
        Accu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input MTD accumulator file");
        return(false);
    }

    Accu->PrintInfo(vout);
    EneProxy->Init(Accu);

// -----------------------------------------------------------------------------
// calculate FES
    vout << endl;
    vout << format("%02d:Calculating FES ...")%State << endl;
    State++;

// print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(OutputFile);
        Accu->PrintInfo(OutputFile);
    }

// allocate surface
    FES->Allocate(Accu);

// calculate energy
    for(int i=0; i < Accu->GetNumOfBins(); i++){
        FES->SetNumOfSamples(i,1);
        FES->SetEnergy(i, EneProxy->GetValue(i,E_PROXY_VALUE) );
    }

// should we apply offset?
    if(Options.GetOptKeepFloating() == false) {
        // get value of global minumum
        double offset = Options.GetOptOffset() - FES->GetGlobalMinimumValue();
        FES->ApplyOffset(offset);
    }

    vout << endl;
    vout << "Final FES SigmaF2     = " << setprecision(5) << FES->GetSigmaF2() << endl;

// -----------------------------------------------------------------------------
// print energy surface

    vout << endl;
    vout << format("%02d:Writing results to file: %s")%State%string(Options.GetArgOutput()) << endl;
    State++;

    CESPrinter printer;

    printer.SetPrintedES(FES);
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
    InputFile.Close();
    OutputFile.Close();

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




