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

// load accumulator
    vout << endl;
    vout << format("%02d:Loading MTD accumulator/trajectory: %s")%State%string(Options.GetArgInput()) << endl;
    State++;

    try {
        Accus = CPMFAccumulator::LoadFinalSnapshots(InputFile,Options.GetOptSmooth());
    } catch(...) {
        ES_ERROR("unable to load the input MTD accumulator file");
        return(false);
    }

    vout << format("   Number of PMF accumulators read:      %ld")%Accus.size() << endl;

    CPMFAccumulatorPtr accu = Accus.back();
    if( accu == NULL ) {
        RUNTIME_ERROR("no PMF accumulators were read");
    }

    accu->PrintInfo(vout);

// -----------------------------------------------------------------------------
// calculate FES
    vout << endl;
    vout << format("%02d:Calculating FES ...")%State << endl;
    State++;

    CEnergySurfacePtr tmp_FES   = CEnergySurfacePtr(new CEnergySurface);
    CMTDProxy_dG_Ptr  eneproxy  = CMTDProxy_dG_Ptr(new CMTDProxy_dG);
    FES                         = CEnergySurfacePtr(new CEnergySurface);

// allocate surface
    FES->Allocate(accu);
    tmp_FES->Allocate(accu);

    std::list<CPMFAccumulatorPtr>::iterator it = Accus.begin();
    std::list<CPMFAccumulatorPtr>::iterator ie = Accus.end();

    int nf = 1;
    while( it != ie ){
        accu = *it;
        eneproxy->Init(accu);
    // calculate energy
        for(int j=0; j < accu->GetNumOfBins(); j++){
            tmp_FES->SetNumOfSamples(j,eneproxy->GetNumOfSamples(j));
            tmp_FES->SetEnergy(j, eneproxy->GetValue(j,E_PROXY_VALUE) );
        }
        int numofhills = accu->GetTotalNumOfSamples();
        CSmallString type = "MTD";
        if( eneproxy->IsWTMeta() ){
            type << "-WT";
        }
        vout << format("   #%03d FES SigmaF2     = %12.5lf, #hills = %10d, Type: %s")%nf%tmp_FES->GetSigmaF2All()%numofhills%type << endl;
        FES->AddFES(tmp_FES);
        it++;
        nf++;
    }
    vout << "   -----------------------------------" << endl;
    double num = Accus.size();
    if( num > 0 ){
        FES->DivideFES(num);
    }

    vout << format("   Smoothed FES SigmaF2 = %12.5f")%FES->GetSigmaF2All() << endl;
    vout << "   >> FES post processing ..." << endl;

// should we apply offset?
    if(Options.GetOptKeepFloating() == false) {
        // get value of global minimum
        double offset = - FES->GetGlobalMinimumValue();
        FES->ApplyOffset(offset);
    }

// post-processing
    FES->ApplyOffset(Options.GetOptOffset());

    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            FES->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            FES->AdaptUnsampledToMaxEnergy();
        }
    }

    vout << format("   Final FES SigmaF2    = %12.5f")%FES->GetSigmaF2All() << endl;

// -----------------------------------------------------------------------------
// print energy surface

    vout << endl;
    vout << format("%02d:Writing results to file: %s")%State%string(Options.GetArgOutput()) << endl;
    State++;

// print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(OutputFile);
        accu = Accus.back();
        accu->PrintInfo(OutputFile);
    }

    CESPrinter printer;

    if(Options.GetOptPrintAll()) {
        printer.SetSampleLimit(0);
    } else {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());

    printer.SetPrintedES(FES);
    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());

    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
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




