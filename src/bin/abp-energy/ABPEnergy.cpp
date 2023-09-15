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
#include <ABPProxy_dG.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

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
    int State = 1;

// -----------------------------------------------------------------------------
// setup accu, energy proxy, and output FES
    Accu        = CPMFAccumulatorPtr(new CPMFAccumulator);
    FES         = CEnergySurfacePtr(new CEnergySurface);
    EneProxy    = CABPProxy_dG_Ptr(new CABPProxy_dG);

// load ABP accumulator
    vout << endl;
    vout << format("%02d:Loading ABP accumulator: %s")%State%string(Options.GetArgInput()) << endl;
    try {
        Accu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABP accumulator file");
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

// calculate mollified FES
    for(int i=0; i < Accu->GetNumOfBins(); i++){
        FES->SetNumOfSamples(i,EneProxy->GetNSamples(i));
        FES->SetEnergy(i, EneProxy->GetValue(i,E_PROXY_VALUE) );
    }

    vout << format("   Mollified FES SigmaF2         = %10.5f")%FES->GetSigmaF2() << endl;

    if( Options.GetOptMode() == "rl" ){
        //deconvolution
        RunRLDeconvolution();
        vout << "   -------------------------------------" << endl;
    } else if( Options.GetOptMode() == "mollified" ) {
        // nothing to be here
    } else {
        CSmallString error;
        error << "unsupported mode: " << Options.GetOptMode();
        ES_ERROR(error);
        return(false);
    }

// should we apply offset?
    if(Options.GetOptKeepFloating() == false) {
        // get value of global minimum
        double offset = - FES->GetGlobalMinimumValue();
        FES->ApplyOffset(offset);
    }

    vout << format("   Final FES SigmaF2             = %10.5f")%FES->GetSigmaF2() << endl;

// post-processing
    FES->ApplyOffset(Options.GetOptOffset());

    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            FES->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            FES->AdaptUnsampledToMaxEnergy();
        }
    }

// -----------------------------------------------------------------------------
// print energy surface

    vout << endl;
    vout << format("%02d:Writing results to file: %s")%State%string(Options.GetArgOutput()) << endl;
    State++;

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

void CABPEnergy::RunRLDeconvolution(void)
{
    vout << "   Initiating Lucy-Richardson deconvolution ..." << endl;
    // POP is de-convoluted
    // create copy of original data
    CPMFAccuDataPtr upop = Accu->GetSectionData("POP");
    Accu->DeleteSectionData("DPOP"); // this will not be valid after deconvolution

    // normalize POP
    double popsum = 0.0;
    for(int ibin = 0; ibin < upop->GetNumOfBins(); ibin++){
        popsum = popsum + upop->GetData(ibin);
    }
    if( popsum != 0.0 ){
        for(int ibin = 0; ibin < upop->GetNumOfBins(); ibin++){
            double pop = upop->GetData(ibin) / popsum;
            upop->SetData(ibin,pop);
        }
    }

    CPMFAccuDataPtr dpop = upop->Duplicate();

    for(int i=0; i < Options.GetOptRLIter(); i++){

    // run iteration
        for(int jbin=0; jbin < upop->GetNumOfBins(); jbin++){
            double u  = upop->GetData(jbin);
            double f = 0.0;
            for(int ibin=0; ibin < upop->GetNumOfBins(); ibin++){
                double di  = dpop->GetData(ibin);
                double pij = PSF(ibin,jbin);
                double ci  = 0.0;
                for(int kbin=0; kbin < upop->GetNumOfBins(); kbin++){
                    double uk  = upop->GetData(kbin);
                    double pik = PSF(kbin,ibin);
                    ci = ci + uk * pik;
                }
                f = f + di * pij / ci;
            }
            u = u * f;
            upop->SetData(jbin,u);
        }

    // calculate deconvoluted FES
        for(int i=0; i < Accu->GetNumOfBins(); i++){
            FES->SetNumOfSamples(i,EneProxy->GetNSamples(i));
            FES->SetEnergy(i, EneProxy->GetValue(i,E_PROXY_VALUE) );
        }
        vout << format("   #%03d Deconvoluted FES SigmaF2 = %10.5f")%(i+1)%FES->GetSigmaF2() << endl;
    }
}

//------------------------------------------------------------------------------

double CABPEnergy::PSF(int ibin,int jbin)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;

    ipos.CreateVector(Accu->GetNumOfCVs());
    Accu->GetPoint(ibin,ipos);

    jpos.CreateVector(Accu->GetNumOfCVs());
    Accu->GetPoint(jbin,jpos);

    CPMFAccuDataPtr widths = Accu->GetSectionData("WIDTHS");

    double arg = 0.0;
    double dnorm = 1.0;

    for(int icv=0; icv < Accu->GetNumOfCVs(); icv++){
        double diff     = Accu->GetCV(icv)->GetDifference(ipos[icv],jpos[icv]);
        double width    = widths->GetData(icv);
        arg             = arg + diff*diff / (width*width);
        dnorm           = dnorm / (width * sqrt(2.0*M_PI));
    }

    return( dnorm * exp( - 0.5*arg) );
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




