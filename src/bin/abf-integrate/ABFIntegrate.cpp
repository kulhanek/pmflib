// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Martin Petrek, petrek@chemi.muni.cz
//                       Petr Kulhanek, kulhanek@enzim.hu
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

#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <ABFIntegratorFD.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include "ABFIntegrate.hpp"
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegrate::CABFIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFIntegrate::Init(int argc,char* argv[])
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
    vout << "# abf-integrate (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accu file (in)    : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if(Options.GetArgFEOutputName() != "-") {
        vout << "# Free energy file (out): " << Options.GetArgFEOutputName() << endl;
    } else {
        vout << "# Free energy file (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Limit                 : all bins will be taken into account" << endl;
    } else {
        vout << "# Limit                 : " << Options.GetOptLimit() << endl;
    }
    vout << "# FD order              : " << Options.GetOptFDOrder() << endl;
    vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# Periodicity           : " << bool_to_str(Options.GetOptPeriodicity()) << endl;
    vout << "# Integrate errors      : " << bool_to_str(Options.GetOptErrors()) << endl;
    vout << "# ------------------------------------------------" << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgFEOutputName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::Run(void)
{
// load accumulator
    vout << endl;
    vout << "1) Loading ABF accumulator: " << Options.GetArgABFAccuName() << endl;
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done" << endl;

    // print CVS info
    Accumulator.PrintCVSInfo(vout);

// prepare accumulator --------------------------
    vout << endl;
    vout << "2) Preparing ABF accumulator for integration"<< endl;
    PrepareAccumulator();
    vout << "   Done" << endl;

// integrate data ------------------------------
    vout << endl;
    vout << "3) ABF accumulator integration"<< endl;
    CABFIntegratorFD   integrator;
    CEnergySurface     fes;

    integrator.SetVerbosity(Options.GetOptVerbose());
    integrator.SetPeriodicity(Options.GetOptPeriodicity());
    integrator.SetFDOrder(Options.GetOptFDOrder());

    integrator.SetInputABFAccumulator(&Accumulator);
    integrator.SetOutputFESurface(&fes);

    if(integrator.Integrate() == false) {
        ES_ERROR("unable to prepare ABF accumulator");
        return(false);
    }
    vout << "   Done" << endl;

// print result ---------------------------------
    vout << endl;
    vout << "4) Writing results to file: " << Options.GetArgFEOutputName() << endl;
    CESPrinter printer;

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

    if(Options.GetOptPrintWithLimit()) {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetPrintedES(&fes);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output free energy file");
        return(false);
    }
    vout << "   Done" << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero
//    b) accumulated sums/sum squares will be recalculated to averages

void CABFIntegrate::PrepareAccumulator(void)
{
// print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        fprintf(OutputFile,"# data integrated by reverse finite difference method\n");
        fprintf(OutputFile,"# Number of coordinates : %d\n",Accumulator.GetNumberOfCoords());
        fprintf(OutputFile,"# Total number of bins  : %d\n",Accumulator.GetNumberOfBins());
        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        fprintf(OutputFile,"# Periodicity           : %s\n",(const char*)bool_to_str(Options.GetOptPeriodicity()));
        fprintf(OutputFile,"# Integrate errors      : %s\n",(const char*)bool_to_str(Options.GetOptErrors()));
        fprintf(OutputFile,"# FD order              : %d\n",Options.GetOptFDOrder());
    }

    for(int icoord=0; icoord < Accumulator.GetNumberOfCoords(); icoord++) {
        for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {

            double value = 0.0;
            double sum = 0.0;
            double sum_square = 0.0;
            int    nsamples = 0;
            int    newsamples = 0;

            nsamples = Accumulator.GetNumberOfABFSamples(ibin);
            // abf force
            sum = Accumulator.GetABFForceSum(icoord,ibin);
            sum_square = Accumulator.GetABFForceSquareSum(icoord,ibin);

            if((nsamples > 0) && (nsamples > Options.GetOptLimit())) {
                if(Options.GetOptErrors() == false) {
                    // calculate average
                    value = sum / nsamples;
                } else {
                    // calculate error of average
                    double sq = nsamples*sum_square - sum*sum;
                    if(sq > 0) {
                        sq = sqrt(sq) / nsamples;
                    } else {
                        sq = 0.0;
                    }
                    value = sq / sqrt((double)nsamples);
                }
                newsamples = nsamples;
            }
            Accumulator.SetABFForceSum(icoord,ibin,value);
            Accumulator.SetNumberOfABFSamples(ibin,newsamples);
        }
    }

    // calculate sampled area
    double maxbins = Accumulator.GetNumberOfBins();
    double sampled = Accumulator.GetNumberOfBinsWithABFLimit(Options.GetOptLimit());
    if( maxbins > 0 ){
        vout << "   Sampled area: " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegrate::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-integrate terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

