// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "ABFMask.hpp"
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFMask)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFMask::CABFMask(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFMask::Init(int argc,char* argv[])
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
    vout << "# abf-mask (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accu file (in)    : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if(Options.GetArgMaskOutputName() != "-") {
        vout << "# Mask file (out)       : " << Options.GetArgMaskOutputName() << endl;
    } else {
        vout << "# Mask file (out)       : - (standard output)" << endl;
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
    vout << "# ------------------------------------------------" << endl;
    vout << "# Output mask format    : " << Options.GetOptOutputFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgMaskOutputName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFMask::Run(void)
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

// initializing mask ----------------------------
    vout << endl;
    vout << "4) Initializing mask" << endl;
    vout << "   Sample limit: " << Options.GetOptLimit() << endl;
    vout << "   Energy limit: " << Options.GetOptMaxEnergy() << endl;
    int masked = 0;
    for(int i = 0; i < Accumulator.GetNumberOfBins(); i++){
        double energy = fes.GetEnergy(i);
        if( (energy > Options.GetOptMaxEnergy()) || (Accumulator.GetNumberOfABFSamples(i) < Options.GetOptLimit()) ){
            Accumulator.SetMaskWeight(i,0.0);
            masked++;
        }
    }
    double maxbins = Accumulator.GetNumberOfBins();
    if( maxbins > 0 ){
        vout << "   Masked area: " << setw(5) << setprecision(1) << fixed << masked/maxbins*100 <<"%" << endl;
    }
    vout << "   Done" << endl;

// print result ---------------------------------
    vout << endl;
    vout << "5) Writing mask to file: " << Options.GetArgMaskOutputName() << endl;

    if(Options.GetOptOutputFormat() == "mask") {
        Accumulator.SaveMask(OutputFile);
    }

    if(Options.GetOptOutputFormat() == "gnuplot") {
        Accumulator.SaveMaskGnuPlot(OutputFile);
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

void CABFMask::PrepareAccumulator(void)
{
    for(int icoord=0; icoord < Accumulator.GetNumberOfCoords(); icoord++) {
        for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {

            double value = 0.0;
            double sum = 0.0;
            int    nsamples = 0;
            int    newsamples = 0;

            nsamples = Accumulator.GetNumberOfABFSamples(ibin);
            // abf force
            sum = Accumulator.GetABFForceSum(icoord,ibin);

            if((nsamples > 0) && (nsamples > Options.GetOptLimit())) {
                // calculate average
                value = sum / nsamples;
                newsamples = nsamples;
            }
            Accumulator.SetABFForceSum(icoord,ibin,value);
            Accumulator.SetNumberOfABFSamples(ibin,newsamples);
        }
    }

    // calculate sampled area
    double maxbins = Accumulator.GetNumberOfBins();
    double sampled = Accumulator.GetNumberOfBinsWithABFLimit(1);
    if( maxbins > 0 ){
        vout << "   Sampled area: " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFMask::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-mask terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

