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

#include <math.h>
#include <errno.h>
#include "ABFDerivatives.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFDerivatives)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFDerivatives::CABFDerivatives(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFDerivatives::Init(int argc,char* argv[])
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
    vout << "# abf-derivatives (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accumulator file (in) : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# ABF accumulator file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgOutputName() != "-") {
        vout << "# Derivatives file (out)    : " << Options.GetArgOutputName() << endl;
    } else {
        vout << "# Derivatives file (out)    : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Limit                     : all bins will be printed" << endl;
    } else {
        vout << "# Limit                     : " << Options.GetOptLimit() << endl;
    }
    vout << "# CV item                   : " << Options.GetOptItem() << endl;
    vout << "# Print sigmas              : " << bool_to_str(Options.GetOptSigma()) << endl;
    vout << "# Print errors              : " << bool_to_str(Options.GetOptError()) << endl;
    vout << "# Number of corr. samples   : " << Options.GetOptNCorr() << endl;
    vout << "# ------------------------------------------------" << endl;
    vout << "# No GNUPlot delimiters     : " << bool_to_str(Options.GetOptNoGNUPlot()) << endl;
    vout << "# No header to output       : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# X format                  : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format                  : " << Options.GetOptOSFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
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

bool CABFDerivatives::Run(void)
{
// load accumulators
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }

// print header
    if(Options.GetOptNoHeader() == false) {
        fprintf(OutputFile,"# derivatives\n");
        fprintf(OutputFile,"# Number of coordinates : %d\n",Accumulator.GetNumberOfCoords());
        fprintf(OutputFile,"# Total number of bins  : %d\n",Accumulator.GetNumberOfBins());
        fprintf(OutputFile,"# CV item               : %d\n",Options.GetOptItem());
        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        fprintf(OutputFile,"# Include std. dev.     : %s\n",(const char*)bool_to_str(Options.GetOptSigma()));
        fprintf(OutputFile,"# Include std. err.     : %s\n",(const char*)bool_to_str(Options.GetOptError()));
        fprintf(OutputFile,"# Number of corr. sam.  : %5.3f\n",Options.GetOptNCorr());
    }

// print samples
    if(PrintDerivatives() == false) {
        ES_ERROR("unable to print derivatives to output");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFDerivatives::PrintDerivatives(void)
{
    if( (Options.GetOptItem() < 0) || (Options.GetOptItem() > Accumulator.GetNumberOfCoords())){
        ES_ERROR("requested CV is out of range");
        return(false);
    }

    CSimpleVector<double>   pos;
    CSimpleVector<int>      ipos;

    pos.CreateVector(Accumulator.GetNumberOfCoords());
    ipos.CreateVector(Accumulator.GetNumberOfCoords());

    int last_cv = -1;

    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++){

        // do we have enough samples?
        double nsamples = Accumulator.GetNumberOfABFSamples(ibin);
        if( nsamples < Options.GetOptLimit() ) continue;

    // write block delimiter - required by GNUPlot
        if(Options.GetOptNoGNUPlot() == false) {
            int ncvs = Accumulator.GetNumberOfCoords();
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
        for(int i=0; i < Accumulator.GetNumberOfCoords(); i++) {
            double xvalue = pos[i];
            if(fprintf(OutputFile,xformat,xvalue) <= 0) {
                CSmallString error;
                error << "unable to write to output (" << strerror(errno) << ")";
                ES_ERROR(error);
                return(false);
            }
        }

        for(int i=0; i < Accumulator.GetNumberOfCoords(); i++) {
            double value = 0.0;
            double sigma = 0.0;
            double error = 0.0;
            double sum = 0.0;
            double sum_square = 0.0;

            // abf force
            sum = Accumulator.GetABFForceSum(i,ibin);
            sum_square = Accumulator.GetABFForceSquareSum(i,ibin);

            if(nsamples > 0) {
                // calculate average
                value = sum / nsamples;
                if( (Options.GetOptSigma() == true) || (Options.GetOptError() == true)  ){
                    // calculate sigma of samples
                    double sq = nsamples*sum_square - sum*sum;
                    if(sq > 0) {
                        sq = sqrt(sq) / nsamples;
                    } else {
                        sq = 0.0;
                    }
                    sigma = sq;
                }
                if( Options.GetOptError() == true ){
                    // calculate error of average
                    error = sigma / sqrt(nsamples);
                }
            }

            if( (Options.GetOptItem() == 0) || (Options.GetOptItem() == i+1) ){

                sformat = Options.GetOptOSFormat() + " ";
                // and value and optionaly sigma and error
                if(fprintf(OutputFile,sformat,value) <= 0) {
                    CSmallString error;
                    error << "unable to write to output (" << strerror(errno) << ")";
                    ES_ERROR(error);
                    return(false);
                }
                if( Options.GetOptSigma() == true ) {
                    if(fprintf(OutputFile,sformat,sigma) <= 0) {
                        CSmallString error;
                        error << "unable to write to output (" << strerror(errno) << ")";
                        ES_ERROR(error);
                        return(false);
                    }
                }
                if( Options.GetOptError() == true ) {
                    if(fprintf(OutputFile,sformat,error) <= 0) {
                        CSmallString error;
                        error << "unable to write to output (" << strerror(errno) << ")";
                        ES_ERROR(error);
                        return(false);
                    }
                }
            }
        }
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

void CABFDerivatives::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-derivatives terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

