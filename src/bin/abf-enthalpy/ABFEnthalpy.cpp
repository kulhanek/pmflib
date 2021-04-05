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
#include "ABFEnthalpy.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <ESPrinter.hpp>
#include <iomanip>
#include <ABFEnthalpyGPR.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFEnthalpy)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFEnthalpy::CABFEnthalpy(void)
{
    State = 1;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFEnthalpy::Init(int argc,char* argv[])
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

    ABFAccuName = Options.GetProgArg(0);
    HEOutputName = Options.GetProgArg(1);

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-enthalpy (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(ABFAccuName != "-") {
        vout << "# ABF accumulator file (in) : " << ABFAccuName << endl;
    } else {
        vout << "# ABF accumulator file (in) : - (standard input)" << endl;
    }
    if(HEOutputName != "-") {
        vout << "# Enthalpy file (out)       : " << HEOutputName << endl;
    } else {
        vout << "# Enthalpy file (out)       : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptMethod() == "raw" ) {
        vout << "# Method                    :     RAW (raw data from ABF accumulator)" << endl;
    } else if( Options.GetOptMethod() == "gpr" ) {
        vout << "# Method                    :     GPR (gaussian process filtered data)" << endl;
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }
    if(Options.GetOptLimit() == 0) {
        vout << "# Limit                     : all bins will be printed" << endl;
    } else {
        vout << "# Limit                     : " << Options.GetOptLimit() << endl;
    }
    vout << "# Print errors              : " << bool_to_str(Options.GetOptWithError()) << endl;
    vout << "# Number of corr. samples   : " << Options.GetOptNCorr() << endl;
    vout << "# ------------------------------------------------" << endl;
    vout << "# No header to output       : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# X format                  : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format                  : " << Options.GetOptOEFormat() << endl;
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

bool CABFEnthalpy::Run(void)
{
// load accumulator
    State = 1;

    vout << endl;
    vout << format("%02d:Loading ABF accumulator: %s")%State%string(ABFAccuName) << endl;
    State++;
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    // print CVS info
    Accumulator.PrintCVSInfo(vout);
    // DO NOT SET IT HERE, Ncorr can be GPR hyperparameter
    // Accumulator.SetNCorr(Options.GetOptNCorr());
    HES.Allocate(&Accumulator);
    HES.SetSLevel(Options.GetOptSLevel());

    vout << endl;
    vout << format("%02d:Statistics of input ABF accumulator")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done." << endl;

// sampling limit -------------------------------
    vout << endl;
    vout << format("%02d:Preparing ABF accumulator for analysis (sampling limit)")%State << endl;
    State++;
    PrepareAccumulatorI();
    PrintSampledStat();
    vout << "   Done." << endl;

// -----------------------------------------------------------------------------

    if( Options.GetOptMethod() == "raw" ){
        Accumulator.SetNCorr(Options.GetOptNCorr());
        vout << endl;
        vout << format("%02d:Raw absolute enthalpy")%State << endl;
        GetRawEnthalpy();
        vout << "      SigmaF2   = " << setprecision(5) << HES.GetSigmaF2() << endl;
        vout << "      SigmaF2+  = " << setprecision(5) << HES.GetSigmaF2p() << endl;
        vout << "      SigmaF2-  = " << setprecision(5) << HES.GetSigmaF2m() << endl;
        State++;
        vout << "   Done." << endl;
    } else if ( Options.GetOptMethod() == "gpr" ){
        vout << endl;
        vout << format("%02d:GPR interpolated enthalpy")%State << endl;
        CABFEnthalpyGPR   entgpr;

        entgpr.SetInputABFAccumulator(&Accumulator);
        entgpr.SetOutputHESurface(&HES);

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(entgpr);
        } else {
            entgpr.SetSigmaF2(Options.GetOptSigmaF2());
            entgpr.SetNCorr(Options.GetOptNCorr());
            entgpr.SetWFac(Options.GetOptWFac());
        }

        entgpr.SetIncludeError(Options.GetOptWithError());

        if( Options.IsOptGlobalMinSet() ){
            entgpr.SetGlobalMin(Options.GetOptGlobalMin());
        }

        entgpr.SetRCond(Options.GetOptRCond());
        entgpr.SetLAMethod(Options.GetOptLAMethod());
        entgpr.SetKernel(Options.GetOptGPRKernel());

        if(entgpr.Interpolate(vout) == false) {
            ES_ERROR("unable to interpolate enthalpy");
            return(false);
        }
        State++;
        vout << "   Done." << endl;
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    if( Options.GetOptAbsolute() == false ){
        if( ! Options.IsOptGlobalMinSet() ){
            HES.ApplyOffset(Options.GetOptOffset() - HES.GetGlobalMinimumValue());
        } else {
            HES.ApplyOffset(Options.GetOptOffset());
        }
    }

// print enthalpy
    if(PrintHES() == false) {
        ES_ERROR("unable to print enthalpy");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFEnthalpy::LoadGPRHyprms(CABFEnthalpyGPR& gpr)
{
    ifstream fin;
    fin.open(Options.GetOptLoadHyprms());
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << Options.GetOptLoadHyprms();
        RUNTIME_ERROR(error);
    }

    string line;
    while( getline(fin,line) ){
        // is it comment?
        if( (line.size() > 0) && (line[0] == '#') ) continue;

        // parse line
        stringstream str(line);
        string key, buf;
        double value;
        str >> key >> buf >> value;
        if( ! str ){
            CSmallString error;
            error << "GPR hyperparameters file, unable to decode line: " << line.c_str();
            RUNTIME_ERROR(error);
        }
        if( key == "SigmaF2" ){
            gpr.SetSigmaF2(value);
        } else if( key == "NCorr" ){
            gpr.SetNCorr(value);
        } else if( key.find("NCorr#") != string::npos ) {
            CSmallString error;
            error << "GPR hyperparameters file, unsupported splitted NCorr";
            RUNTIME_ERROR(error);

        } else if( key.find("WFac#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode wfac key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            gpr.SetWFac(cvind,value);
        } else {
            CSmallString error;
            error << "GPR hyperparameters file, unrecognized key: " << key.c_str();
            RUNTIME_ERROR(error);
        }
    }
}

//------------------------------------------------------------------------------

void CABFEnthalpy::GetRawEnthalpy(void)
{
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++){
        int    nsamples = Accumulator.GetNumOfSamples(ibin);
        double ent = Accumulator.GetValue(ibin,EABF_H_VALUE);
        double error = Accumulator.GetValue(ibin,EABF_H_ERROR);
        HES.SetNumOfSamples(ibin,nsamples);
        HES.SetEnergy(ibin,ent);
        HES.SetError(ibin,error);
    }
}

//------------------------------------------------------------------------------

bool CABFEnthalpy::PrintHES(void)
{
    vout << endl;
    vout << format("%02d:Writing results to file: %s")%State%string(HEOutputName) << endl;

    if( OutputFile.Open(HEOutputName,"w") == false ){
        ES_ERROR("unable to open output file");
        return(false);
    }

    State++;
    CESPrinter printer;

    WriteHeader();

    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else if(Options.GetOptOutputFormat() == "fes") {
        printer.SetOutputFormat(EESPF_PMF_FES);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    if(Options.GetOptPrintAll()) {
        printer.SetSampleLimit(0);
    } else {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetPrintedES(&HES);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output enthalpy file");
        return(false);
    }
    vout << "   Done." << endl;

    return(true);
}

//------------------------------------------------------------------------------

void CABFEnthalpy::WriteHeader(void)
{
    // print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        fprintf(OutputFile,"# PMFLib version        : %s\n",LibBuildVersion_PMF);
        fprintf(OutputFile,"# Method                : ");
        if(Options.GetOptMethod() == "raw" ) {
            fprintf(OutputFile,"RAW (raw data from ABF accumulator)\n");
        } else if( Options.GetOptMethod() == "gpr" ) {
            fprintf(OutputFile,"GPR (gaussian process filtered data)\n");
        } else {
            INVALID_ARGUMENT("method - not implemented");
        }
        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        if( Options.GetOptMethod() == "raw" ){
            // nothing to be here
        } else if( Options.GetOptMethod() == "gpr" ){
            fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
            fprintf(OutputFile,"# Linear algebra        : %s\n", (const char*)Options.GetOptLAMethod());
            fprintf(OutputFile,"# SVD rcond             : %5.4e\n",Options.GetOptRCond());
        } else {
            ES_ERROR("not implemented method");
        }

        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        if( Options.IsOptGlobalMinSet() ){
        fprintf(OutputFile,"# Global FES minimum    : %s\n",(const char*)Options.GetOptGlobalMin());
        } else {
        fprintf(OutputFile,"# Global FES minimum    : -auto-\n");
        }
        fprintf(OutputFile,"# Energy offset         : %5.3f\n", Options.GetOptOffset());
        fprintf(OutputFile,"# Number of coordinates : %d\n",Accumulator.GetNumOfCVs());
        fprintf(OutputFile,"# Total number of bins  : %d\n",Accumulator.GetNumOfBins());
    }
}

//------------------------------------------------------------------------------

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero

void CABFEnthalpy::PrepareAccumulatorI(void)
{
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        // erase datapoints not properly sampled, preserve glueing
        if( (Accumulator.GetNumOfSamples(ibin) >= 0) && (Accumulator.GetNumOfSamples(ibin) <= Options.GetOptLimit()) ) {
            Accumulator.SetNumberOfABFSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CABFEnthalpy::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = Accumulator.GetNumOfBins();
    int    sampled = 0;
    int    holes = 0;
    int    glued = 0;
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( Accumulator.GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( Accumulator.GetNumOfSamples(ibin) < 0 ) {
            glued++;
        }
        if( Accumulator.GetNumOfSamples(ibin) == -1 ) {
            holes++;
        }
    }
    if( (maxbins > 0) && (glued != 0) ){
        vout << "      Sampled area:               "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
    if( glued > 0 ){
        vout << "      All inter/extrapolated area:"
             << setw(6) << glued << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << glued/maxbins*100 <<"%" << endl;
    }
    if( holes > 0 ){
        vout << "         Interpolated area:       "
             << setw(6) << holes << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << holes/maxbins*100 <<"%" << endl;
    }
    if( (glued-holes) > 0 ){
        vout << "         Extrapolated area:       "
             << setw(6) << (glued-holes) << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << (glued-holes)/maxbins*100 <<"%" << endl;
    }
    if( glued+sampled > 0 ){
        vout << "      Total area:                 "
             << setw(6) << (glued+sampled) << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << (glued+sampled)/maxbins*100 <<"%" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFEnthalpy::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-enthalpy terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

