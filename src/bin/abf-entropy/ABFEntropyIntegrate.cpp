// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <ABFIntegratorGPR.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include "ABFEntropyIntegrate.hpp"
#include <iomanip>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFEntropyIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFEntropyIntegrate::CABFEntropyIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFEntropyIntegrate::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFEntropyIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

    if( (Options.GetNumberOfProgArgs() != 2) && (Options.GetNumberOfProgArgs() != 3) ){
        ES_ERROR("two or three arguments are expected");
        return(SO_OPTS_ERROR);
    }

    ABFAccuName = Options.GetProgArg(0);
    FEOutputName = Options.GetProgArg(1);
    if( Options.GetNumberOfProgArgs() == 3 ) {
        // optional
        FullFEOutputName = Options.GetProgArg(2);
    }

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    StartTime.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-entropy (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if( ABFAccuName != "-") {
        vout << "# ABF accu file (in)    : " << ABFAccuName << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if( FEOutputName != "-") {
        vout << "# Free energy file (out): " << FEOutputName << endl;
    } else {
        vout << "# Free energy file (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    vout << "# Integration method    : GPR (gaussian process)" << endl;

    if( Options.GetOptWithError() ) {
    vout << "# Integrated domains    : force+error" << endl;
    } else {
    vout << "# Integrated domains    : force only" << endl;
    }
// ---------------------------------------
    vout << "# Linear algebra        : " << Options.GetOptLAMethod() << endl;
    if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
    vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
    }
// ---------------------------------------
    vout << "# ------------------------------------------------" << endl;
        vout << "# Split NCorr mode      : " << bool_to_str(Options.GetOptSplitNCorr()) << endl;
    if( Options.IsOptLoadHyprmsSet() ){
        vout << "# GPR hyperprms file    : " << Options.GetOptLoadHyprms() << endl;
        // actual values are printed in detailed output from integrator
    } else {
        vout << "# SigmaF2               : " << setprecision(3) << Options.GetOptSigmaF2() << endl;
        vout << "# NCorr                 : " << setprecision(3) << Options.GetOptNCorr() << endl;
        vout << "# Width factor wfac     : " << Options.GetOptWFac() << endl;
    }

    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Sampling limit        : all bins will be taken into account" << endl;
    } else {
        vout << "# Sampling limit        : " << Options.GetOptLimit() << endl;
    }

    vout << "# ------------------------------------------------" << endl;

    if( Options.IsOptGlobalMinSet() ){
    vout << "# Global FES minimum    : " << Options.GetOptGlobalMin() << endl;
    } else {
    vout << "# Global FES minimum    : -auto-" << endl;
    }
    vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# Include bin statuses  : " << bool_to_str(Options.GetOptIncludeBinStat()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;
    if( Options.IsOptKeepCVsSet() ){
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << "# Keep CVs              : " << Options.GetOptKeepCVs() << endl;
    vout << "# Reduced FES file      : " << Options.GetOptReducedFES() << endl;
    }

    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(ABFAccuName,"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFEntropyIntegrate::Run(void)
{
// load accumulator
    int state = 1;

    vout << endl;
    vout << format("%02d:Loading ABF accumulator: %s")%state%string(ABFAccuName) << endl;
    state++;
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    // print CVS info
    Accumulator.PrintCVSInfo(vout);
    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accumulator.SetNCorr(Options.GetOptNCorr());
    FES.Allocate(&Accumulator);
    FES.SetSLevel(Options.GetOptSLevel());

// reduced FES options
    if( Options.IsOptKeepCVsSet() ){
        DecodeEList(Options.GetOptKeepCVs(),KeepCVs,"--keepcvs");
    }

    vout << endl;
    vout << format("%02d:Statistics of input ABF accumulator")%state << endl;
    state++;
    PrintSampledStat();
    vout << "   Done." << endl;

    // test early stage parsing of --globalmin
    CABFIntegratorGPR   integrator;

    integrator.SetInputABFAccumulator(&Accumulator);
    integrator.SetOutputFESurface(&FES);

    if( Options.IsOptGlobalMinSet() ){
        integrator.SetGlobalMin(Options.GetOptGlobalMin());
    }

// sampling limit -------------------------------
    vout << endl;
    vout << format("%02d:Preparing ABF accumulator for integration (sampling limit)")%state << endl;
    state++;
    PrepareAccumulatorI();
    if( ! Options.GetOptSkipFFTest() ){
        FloodFillTest();
    }
    PrintSampledStat();
    vout << "   Done." << endl;

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:ABF accumulator integration (GPR)")%state << endl;
    state++;
    if( Integrate() == false ) return(false);
    vout << "   Done." << endl;

 // apply offset
    if( ! Options.IsOptGlobalMinSet() ){
        FES.ApplyOffset(Options.GetOptOffset() - FES.GetGlobalMinimumValue());
    } else {
        FES.ApplyOffset(Options.GetOptOffset());
    }

// post-processing
    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            FES.AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            FES.AdaptUnsampledToMaxEnergy();
        }
    }

// reduce FES ------------------------------
    if( Options.IsOptReducedFESSet() ){
        vout << endl;
        vout << format("%02d:Reducing FES by statistical reweighting")%state << endl;
        state++;
        if( ReduceFES() == false ) return(false);
        vout << "   Done." << endl;
    }

// print result ---------------------------------
    vout << endl;
    vout << format("%02d:Writing results to file: %s")%state%string(FEOutputName) << endl;

    if( OutputFile.Open(FEOutputName,"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    state++;
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
    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
    printer.SetPrintedES(&FES);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output free energy file");
        return(false);
    }
    vout << "   Done." << endl;

    if( Options.GetNumberOfProgArgs() == 3 ){
        vout << endl;
        vout << format("%02d:Writing results to file: %s (full version, --printall)")%state%string(FullFEOutputName) << endl;
        state++;

        if( OutputFile.Open(FullFEOutputName,"w") == false ){
            ES_ERROR("unable to open output file");
            return(SO_USER_ERROR);
        }

        WriteHeader();

        CESPrinter printer;

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

        // print all
        printer.SetSampleLimit(0);
        printer.SetIncludeError(Options.GetOptWithError());
        printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
        printer.SetPrintedES(&FES);

        try {
            printer.Print(OutputFile);
        } catch(...) {
            ES_ERROR("unable to save the output free energy file");
            return(false);
        }
        vout << "   Done." << endl;
    }

// save accumulator if requested
    if( Options.GetOptSaveABF() != NULL ){
        vout << endl;
        vout << format("%02d:Saving ABF accumulator to : %s")%state%string(Options.GetOptSaveABF()) << endl;
        state++;
        try {
            Accumulator.Save(Options.GetOptSaveABF());
        } catch(...) {
            ES_ERROR("unable to save the ABF accumulator file");
            return(false);
        }
        vout << "   Done." << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::WriteHeader()
{
    // print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        fprintf(OutputFile,"# PMFLib version        : %s\n",LibBuildVersion_PMF);
        fprintf(OutputFile,"# data integrated by    : GPR (gaussian process)\n");
        if( Options.GetOptWithError() ) {
        fprintf(OutputFile,"# Integrated domains    : force+error\n");
        } else {
        fprintf(OutputFile,"# Integrated domains    : force only\n");
        }
        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
            fprintf(OutputFile,"# Linear algebra        : %s\n", (const char*)Options.GetOptLAMethod());
            if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
            fprintf(OutputFile,"# SVD rcond             : %5.4e\n",Options.GetOptRCond());
            }
        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        if( Options.IsOptLoadHyprmsSet() ){
            fprintf(OutputFile,"# GPR hyperprms file    : %s\n", (const char*)Options.GetOptLoadHyprms());
            PrintGPRHyprms(OutputFile);
        } else {
            fprintf(OutputFile,"# SigmaF2               : %10.4f\n", Options.GetOptSigmaF2());
            fprintf(OutputFile,"# NCorr                 : %s\n", (const char*)Options.GetOptNCorr());
            fprintf(OutputFile,"# Width factor wfac     : %s\n", (const char*)Options.GetOptWFac());
            fprintf(OutputFile,"# Split NCorr mode      : %s\n", bool_to_str(Options.GetOptSplitNCorr()));

        }

        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        if( Options.IsOptGlobalMinSet() ){
        fprintf(OutputFile,"# Global FES minimum    : %s\n",(const char*)Options.GetOptGlobalMin());
        } else {
        fprintf(OutputFile,"# Global FES minimum    : -auto-\n");
        }
        fprintf(OutputFile,"# Integration offset    : %5.3f\n", Options.GetOptOffset());
        fprintf(OutputFile,"# Include bin statuses  : %s\n",bool_to_str(Options.GetOptIncludeBinStat()));
        fprintf(OutputFile,"# Number of coordinates : %d\n",Accumulator.GetNumOfCVs());
        fprintf(OutputFile,"# Total number of bins  : %d\n",Accumulator.GetNumOfBins());
    }
}

//------------------------------------------------------------------------------

bool CABFEntropyIntegrate::Integrate(void)
{
    CABFIntegratorGPR   integrator;

    integrator.SetInputABFAccumulator(&Accumulator);
    integrator.SetOutputFESurface(&FES);
    integrator.SetIntegratedRealm(EGPR_REALM_TDS);

    integrator.SetSplitNCorr(Options.GetOptSplitNCorr());

    if( Options.IsOptLoadHyprmsSet() ){
        LoadGPRHyprms(integrator);
    } else {
        integrator.SetSigmaF2(Options.GetOptSigmaF2());
        integrator.SetNCorr(Options.GetOptNCorr());
        integrator.SetWFac(Options.GetOptWFac());
    }

    integrator.SetIncludeError(Options.GetOptWithError());
    integrator.SetNoEnergy(Options.GetOptNoEnergy());

    if( Options.GetOptUseRealGlobalMin() == false ){
        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }
    }

    integrator.SetRCond(Options.GetOptRCond());
    integrator.SetLAMethod(Options.GetOptLAMethod());
    integrator.SetKernel(Options.GetOptGPRKernel());
    integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());

    if(integrator.Integrate(vout) == false) {
        ES_ERROR("unable to integrate ABF accumulator");
        return(false);
    }

    if( Options.IsOptMFInfoSet() ){
        if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
    }

    GPos = integrator.GetGlobalMin();

    return(true);
}

//------------------------------------------------------------------------------

bool CABFEntropyIntegrate::ReduceFES(void)
{
    vout << format("   Reduced FES : %s")%string(Options.GetOptReducedFES()) << endl;

    size_t nrcvs = 0;
           vout << "   Kept CVs    : ";
    for(size_t i=0; i < KeepCVs.size(); i++){
        if( KeepCVs[i] ){
            vout << "T";
            nrcvs++;
        } else {
            vout << "F";
        }
        if( (i+1) < KeepCVs.size() ) vout << "x";
    }
    vout << endl;
    if( nrcvs == (size_t)Accumulator.GetNumOfCVs() ){
        vout << "   No reduction specified, skipping ..." << endl;
        return(true);
    }
    if( nrcvs == 0 ){
        vout << "   Too large reduction specified, skipping ..." << endl;
        return(true);
    }

    vout << format("   Temperature : %.1f K")%Options.GetOptTemperature() << endl;

    CEnergySurface reducedFES;

    if( Options.GetOptWithError() ){
        // need to run another integration
        CABFIntegratorGPR   integrator;

        // FES is destroyed during reduction by CABFIntegratorGPR, thus use some temp version
        CEnergySurface tmp_FES;
        tmp_FES.Allocate(&Accumulator);

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&tmp_FES);
        integrator.SetIntegratedRealm(EGPR_REALM_TDS);

        integrator.SetSplitNCorr(Options.GetOptSplitNCorr());

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetFastError(true);
        integrator.SetIncludeError(true);
        integrator.SetNoEnergy(false);
        integrator.SetGlobalMin(GPos);
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetLAMethod(Options.GetOptLAMethod());
        integrator.SetKernel(Options.GetOptGPRKernel());
        integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
        if( integrator.ReduceFES(KeepCVs,Options.GetOptTemperature(),&reducedFES) == false ){
            ES_ERROR("unable to reduce FES");
            return(false);
        }

    } else {
        if( FES.ReduceFES(KeepCVs,Options.GetOptTemperature(),&reducedFES) == false ){
            ES_ERROR("unable to reduce FES");
            return(false);
        }
    }

// post-processing
    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            reducedFES.AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            reducedFES.AdaptUnsampledToMaxEnergy();
        }
    }

    CESPrinter printer;

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

    printer.SetSampleLimit(0);
    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
    printer.SetPrintedES(&reducedFES);

    try {
        printer.Print(Options.GetOptReducedFES());
    } catch(...) {
        ES_ERROR("unable to save the reduced free energy file");
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFEntropyIntegrate::PrintGPRHyprms(FILE* p_fout)
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
            fprintf(p_fout,"# SigmaF2               : %10.4f\n",value);
        } else if( key == "NCorr" ){
            fprintf(p_fout,"# NCorr                 : %10.4f\n",value);
        } else if( key.find("NCorr#") != string::npos ) {
            fprintf(p_fout,"# %-7s               : %10.4f\n",key.c_str(),value);
        } else if( key.find("WFac#") != string::npos ) {
            fprintf(p_fout,"# %-7s               : %10.4f\n",key.c_str(),value);
        } else {
            CSmallString error;
            error << "GPR hyperparameters file, unrecognized key: " << key.c_str();
            RUNTIME_ERROR(error);
        }
    }
}

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::LoadGPRHyprms(CABFIntegratorGPR& gpr)
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
            gpr.SetNCorr(0,value);
            if( Options.GetOptSplitNCorr() ){
                vout << "   >> WARNING: Expanding NCorr due to --splitncorr" << endl;
                for(int i=0; i < Accumulator.GetNumOfCVs();i++){
                    gpr.SetNCorr(i,value);
                }
            }
        } else if( key.find("NCorr#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string sncorr;
            int    cvind;
            kstr >> sncorr >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode ncorr key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( Options.GetOptSplitNCorr() ){
                gpr.SetNCorr(cvind,value);
            } else {
                if( cvind > 0 ){
                    vout << "   >> WARNING: NCorr with index higher than one ignored (no --splitncorr)!" << endl;
                } else {
                    gpr.SetNCorr(cvind,value);
                }
            }

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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero

void CABFEntropyIntegrate::PrepareAccumulatorI(void)
{
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        // erase datapoints not properly sampled, preserve glueing
        if( (Accumulator.GetNumOfSamples(ibin) >= 0) && (Accumulator.GetNumOfSamples(ibin) <= Options.GetOptLimit()) ) {
            Accumulator.SetNumberOfABFSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::SyncFESWithACCU(void)
{
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( Accumulator.GetNumOfSamples(ibin) == 0 ) {
            FES.SetEnergy(ibin,0.0);
            FES.SetEnergy(ibin,0.0);
            FES.SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::PrintSampledStat(void)
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

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::FloodFillTest(void)
{
    vout << "   Searching for discontinuous regions ..." << endl;
    int seedid = 1;

    FFSeeds.CreateVector(Accumulator.GetNumOfBins());
    FFSeeds.SetZero();
    IPos.CreateVector(Accumulator.GetNumOfCVs());
    TPos.CreateVector(Accumulator.GetNumOfCVs());

    double maxbins = Accumulator.GetNumOfBins();
    int    maxseedid = 0;
    int    maxsampled = 0;
    bool   first = true;

    while( InstallNewSeed(seedid,false) ){
        int sampled = 1;    // for initial seed set by InstallNewSeed
        int newsamples = 0;

        while( (newsamples = FillSeed(seedid,false)) > 0 ){
            sampled += newsamples;
        }

        if( maxbins > 0 ){
            vout << "   Region: " << setw(6) << seedid << " - sampled area: "
                 << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
        }

        if( first || (maxsampled < sampled) ){
            first = false;
            maxsampled = sampled;
            maxseedid = seedid;
        }

        seedid++;
    }
    seedid--;

    // quit if one or none region
    if( seedid <= 1 ){
        vout << "   -- All is continuous." << endl;
        return;
    }

        vout << "   -- Clearing all except region: " << maxseedid <<  endl;

    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] != maxseedid ) {
            Accumulator.SetNumberOfABFSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

bool CABFEntropyIntegrate::InstallNewSeed(int seedid,bool unsampled)
{
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( unsampled ){
            if( (FFSeeds[ibin] == 0) && ( Accumulator.GetNumOfSamples(ibin) == 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        } else {
            if( (FFSeeds[ibin] == 0) && ( Accumulator.GetNumOfSamples(ibin) != 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

int CABFEntropyIntegrate::FillSeed(int seedid,bool unsampled)
{
    int newsamples = 0;
    int ndir = 1;
    for(int j=0; j < Accumulator.GetNumOfCVs(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( unsampled ){
            if( Accumulator.GetNumOfSamples(ibin) > 0 ) continue; // skip sampled regions
        } else {
            if( Accumulator.GetNumOfSamples(ibin) == 0 ) continue; // skip unsampled regions
        }
        if( FFSeeds[ibin] != seedid ) continue; // skip different regions

        // convert to ipont
        Accumulator.GetIPoint(ibin,IPos);

        // in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accumulator.GetGlobalIndex(TPos);
            if( tbin >= 0 ){
                if( FFSeeds[tbin] == 0 ){
                    if( unsampled ){
                        if( Accumulator.GetNumOfSamples(tbin) == 0 ){
                            FFSeeds[tbin] = seedid;
                            newsamples++;
                        }
                    } else {
                        if( Accumulator.GetNumOfSamples(tbin) != 0 ){
                            FFSeeds[tbin] = seedid;
                            newsamples++;
                        }
                    }
                }
            }
        }
    }

    return(newsamples);
}

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::GetTPoint(CSimpleVector<int>& ipos,int d,CSimpleVector<int>& tpos)
{
    for(int k=Accumulator.GetNumOfCVs()-1; k >= 0; k--) {
        int ibin = d % 3 - 1;
        tpos[k] = ibin + ipos[k];
        d = d / 3;
    }
}

//------------------------------------------------------------------------------

int CABFEntropyIntegrate::SeedSampled(int seedid)
{
    int sampled = 0;

    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( Accumulator.GetNumOfSamples(ibin) > 0 ) {
            FFSeeds[ibin] = seedid;
            sampled++;
        }
    }

    return(sampled);
}

//------------------------------------------------------------------------------

bool CABFEntropyIntegrate::IsHole(int seedid)
{
    IPos.CreateVector(Accumulator.GetNumOfCVs());
    TPos.CreateVector(Accumulator.GetNumOfCVs());

    int ndir = 1;
    for(int j=0; j < Accumulator.GetNumOfCVs(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] != seedid ) continue;

        // convert to ipont
        Accumulator.GetIPoint(ibin,IPos);

        // test neighbouring of the point in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accumulator.GetGlobalIndex(TPos);
            if( tbin < 0 ){
                // outside of ABF accumulator - it is not a hole
                return(false);
            }
        }
    }

    // it is a hole - all points are confined in sampled regions.
    return(true);
}

//------------------------------------------------------------------------------

void CABFEntropyIntegrate::MarkAsHole(int seedid)
{
    for(int ibin=0; ibin < Accumulator.GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] == seedid ){
            Accumulator.SetNumberOfABFSamples(ibin,-1);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFEntropyIntegrate::DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname)
{
    int ncvs = Accumulator.GetNumOfCVs();

    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( (int)slist.size() > ncvs ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << ncvs << ")";
        RUNTIME_ERROR(error);
    }

    elist.resize(ncvs);

    // parse values
    bool last_st = false;
    for(int i=0; i < (int)slist.size(); i++){
        stringstream str(slist[i]);
        char letter;
        str >> letter;
        if( ! str ){
            CSmallString error;
            error << "unable to decode value for " << optionname << " at position: " << i+1;
            RUNTIME_ERROR(error);
        }
        if( (letter == 'T') || (letter == 't') ){
            last_st = true;
        } else {
            last_st = false;
        }
        elist[i] = last_st;
    }

    // pad the rest with the last value
    for(int i=slist.size(); i < ncvs; i++){
        elist[i] = last_st;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFEntropyIntegrate::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    CSmallTime dur;
    dur = dt - StartTime;

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-entropy terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

