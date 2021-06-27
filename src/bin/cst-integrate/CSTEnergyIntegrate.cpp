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

#include "CSTEnergyIntegrate.hpp"
#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <IntegratorRFD.hpp>
#include <IntegratorRBF.hpp>
#include <IntegratorGPR.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include <iomanip>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <CSTProxy_dG.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CCSTEnergyIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCSTEnergyIntegrate::CCSTEnergyIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCSTEnergyIntegrate::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CCSTIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

    if( (Options.GetNumberOfProgArgs() != 2) && (Options.GetNumberOfProgArgs() != 3) ){
        ES_ERROR("two or three arguments are expected");
        return(SO_OPTS_ERROR);
    }

    CSTAccuName = Options.GetProgArg(0);
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
    vout << "# abf-integrate (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if( CSTAccuName != "-") {
        vout << "# CST accu file (in)    : " << CSTAccuName << endl;
    } else {
        vout << "# CST accu file (in)    : - (standard input)" << endl;
    }
    if( FEOutputName != "-") {
        vout << "# Free energy file (out): " << FEOutputName << endl;
    } else {
        vout << "# Free energy file (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
        if(Options.GetOptMethod() == "rfd" ) {
            vout << "# Integration method    : RFD (reverse finite differences via csparse)" << endl;
        } else if( Options.GetOptMethod() == "rbf" ){
            vout << "# Integration method    : RBF (radial basis functions)" << endl;
        } else if( Options.GetOptMethod() == "gpr" ) {
            vout << "# Integration method    : GPR (gaussian process)" << endl;
        } else {
            INVALID_ARGUMENT("method - not implemented");
        }
        if( Options.GetOptWithError() ) {
        vout << "# Integrated domains    : force+error" << endl;
        } else {
        vout << "# Integrated domains    : force only" << endl;
        }
// ---------------------------------------
        vout << "# Linear algebra        : " << Options.GetOptLAMethod() << endl;
    if ( Options.GetOptMethod() == "rbf" ){
        if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "default") ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        }
   } else if ( Options.GetOptMethod() == "gpr"  ) {
        if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        }
    }
// ---------------------------------------
        vout << "# ------------------------------------------------" << endl;
    if( Options.GetOptMethod() == "rfd" ){
        vout << "# FD number of points   : " << Options.GetOptFDPoints() << endl;
        vout << "# Periodicity           : " << bool_to_str(Options.GetOptPeriodicity()) << endl;
    } else if ( Options.GetOptMethod() == "rbf" ){
        vout << "# Reduction factor rfac : " << Options.GetOptRFac() << endl;
        vout << "# Width factor wfac     : " << Options.GetOptWFac() << endl;
        vout << "# RBF overhang          : " << Options.GetOptOverhang() << endl;
    } else if ( Options.GetOptMethod() == "gpr"  ) {
            vout << "# Split NCorr mode      : " << bool_to_str(Options.GetOptSplitNCorr()) << endl;
        if( Options.IsOptLoadHyprmsSet() ){
            vout << "# GPR hyperprms file    : " << Options.GetOptLoadHyprms() << endl;
            // actual values are printed in detailed output from integrator
        } else {
            vout << "# SigmaF2               : " << setprecision(3) << Options.GetOptSigmaF2() << endl;
            vout << "# NCorr                 : " << setprecision(3) << Options.GetOptNCorr() << endl;
            vout << "# Width factor wfac     : " << Options.GetOptWFac() << endl;
        }
    } else {
        ES_ERROR("not implemented method");
        return(SO_USER_ERROR);
    }

    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptEnergyLimit() == -1) {
        vout << "# Energy limit          : not applied" << endl;
    } else {
        vout << "# Energy limit          : " << Options.GetOptEnergyLimit() << endl;
    }
        vout << "# Skip last energy limit: " << bool_to_str(Options.GetOptSkipLastEnergyLimit()) << endl;

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
    if( InputFile.Open(CSTAccuName,"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CCSTEnergyIntegrate::Run(void)
{
// load accumulator
    int state = 1;

// -----------------------------------------------------------------------------
// setup accu, energy proxy, and output FES
    Accu        = CPMFAccumulatorPtr(new CPMFAccumulator);
    FES         = CEnergySurfacePtr(new CEnergySurface);
    DerProxy    = CCSTProxy_dG_Ptr(new CCSTProxy_dG);

    vout << endl;
    vout << format("%02d:Loading CST accumulator: %s")%state%string(CSTAccuName) << endl;
    state++;
    try {
        Accu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input CST accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    // print CVS info
    Accu->PrintInfo(vout);
    DerProxy->Init(Accu);

    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accu->SetNCorr(Options.GetOptNCorr());

    FES->Allocate(Accu);
    FES->SetSLevel(Options.GetOptSLevel());

// reduced FES options
    if( Options.IsOptKeepCVsSet() ){
        DecodeEList(Options.GetOptKeepCVs(),KeepCVs,"--keepcvs");
    }

    vout << endl;
    vout << format("%02d:Statistics of input CST accumulator")%state << endl;
    state++;
    PrintSampledStat();
    vout << "   Done." << endl;

    if( (Options.GetOptMethod() == "gpr") || (Options.GetOptMethod() == "rbf") ){
        // test early stage parsing of --globalmin
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }
    }

// energy limit --------------------------------

    if( Options.GetOptEnergyLimit() > 0.0 ){
        vout << endl;
        vout << format("%02d:CST accumulator integration (%s) for energy limit")%state%string(Options.GetOptEcutMethod()) << endl;
        state++;
        if( IntegrateForEcut() == false ) return(false);

        vout << endl;
        vout << format("%02d:Preparing CST accumulator for integration (energy limit)")%state << endl;
        state++;
        PrepareAccumulatorII();
        PrintSampledStat();
        vout << "   Done." << endl;

        FES->Clear();
    }

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:CST accumulator integration (%s)")%state%string(Options.GetOptMethod()) << endl;
    state++;
    if( Integrate() == false ) return(false);
    vout << "   Done." << endl;

 // apply offset
    if( ! Options.IsOptGlobalMinSet() ){
        FES->ApplyOffset(Options.GetOptOffset() - FES->GetGlobalMinimumValue());
    } else {
        FES->ApplyOffset(Options.GetOptOffset());
    }

// final energy limit --------------------------------

    if( (Options.GetOptEnergyLimit() > 0.0) && (Options.GetOptSkipLastEnergyLimit() == false) ){
        vout << endl;
        vout << format("%02d:Cleaning FES (energy limit)")%state << endl;
        state++;
        PrepareAccumulatorII();
        PrintSampledStat();
        SyncFESWithACCU();
        vout << "   Done." << endl;
    }

// post-processing
    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            FES->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            FES->AdaptUnsampledToMaxEnergy();
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
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    if(Options.GetOptPrintAll()) {
        printer.SetSampleLimit(0);
    }

    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
    printer.SetPrintedES(FES);

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
        } else {
            INVALID_ARGUMENT("output format - not implemented");
        }

        // print all
        printer.SetSampleLimit(0);
        printer.SetIncludeError(Options.GetOptWithError());
        printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
        printer.SetPrintedES(FES);

        try {
            printer.Print(OutputFile);
        } catch(...) {
            ES_ERROR("unable to save the output free energy file");
            return(false);
        }
        vout << "   Done." << endl;
    }

// save accumulator if requested
    if( Options.GetOptSaveCST() != NULL ){
        vout << endl;
        vout << format("%02d:Saving CST accumulator to : %s")%state%string(Options.GetOptSaveCST()) << endl;
        state++;
        try {
            Accu->Save(Options.GetOptSaveCST());
        } catch(...) {
            ES_ERROR("unable to save the CST accumulator file");
            return(false);
        }
        vout << "   Done." << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CCSTEnergyIntegrate::WriteHeader()
{
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(OutputFile);
        Accu->PrintInfo(OutputFile);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CCSTEnergyIntegrate::IntegrateForEcut(void)
{
    if(Options.GetOptEcutMethod() == "rfd" ) {
        CIntegratorRFD   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "lu" ){
                // nothing to do - LU is default
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "rbf" ){
        CIntegratorRBF   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac(Options.GetOptRFac());
        integrator.SetOverhang(Options.GetOptOverhang());

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            integrator.SetLLSMethod(Options.GetOptLAMethod());
        }

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "gpr" ){
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        integrator.SetSplitNCorr(Options.GetOptSplitNCorr());

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetUseNumDiff(Options.GetOptGPRNumDiff());
        integrator.SetIncludeError(false);

        integrator.SetRCond(Options.GetOptRCond());
        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            integrator.SetLAMethod(Options.GetOptLAMethod());
        }
        integrator.SetUseInv(Options.GetOptGPRUseInv());
        integrator.SetKernel(Options.GetOptGPRKernel());

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CCSTEnergyIntegrate::Integrate(void)
{
    if(Options.GetOptMethod() == "rfd" ) {
        CIntegratorRFD   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        if( Options.GetOptLAMethod() == "lu" ){
            // nothing to do - LU is default
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetUseOldRFDMode(Options.GetOptUseOldRFD());

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "rbf" ){
        CIntegratorRBF   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac(Options.GetOptRFac());
        integrator.SetOverhang(Options.GetOptOverhang());

        integrator.SetLLSMethod(Options.GetOptLAMethod());

        if( Options.GetOptUseRealGlobalMin() == false ){
            if( Options.IsOptGlobalMinSet() ){
                integrator.SetGlobalMin(Options.GetOptGlobalMin());
            }
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

    } else if( Options.GetOptMethod() == "gpr" ){
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        integrator.SetSplitNCorr(Options.GetOptSplitNCorr());

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetFastError(!Options.GetOptGPRNoFastError());
        integrator.SetIncludeError(Options.GetOptWithError());
        integrator.SetNoEnergy(Options.GetOptNoEnergy());
        integrator.SetUseNumDiff(Options.GetOptGPRNumDiff());

        if( Options.GetOptUseRealGlobalMin() == false ){
            if( Options.IsOptGlobalMinSet() ){
                integrator.SetGlobalMin(Options.GetOptGlobalMin());
            }
        }

        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetLAMethod(Options.GetOptLAMethod());
        integrator.SetUseInv(Options.GetOptGPRUseInv());
        integrator.SetKernel(Options.GetOptGPRKernel());
        integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());
        integrator.SetUseZeroPoint(Options.GetOptGPRIncludeZPE());

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

        GPos = integrator.GetGlobalMin();

    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CCSTEnergyIntegrate::ReduceFES(void)
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
    if( nrcvs == (size_t)Accu->GetNumOfCVs() ){
        vout << "   No reduction specified, skipping ..." << endl;
        return(true);
    }
    if( nrcvs == 0 ){
        vout << "   Too large reduction specified, skipping ..." << endl;
        return(true);
    }

    vout << format("   Temperature : %.1f K")%(Accu->GetTemperature()) << endl;

    // FIXME
    CEnergySurfacePtr reducedFES;

    if( (Options.GetOptMethod() == "gpr") && Options.GetOptWithError() ){
        // need to run another integration
        CIntegratorGPR   integrator;

        // FES is destroyed during reduction by CIntegratorGPR, thus use some temp version
        CEnergySurfacePtr tmp_FES = CEnergySurfacePtr(new CEnergySurface);
        tmp_FES->Allocate(Accu);

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(tmp_FES);

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
        integrator.SetUseNumDiff(Options.GetOptGPRNumDiff());
        integrator.SetGlobalMin(GPos);
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetLAMethod(Options.GetOptLAMethod());
        integrator.SetUseInv(Options.GetOptGPRUseInv());
        integrator.SetKernel(Options.GetOptGPRKernel());
        integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());
        integrator.SetUseZeroPoint(true);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate CST accumulator");
            return(false);
        }
        reducedFES = integrator.ReduceFES(KeepCVs);
        if( reducedFES == NULL ) {
            ES_ERROR("unable to reduce FES");
            return(false);
        }

    } else {
        reducedFES = FES->ReduceFES(KeepCVs);
        if( reducedFES == NULL ) {
            ES_ERROR("unable to reduce FES");
            return(false);
        }
    }

// post-processing
    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            reducedFES->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            reducedFES->AdaptUnsampledToMaxEnergy();
        }
    }

    CESPrinter printer;

    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    printer.SetSampleLimit(0);
    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
    printer.SetPrintedES(reducedFES);

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

void CCSTEnergyIntegrate::PrintGPRHyprms(FILE* p_fout)
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

void CCSTEnergyIntegrate::LoadGPRHyprms(CIntegratorGPR& gpr)
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
                for(int i=0; i < Accu->GetNumOfCVs();i++){
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
//    a) erase data points with large energy

void CCSTEnergyIntegrate::PrepareAccumulatorII(void)
{
    // filter by energy
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( DerProxy->GetNumOfSamples(ibin) != 0 ) {
            // consider only properly sampled data points
            if( FES->GetEnergy(ibin) > Options.GetOptEnergyLimit() ){
                // erase datapoints with too large energy
                DerProxy->SetNumOfSamples(ibin,0);
            }
            if( Options.GetOptEraseNegativeEnergy() ){
                if( FES->GetEnergy(ibin) < 0 ){
                    // erase datapoints with negative energy
                    DerProxy->SetNumOfSamples(ibin,0);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void CCSTEnergyIntegrate::SyncFESWithACCU(void)
{
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( DerProxy->GetNumOfSamples(ibin) == 0 ) {
            FES->SetEnergy(ibin,0.0);
            FES->SetEnergy(ibin,0.0);
            FES->SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CCSTEnergyIntegrate::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = Accu->GetNumOfBins();
    int    sampled = 0;
    int    holes = 0;
    int    glued = 0;
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( DerProxy->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( DerProxy->GetNumOfSamples(ibin) < 0 ) {
            glued++;
        }
        if( DerProxy->GetNumOfSamples(ibin) == -1 ) {
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

void CCSTEnergyIntegrate::DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname)
{
    int ncvs = Accu->GetNumOfCVs();

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

void CCSTEnergyIntegrate::Finalize(void)
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
    vout << "# abf-integrate terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

