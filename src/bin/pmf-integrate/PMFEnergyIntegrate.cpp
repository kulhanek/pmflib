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

#include "PMFEnergyIntegrate.hpp"
#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <IntegratorRFD.hpp>
#include <IntegratorRBF.hpp>
#include <IntegratorGPR.hpp>
#include <SmootherGPR.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include <iomanip>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <ABFProxy_dG.hpp>
#include <ABFProxy_mTdS.hpp>
#include <CSTProxy_dG.hpp>
#include <CSTProxy_mTdS.hpp>
#include <CSTProxy_MTC.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CPMFEnergyIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFEnergyIntegrate::CPMFEnergyIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFEnergyIntegrate::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

    if( (Options.GetNumberOfProgArgs() != 2) && (Options.GetNumberOfProgArgs() != 3) ){
        ES_ERROR("two or three arguments are expected");
        return(SO_OPTS_ERROR);
    }

    PMFAccuName = Options.GetProgArg(0);
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
    vout << "# pmf-integrate (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if( PMFAccuName != "-") {
        vout << "# ABF accu file (in)    : " << PMFAccuName << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if( FEOutputName != "-") {
        vout << "# Free energy file (out): " << FEOutputName << endl;
    } else {
        vout << "# Free energy file (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
        vout << "# Integrated realm      : " << Options.GetOptRealm() << endl;
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
    if(Options.GetOptLimit() == 0) {
        vout << "# Sampling limit        : all bins will be taken into account" << endl;
    } else {
        vout << "# Sampling limit        : " << Options.GetOptLimit() << endl;
    }
    if( (Options.GetOptEcutMethod() == "gpr") || (Options.GetOptEcutMethod() == "rbf") ){
    if(Options.GetOptMFMaxZScore() == -1) {
        vout << "# Max MF error Z-score  : not applied" << endl;
    } else {
        vout << "# Max MF error Z-score  : " << Options.GetOptMFMaxZScore() << endl;
        vout << "# Number of MF Z-tests  : " << Options.GetOptMFZTestPasses() << endl;
    }
    }
        vout << "# Glueing FES factor    : " << Options.GetOptGlueingFactor() << endl;
        vout << "# Glue holes on FES     : " << bool_to_str(Options.GetOptGlueHoles()) << endl;
    if(Options.GetOptEnergyLimit() == -1) {
        vout << "# Energy limit          : not applied" << endl;
    } else {
        vout << "# Energy limit          : " << Options.GetOptEnergyLimit() << endl;
    }
        vout << "# Skip last energy limit: " << bool_to_str(Options.GetOptSkipLastEnergyLimit()) << endl;
        vout << "# Skip flood fill test  : " << bool_to_str(Options.GetOptNoHeader()) << endl;

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
    if( InputFile.Open(PMFAccuName,"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CPMFEnergyIntegrate::Run(void)
{
// load accumulator
    int state = 1;

// -----------------------------------------------------------------------------
// setup accu, energy proxy, and output FES
    Accu        = CPMFAccumulatorPtr(new CPMFAccumulator);
    FES         = CEnergySurfacePtr(new CEnergySurface);

    vout << endl;
    vout << format("%02d:Loading ABF accumulator: %s")%state%string(PMFAccuName) << endl;
    state++;
    try {
        Accu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    // print CVS info
    Accu->PrintInfo(vout);


// init energyder proxy
    if( Options.GetOptRealm() == "dG" ){
        if( CABFProxy_dG::IsCompatible(Accu) ){
            DerProxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        } else if (CCSTProxy_dG::IsCompatible(Accu) ) {
            DerProxy    = CCSTProxy_dG_Ptr(new CCSTProxy_dG);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm: " <<  Options.GetOptRealm();
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "dG_p" ) {
        CABFProxy_dG_Ptr proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        proxy->SetType(ABF_MICF_POT);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "dG_k" ) {
        CABFProxy_dG_Ptr proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        proxy->SetType(ABF_MICF_KIN);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS" ) {
        if( CABFProxy_mTdS::IsCompatible(Accu) ){
            DerProxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        } else if (CCSTProxy_mTdS::IsCompatible(Accu) ) {
            DerProxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm: " <<  Options.GetOptRealm();
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_PP" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_PP);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_PK" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_PK);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_PR" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_PR);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_KP" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_KP);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_KK" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_KK);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_KR" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_KR);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_HP" ) {
        CCSTProxy_mTdS_Ptr proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        proxy->SetType(CST_C11HP);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_HK" ) {
        CCSTProxy_mTdS_Ptr proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        proxy->SetType(CST_C11HK);
        DerProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "-TdS_HR" ) {
        CCSTProxy_mTdS_Ptr proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        proxy->SetType(CST_C11HR);
        DerProxy = proxy;
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << Options.GetOptRealm();
        RUNTIME_ERROR(error);
    }

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
    vout << format("%02d:Statistics of input ABF accumulator")%state << endl;
    state++;
    PrintSampledStat();
    vout << "   Done." << endl;

    if( (Options.GetOptMethod() == "rfd") || (Options.GetOptMethod() == "gpr") || (Options.GetOptMethod() == "rbf") ){
        // test early stage parsing of --globalmin
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }
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

    if( Options.GetOptMFMaxZScore() > 0.0 ){
        for(int i=1; i <= Options.GetOptMFZTestPasses(); i++ ){
            vout << endl;
            vout << format("%02d:ABF accumulator integration (%s) for mean force error Z-score test #%d")%state%string(Options.GetOptEcutMethod())%i << endl;
            state++;
            if( IntegrateForMFZScore(i) == false ) return(false);

            vout << endl;
            vout << format("%02d:Preparing ABF accumulator for integration (mean force error z-score test #%d")%state%i << ")"<< endl;
            state++;
            PrepareAccumulatorI();
            if( ! Options.GetOptSkipFFTest() ){
                FloodFillTest();
            }
            PrintSampledStat();
            vout << "   Done." << endl;

            FES->Clear();
        }
    }

// glue fes ------------------------------------
    if( Options.GetOptGlueHoles() ){
        vout << endl;
        vout << format("%02d:Preparing ABF accumulator for integration (glue holes on FES)")%state << endl;
        state++;
        GlueHoles();
        PrintSampledStat();
        vout << "   Done." << endl;
        FES->Clear();
    }

    if( Options.GetOptGlueingFactor() > 0 ){
        vout << endl;
        vout << format("%02d:Preparing ABF accumulator for integration (glueing FES)")%state << endl;
        state++;
        vout << "   Searching for border regions in close vicinity of sampled areas ..." << endl;
        int tg = 0;
        for(int i=1; i <= Options.GetOptGlueingFactor(); i++ ){
            tg += GlueingFES(i);
        }
        vout << "   -- Total glued bins: " << tg << endl;
        PrintSampledStat();
        vout << "   Done." << endl;

        FES->Clear();
    }

// energy limit --------------------------------

    if( Options.GetOptEnergyLimit() > 0.0 ){
        vout << endl;
        vout << format("%02d:ABF accumulator integration (%s) for energy limit")%state%string(Options.GetOptEcutMethod()) << endl;
        state++;
        if( IntegrateForEcut() == false ) return(false);

        vout << endl;
        vout << format("%02d:Preparing ABF accumulator for integration (energy limit)")%state << endl;
        state++;
        PrepareAccumulatorII();
        if( ! Options.GetOptSkipFFTest() ){
            FloodFillTest();
        }
        PrintSampledStat();
        vout << "   Done." << endl;

        FES->Clear();
    }

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:ABF accumulator integration (%s)")%state%string(Options.GetOptMethod()) << endl;
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
        if( ! Options.GetOptSkipFFTest() ){
            FloodFillTest();
        }
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
    } else {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());
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
        printer.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());
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
    if( Options.GetOptSaveABF() != NULL ){
        vout << endl;
        vout << format("%02d:Saving ABF accumulator to : %s")%state%string(Options.GetOptSaveABF()) << endl;
        state++;
        try {
            Accu->Save(Options.GetOptSaveABF());
        } catch(...) {
            ES_ERROR("unable to save the ABF accumulator file");
            return(false);
        }
        vout << "   Done." << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CPMFEnergyIntegrate::WriteHeader()
{
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(OutputFile);
        Accu->PrintInfo(OutputFile);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFEnergyIntegrate::IntegrateForMFZScore(int pass)
{
    if(Options.GetOptEcutMethod() == "rfd" ) {
        ES_ERROR("illegal combination: --emethod=rfd and --maxzscore");
        return(false);
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

        integrator.SetNoEnergy(true);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            CSmallString mfinfo = Options.GetOptMFInfo();
            mfinfo << ".mflimit" << pass;
            if( integrator.WriteMFInfo(mfinfo) == false ) return(false);
        }

        // apply mean force limit
        integrator.FilterByMFZScore(Options.GetOptMFMaxZScore(),vout);

    } else if( Options.GetOptEcutMethod() == "gpr" ){
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

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

        integrator.SetNoEnergy(true);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            CSmallString mfinfo = Options.GetOptMFInfo();
            mfinfo << ".mflimit" << pass;
            if( integrator.WriteMFInfo(mfinfo) == false ) return(false);
        }

        // apply mean force limit
        integrator.FilterByMFZScore(Options.GetOptMFMaxZScore(),vout);

    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    // add metric tensor correction
    if( Options.GetOptRealm() == "dG" ){
        AddMTCorr();
    } else if ( Options.GetOptRealm() == "-TdS" ) {
        AddMTCorr();
    } else {
        // do not add MTC
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFEnergyIntegrate::IntegrateForEcut(void)
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

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
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
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            integrator.SetLLSMethod(Options.GetOptLAMethod());
        }

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "gpr" ){
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetUseNumDiff(Options.GetOptGPRNumDiff());
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());
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
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    // add metric tensor correction
    if( Options.GetOptRealm() == "dG" ){
        AddMTCorr();
    } else if ( Options.GetOptRealm() == "-TdS" ) {
        AddMTCorr();
    } else {
        // do not add MTC
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFEnergyIntegrate::Integrate(void)
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

        if( Options.GetOptUseRealGlobalMin() == false ){
            if( Options.IsOptGlobalMinSet() ){
                integrator.SetGlobalMin(Options.GetOptGlobalMin());
            }
        }

        integrator.SetUseOldRFDMode(Options.GetOptUseOldRFD());

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        GPos = integrator.GetGlobalMin();

    } else if( Options.GetOptMethod() == "rbf" ){
        CIntegratorRBF   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac(Options.GetOptRFac());
        integrator.SetOverhang(Options.GetOptOverhang());
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());

        integrator.SetLLSMethod(Options.GetOptLAMethod());

        if( Options.GetOptUseRealGlobalMin() == false ){
            if( Options.IsOptGlobalMinSet() ){
                integrator.SetGlobalMin(Options.GetOptGlobalMin());
            }
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

        GPos = integrator.GetGlobalMin();

    } else if( Options.GetOptMethod() == "gpr" ){
        CIntegratorGPR   integrator;

        integrator.SetInputEnergyDerProxy(DerProxy);
        integrator.SetOutputES(FES);

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
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());

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
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

        GPos = integrator.GetGlobalMin();

    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    // add metric tensor correction
    if( Options.GetOptRealm() == "dG" ){
        AddMTCorr();
    } else if ( Options.GetOptRealm() == "-TdS" ) {
        AddMTCorr();
    } else {
        // do not add MTC
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFEnergyIntegrate::ReduceFES(void)
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
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptGlueHoles()||Options.GetOptIncludeGluedRegions());
        integrator.SetGlobalMin(GPos);
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetLAMethod(Options.GetOptLAMethod());
        integrator.SetUseInv(Options.GetOptGPRUseInv());
        integrator.SetKernel(Options.GetOptGPRKernel());
        integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());
        integrator.SetUseZeroPoint(true);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
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

void CPMFEnergyIntegrate::PrintGPRHyprms(FILE* p_fout)
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

void CPMFEnergyIntegrate::LoadGPRHyprms(CIntegratorGPR& gpr)
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

void CPMFEnergyIntegrate::LoadGPRHyprms(CSmootherGPR& gpr)
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

void CPMFEnergyIntegrate::PrepareAccumulatorI(void)
{
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        // erase datapoints not properly sampled, preserve glueing
        if( (DerProxy->GetNumOfSamples(ibin) >= 0) && (DerProxy->GetNumOfSamples(ibin) <= Options.GetOptLimit()) ) {
            DerProxy->SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

// this part performs following tasks:
//    a) erase data points with large energy

void CPMFEnergyIntegrate::PrepareAccumulatorII(void)
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

void CPMFEnergyIntegrate::SyncFESWithACCU(void)
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

void CPMFEnergyIntegrate::PrintSampledStat(void)
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

//------------------------------------------------------------------------------

void CPMFEnergyIntegrate::FloodFillTest(void)
{
    vout << "   Searching for discontinuous regions ..." << endl;
    int seedid = 1;

    FFSeeds.CreateVector(Accu->GetNumOfBins());
    FFSeeds.SetZero();
    IPos.CreateVector(Accu->GetNumOfCVs());
    TPos.CreateVector(Accu->GetNumOfCVs());

    double maxbins = Accu->GetNumOfBins();
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

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] != maxseedid ) {
            DerProxy->SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

bool CPMFEnergyIntegrate::InstallNewSeed(int seedid,bool unsampled)
{
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( unsampled ){
            if( (FFSeeds[ibin] == 0) && ( DerProxy->GetNumOfSamples(ibin) == 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        } else {
            if( (FFSeeds[ibin] == 0) && ( DerProxy->GetNumOfSamples(ibin) != 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

int CPMFEnergyIntegrate::FillSeed(int seedid,bool unsampled)
{
    int newsamples = 0;
    int ndir = 1;
    for(int j=0; j < Accu->GetNumOfCVs(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( unsampled ){
            if( DerProxy->GetNumOfSamples(ibin) > 0 ) continue; // skip sampled regions
        } else {
            if( DerProxy->GetNumOfSamples(ibin) == 0 ) continue; // skip unsampled regions
        }
        if( FFSeeds[ibin] != seedid ) continue; // skip different regions

        // convert to ipont
        Accu->GetIPoint(ibin,IPos);

        // in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accu->GetGlobalIndex(TPos);
            if( tbin >= 0 ){
                if( FFSeeds[tbin] == 0 ){
                    if( unsampled ){
                        if( DerProxy->GetNumOfSamples(tbin) == 0 ){
                            FFSeeds[tbin] = seedid;
                            newsamples++;
                        }
                    } else {
                        if( DerProxy->GetNumOfSamples(tbin) != 0 ){
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

void CPMFEnergyIntegrate::GetTPoint(CSimpleVector<int>& ipos,int d,CSimpleVector<int>& tpos)
{
    for(int k=Accu->GetNumOfCVs()-1; k >= 0; k--) {
        int ibin = d % 3 - 1;
        tpos[k] = ibin + ipos[k];
        d = d / 3;
    }
}

//------------------------------------------------------------------------------

int CPMFEnergyIntegrate::GlueingFES(int factor)
{
    IPos.CreateVector(Accu->GetNumOfCVs());
    TPos.CreateVector(Accu->GetNumOfCVs());

    int ndir = 1;
    for(int j=0; j < Accu->GetNumOfCVs(); j++){
        ndir *= 3;
    }

    vout << "   Gluing FES: factor = " << factor;

    int glued = 0;

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( DerProxy->GetNumOfSamples(ibin) != 0 ) continue; // skip glued or sampled bins

        // convert to ipont
        Accu->GetIPoint(ibin,IPos);

        // is sampled or glued region in close vicinty?

        // in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accu->GetGlobalIndex(TPos);
            if( tbin >= 0 ){
                if( factor == 1 ){
                    if( DerProxy->GetNumOfSamples(tbin) > 0 ){
                        DerProxy->SetNumOfSamples(ibin,-(factor+1));
                        glued++;
                        break;
                    }
                } else {
                    if( DerProxy->GetNumOfSamples(tbin) == -factor ){
                        DerProxy->SetNumOfSamples(ibin,-(factor+1));
                        glued++;
                        break;
                    }
                }
            }
        }
    }

    vout << ", glued bins = " << glued << endl;

    return(glued);
}

//------------------------------------------------------------------------------

void CPMFEnergyIntegrate::GlueHoles(void)
{
    vout << "   Searching for holes on FES ..." << endl;
    int seedid = 1;

    FFSeeds.CreateVector(Accu->GetNumOfBins());
    FFSeeds.SetZero();
    IPos.CreateVector(Accu->GetNumOfCVs());
    TPos.CreateVector(Accu->GetNumOfCVs());

    double maxbins = Accu->GetNumOfBins();
    int    numofholes = 0;

    int sampled = SeedSampled(seedid);

    vout << "   Region: " << setw(6) << seedid << " - area:         "
         << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"% sampled." << endl;

    seedid++;

    int tg = 0;

    while( InstallNewSeed(seedid,true) ){
        int sampled = 1;    // for initial seed set by InstallNewSeed
        int newsamples = 0;

        while( (newsamples = FillSeed(seedid,true)) > 0 ){
            sampled += newsamples;
        }

        // detect type of area
        bool hole = false;
        if( IsHole(seedid) ){
            MarkAsHole(seedid);
            numofholes++;
            hole = true;
        }

        if( maxbins > 0 ){
            vout << "   Region: " << setw(6) << seedid << " - area:         "
                 << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%";
            if( hole ){
                vout << " hole - glued." << endl;
                tg += sampled;
            } else {
                vout << " edge." << endl;
            }
        }
        seedid++;
    }

    // print stat
    vout << "   -- Number of holes      : " <<  numofholes << endl;
    vout << "   -- Number of glued bins : " <<  tg << endl;
}

//------------------------------------------------------------------------------

int CPMFEnergyIntegrate::SeedSampled(int seedid)
{
    int sampled = 0;

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( DerProxy->GetNumOfSamples(ibin) > 0 ) {
            FFSeeds[ibin] = seedid;
            sampled++;
        }
    }

    return(sampled);
}

//------------------------------------------------------------------------------

bool CPMFEnergyIntegrate::IsHole(int seedid)
{
    IPos.CreateVector(Accu->GetNumOfCVs());
    TPos.CreateVector(Accu->GetNumOfCVs());

    int ndir = 1;
    for(int j=0; j < Accu->GetNumOfCVs(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] != seedid ) continue;

        // convert to ipont
        Accu->GetIPoint(ibin,IPos);

        // test neighbouring of the point in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accu->GetGlobalIndex(TPos);
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

void CPMFEnergyIntegrate::MarkAsHole(int seedid)
{
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] == seedid ){
            DerProxy->SetNumOfSamples(ibin,-1);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFEnergyIntegrate::DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname)
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

//------------------------------------------------------------------------------

void CPMFEnergyIntegrate::AddMTCorr(void)
{
    if( CCSTProxy_MTC::IsCompatible(Accu) == false ) return;

    vout << "   Adding metric tensor correction ..." << endl;

    CCSTProxy_MTC_Ptr MTCProxy    = CCSTProxy_MTC_Ptr(new CCSTProxy_MTC);
    MTCProxy->Init(Accu);

    for(int i=0; i < Accu->GetNumOfBins(); i++){
        FES->SetNumOfSamples(i, MTCProxy->GetNumOfSamples(i) );
        double f = FES->GetEnergy(i);
        FES->SetEnergy(i, f + MTCProxy->GetValue(i,E_PROXY_VALUE) );
    }

    if( Options.IsOptGlobalMinSet() ){

        CSmootherGPR        mtcgpr;
        CEnergySurfacePtr   mtc = CEnergySurfacePtr(new CEnergySurface);
        mtc->Allocate(Accu);

        mtcgpr.SetInputEnergyProxy(MTCProxy);
        mtcgpr.SetOutputES(mtc);

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(mtcgpr);
        } else {
            mtcgpr.SetSigmaF2(Options.GetOptSigmaF2());
            mtcgpr.SetNCorr(Options.GetOptNCorr());
            mtcgpr.SetWFac(Options.GetOptWFac());
        }

        mtcgpr.SetGlobalMin(Options.GetOptGlobalMin());

        mtcgpr.SetRCond(Options.GetOptRCond());
        mtcgpr.SetLAMethod(Options.GetOptLAMethod());
        mtcgpr.SetKernel(Options.GetOptGPRKernel());

        if(mtcgpr.Interpolate(vout) == false) {
            RUNTIME_ERROR("unable to interpolate MTC");
        }

        vout << "      Adjusting global minimum for metric tensor correction ..." << endl;
        for(int i=0; i < Accu->GetNumOfBins(); i++){
            FES->SetNumOfSamples(i, MTCProxy->GetNumOfSamples(i) );
            double f = FES->GetEnergy(i);
            FES->SetEnergy(i, f - mtcgpr.GetGlobalMinValue() );
        }
    } else {
        // search for global minimum
        CSimpleVector<double> gpos;
        gpos.CreateVector(Accu->GetNumOfCVs());
        bool   first = true;
        double glb_min = 0.0;
        for(int i=0; i < FES->GetNumOfBins(); i++){
            int samples = FES->GetNumOfSamples(i);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = FES->GetEnergy(i);
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                Accu->GetPoint(i,gpos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum found at: ";
        vout << setprecision(5) << Accu->GetCV(0)->GetRealValue(gpos[0]);
        for(int i=1; i < Accu->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << Accu->GetCV(i)->GetRealValue(gpos[i]);
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(int i=0; i < FES->GetNumOfBins(); i++){
            if( FES->GetNumOfSamples(i) != 0 ) {
                double value = FES->GetEnergy(i);
                FES->SetEnergy(i,value-glb_min);
            }
        }
    }

        vout << "      SigmaF2   = " << setprecision(5) << FES->GetSigmaF2() << endl;
    if( Options.GetOptIncludeGluedRegions() ){
        vout << "      SigmaF2 (including glued bins) = " << setprecision(5) << FES->GetSigmaF2(true) << endl;
    }
        vout << "      SigmaF    = " << setprecision(5) << FES->GetSigmaF() << endl;

    if( Options.GetOptWithError() ){
        vout << "      RMSError  = " << setprecision(5) << FES->GetRMSError() << endl;
        vout << "      MaxError  = " << setprecision(5) << FES->GetMaxError() << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFEnergyIntegrate::Finalize(void)
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
    vout << "# pmf-integrate terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

