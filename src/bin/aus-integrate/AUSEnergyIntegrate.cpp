// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "AUSEnergyIntegrate.hpp"
#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <IntegratorRFD.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include <iomanip>
#include <algorithm>
#include <boost/format.hpp>
//#include <boost/algorithm/string/split.hpp>
//#include <boost/algorithm/string/classification.hpp>
#include <StdIOFile.hpp>
// -------------
#include <GPREngineAUSInit.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
//using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CAUSEnergyIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAUSEnergyIntegrate::CAUSEnergyIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CAUSEnergyIntegrate::Init(int argc,char* argv[])
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

    StartTime.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# aus-integrate (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

        vout << "# Input ABF accumulator:                " << Options.GetArgAccuFile() << endl;
        vout << "# Output free energy surface (FEN):     " << Options.GetArgFENFile() << endl;
        vout << "# Output internal energy surface (INT): " << Options.GetArgINTFile() << endl;
        vout << "# Output entropic energy surface (TDS): " << Options.GetArgTDSFile() << endl;

    vout << "# ------------------------------------------------" << endl;
        if( Options.GetOptWithError() ) {
        vout << "# Include errors        : yes" << endl;
        } else {
        vout << "# Include errors        : no" << endl;
        }
// ---------------------------------------
        vout << "# Linear algebra        : " << Options.GetOptLAMethod() << endl;
        if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        }
// ---------------------------------------
        vout << "# ------------------------------------------------" << endl;
        if( Options.IsOptLoadHyprmsSet() ){
            vout << "# GPR hyperprms file    : " << Options.GetOptLoadHyprms() << endl;
            // actual values are printed in detailed output from integrator
        } else {
            vout << "# SigmaF2               : " << setprecision(3) << Options.GetOptSigmaF2() << endl;
            vout << "# Width factor wfac     : " << Options.GetOptWFac() << endl;
            vout << "# SigmaN2               : " << Options.GetOptSigmaN2() << endl;
        }

    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Sampling limit        : all bins will be taken into account" << endl;
    } else {
        vout << "# Sampling limit        : " << Options.GetOptLimit() << endl;
    }
    vout << "# Balance res. errors   : " << bool_to_str(Options.GetOptBalanceResiduals()) << endl;
    vout << "# ------------------------------------------------" << endl;

    if( Options.IsOptGlobalMinSet() ){
    vout << "# Global FEN minimum    : " << Options.GetOptGlobalMin() << endl;
    } else {
    vout << "# Global FEN minimum    : -auto-" << endl;
    }
    vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# Output FEN format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# Include bin statuses  : " << bool_to_str(Options.GetOptIncludeBinStat()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;

    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CAUSEnergyIntegrate::Run(void)
{
// load accumulator
    State = 1;

    vout << endl;
    CSmallString name = Options.GetArgAccuFile();
    vout << format("%02d:Loading PMF accumulator: %s")%State%string(name) << endl;
    State++;
    Accu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        Accu->Load(name);
    } catch(...) {
        CSmallString error;
        error << "unable to load the input PMF accumulator file '" << name << "'";
        ES_ERROR(error);
        return(false);
    }
    vout << "   Done" << endl;

// realms
    vout << endl;
    vout << format("%02d:Initializing AUS realm and output ENE surfaces")%State  << endl;
    State++;

    AUSEngine = CGPREngineAUSInit::InitEngine(Options.GetOptRealm(),Accu);

// -------
    vout << format("   ** FEN surface") << endl;

    FEN = CEnergySurfacePtr(new CEnergySurface);
    FEN->Allocate(Accu);
    FEN->SetSLevel(Options.GetOptSLevel());

    if( Options.IsOptGlobalMinSet() ){
        FEN->SetGlobalMin(Options.GetOptGlobalMin());
    }

// -------
    vout << format("   ** INT surface") << endl;

    INT = CEnergySurfacePtr(new CEnergySurface);
    INT->Allocate(Accu);
    INT->SetSLevel(Options.GetOptSLevel());

// -------
    vout << format("   ** TDS surface") << endl;

    TDS = CEnergySurfacePtr(new CEnergySurface);
    TDS->Allocate(Accu);
    TDS->SetSLevel(Options.GetOptSLevel());

// -------
    if( Options.GetOptResidualsFile() != NULL ){
    vout << format("   ** RES surface") << endl;

    RES = CEnergySurfacePtr(new CEnergySurface);
    RES->Allocate(Accu);
    RES->SetSLevel(Options.GetOptSLevel());
    }

    vout << "   Done." << endl;

// -------
    vout << endl;
    vout << format("%02d:Statistics of the input PMF accumulator")%State << endl;
    State++;
    PrintAccuStat();
    vout << "   Done." << endl;

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:PMF accumulator processing")%State << endl;
    State++;
    RunAUSEngine();

// print result ---------------------------------
    vout << endl;
    vout << format("%02d:Writing results")%State << endl;
    State++;
    vout << format("   ** FEN [dA(x)]   : %s")%string(Options.GetArgFENFile()) << endl;
    WriteES(FEN,Options.GetArgFENFile());
    vout << format("   ** INT [dU(x)]   : %s")%string(Options.GetArgINTFile()) << endl;
    WriteES(INT,Options.GetArgINTFile());
    vout << format("   ** TDS [-TdS(x)] : %s")%string(Options.GetArgTDSFile()) << endl;
    WriteES(TDS,Options.GetArgTDSFile());

    if( Options.GetOptResidualsFile() != NULL ){
    vout << format("   ** RES (dA(x) - [dU(x)-TdS(x)]) : %s")%string(Options.GetOptResidualsFile()) << endl;
    WriteES(RES,Options.GetOptResidualsFile());
    }

    vout << "   Done." << endl;

    return(true);
}

//------------------------------------------------------------------------------

void CAUSEnergyIntegrate::WriteES(CEnergySurfacePtr& surf,const CSmallString& name)
{
 // apply offset
    if( ! Options.IsOptGlobalMinSet() ){
        surf->ApplyOffset(Options.GetOptOffset() - surf->GetGlobalMinimumValue());
    } else {
        surf->ApplyOffset(Options.GetOptOffset());
    }

// post-processing
    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            surf->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            surf->AdaptUnsampledToMaxEnergy();
        }
    }

// write result
    CStdIOFile out;

    if( out.Open(name,"w") == false ){
        RUNTIME_ERROR("unable to open output file");
    }

    CESPrinter printer;

    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(out);
    }

    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    printer.SetSampleLimit(Options.GetOptLimit());
    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
    printer.SetPrintedES(surf);

    try {
        printer.Print(out);
    } catch(...) {
        RUNTIME_ERROR("unable to save the output energy surface");
    }

    out.Close();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAUSEnergyIntegrate::RunAUSEngine(void)
{
    AUSEngine->SetAccumulator(Accu);

    AUSEngine->SetOutputFEN(FEN);
    AUSEngine->SetOutputINT(INT);
    AUSEngine->SetOutputTDS(TDS);

    if( Options.GetOptResidualsFile() != NULL ){
        AUSEngine->SetOutputRES(RES);
    }

    // it must be here - it can be redefined in LoadGPRHyprms
    AUSEngine->SetKernel(Options.GetOptGPRKernel());
    AUSEngine->UseFirstKernelDerivatives(Options.GetOptUseFDKernel());

    if( Options.IsOptLoadHyprmsSet() ){
        AUSEngine->LoadGPRHyprms(Options.GetOptLoadHyprms());
    } else {
        AUSEngine->SetSigmaF2(Options.GetOptSigmaF2());
        AUSEngine->SetWFac(Options.GetOptWFac());
        AUSEngine->SetSigmaN2(Options.GetOptSigmaN2());
    }

    AUSEngine->SetIncludeError(Options.GetOptWithError());
    AUSEngine->SetNoEnergy(Options.GetOptNoEnergy());
    AUSEngine->SetBalanceResiduals(Options.GetOptBalanceResiduals());
    AUSEngine->SetUseNumDiff(Options.GetOptGPRNumDiff());

    AUSEngine->SetRCond(Options.GetOptRCond());
    AUSEngine->SetLAMethod(Options.GetOptLAMethod());
    AUSEngine->SetUseInv(Options.GetOptGPRUseInv());
    AUSEngine->SetCalcLogPL(Options.GetOptGPRCalcLogPL());

    if( Options.IsOptMFInfoSet() ){
       AUSEngine->PrepForMFInfo();
    }

    if(AUSEngine->RunGPR(vout) == false) {
        ES_ERROR("unable to integrate ABF accumulator");
        return(false);
    }
    vout << "   Done." << endl;

    if( Options.IsOptMFInfoSet() ){
    vout << endl;
    vout << format("%02d:MF Info file: %s")%State%string(Options.GetOptMFInfo()) << endl;
    State++;
        for(int task=0; task < AUSEngine->GetNumOfTasks(); task++ ){
            vout << format("   ** GPR task: %d")%(task+1) << endl;
            CSmallString mfinfo;
            mfinfo = Options.GetOptMFInfo();
            mfinfo << ".t" << (task+1);
            if( AUSEngine->WriteMFInfo(mfinfo,task) == false ) return(false);
        }
    vout << "   Done." << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAUSEnergyIntegrate::PrintAccuStat(void)
{
    // calculate sampled area
    double maxbins = Accu->GetNumOfBins();
    int    sampled = 0;
    int    limit = 0;
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( Accu->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( Accu->GetNumOfSamples(ibin) > Options.GetOptLimit() ) {
            limit++;
        } else {
            Accu->SetNumOfSamples(ibin,0);
        }
    }
    if( maxbins > 0 ){
        vout << " Sampled area: "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" ;
        vout << " ... Within limit: "
             << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%";
    }
    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAUSEnergyIntegrate::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    CSmallTime dur;
    dur = dt - StartTime;

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# aus-integrate terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

