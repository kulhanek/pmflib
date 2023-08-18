// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include "GHSEnergyIntegrate.hpp"
#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <IntegratorRFD.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include <iomanip>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <StdIOFile.hpp>
#include <GHSIntegratorGPR0A.hpp>
#include <GHSIntegratorGPR0B.hpp>
#include <ABFProxy_dG.hpp>
#include <ABFProxy_dH.hpp>
#include <ABFProxy_mTdS.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CGHSEnergyIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGHSEnergyIntegrate::CGHSEnergyIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CGHSEnergyIntegrate::Init(int argc,char* argv[])
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
    vout << "# pmf-integrate (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

        vout << "# Input ABF accumulator:                " << Options.GetArgAccuFile() << endl;
        vout << "# Output free energy surface (FES):     " << Options.GetArgFESFile() << endl;
        vout << "# Output enthalpy surface (HES):        " << Options.GetArgHESFile() << endl;
        vout << "# Output entropic energy surface (SES): " << Options.GetArgSESFile() << endl;

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

    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CGHSEnergyIntegrate::Run(void)
{
// load accumulator
    State = 1;

    vout << endl;
    CSmallString name = Options.GetArgAccuFile();
    vout << format("%02d:Loading PMF accumulator: %s")%State%string(name) << endl;
    State++;
    Accumulator = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        Accumulator->Load(name);
    } catch(...) {
        CSmallString error;
        error << "unable to load the input ABF accumulator file '" << name << "'";
        ES_ERROR(error);
        return(false);
    }
    vout << "   Done" << endl;

// realms
    vout << endl;
    vout << format("%02d:Initializing GHS realms")%State  << endl;
    State++;
    vout << format("   ** FES [from ABF dG(x)/dx]") << endl;
    GDerProxy = CABFProxy_dG_Ptr(new CABFProxy_dG());
    GDerProxy->Init(Accumulator);

    FES = CEnergySurfacePtr(new CEnergySurface);
    FES->Allocate(Accumulator);
    FES->SetSLevel(Options.GetOptSLevel());

    vout << format("   ** HES [from ABF dH(x)/dx]") << endl;
    HDerProxy = CABFProxy_dH_Ptr(new CABFProxy_dH());
    HDerProxy->Init(Accumulator);

    HES = CEnergySurfacePtr(new CEnergySurface);
    HES->Allocate(Accumulator);
    HES->SetSLevel(Options.GetOptSLevel());

    vout << format("   ** SES [from ABF -TdS(x)/dx]") << endl;
    SDerProxy = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS());
    SDerProxy->Init(Accumulator);

    SES = CEnergySurfacePtr(new CEnergySurface);
    SES->Allocate(Accumulator);
    SES->SetSLevel(Options.GetOptSLevel());
    vout << "   Done." << endl;

    vout << endl;
    vout << format("%02d:Statistics of the input PMF accumulator")%State << endl;
    State++;
    PrintAccuStat();
    PrintSampledStat();
    vout << "   Done." << endl;

    // test early stage parsing of --globalmin
    CIntegratorRFD  integrator;

    integrator.SetOutputES(FES);
    integrator.SetInputEnergyDerProxy(GDerProxy);

    if( Options.IsOptGlobalMinSet() ){
        integrator.SetGlobalMin(Options.GetOptGlobalMin());
    }

// sampling limit -------------------------------
    vout << endl;
    vout << format("%02d:Preparing PMF accumulators for integration (sampling limit)")%State << endl;
    State++;
    PrepareAccumulatorI();
    if( ! Options.GetOptSkipFFTest() ){
        FloodFillTest();
    }
    PrintAccuStat();
    PrintSampledStat();
    vout << "   Done." << endl;

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:PMF accumulator integration")%State << endl;
    State++;
    if( Integrate() == false ) return(false);
    vout << "   Done." << endl;

// print result ---------------------------------
    vout << endl;
    vout << format("%02d:Writing results")%State << endl;
    State++;
    vout << format("   ** FES [dG(x)]   : %s")%string(Options.GetArgFESFile()) << endl;
    WriteES(FES,Options.GetArgFESFile());
    vout << format("   ** HES [dH(x)]   : %s")%string(Options.GetArgHESFile()) << endl;
    WriteES(HES,Options.GetArgHESFile());
    vout << format("   ** SES [-TdS(x)] : %s")%string(Options.GetArgSESFile()) << endl;
    WriteES(SES,Options.GetArgSESFile());

    vout << "   Done." << endl;

    return(true);
}

//------------------------------------------------------------------------------

void CGHSEnergyIntegrate::WriteES(CEnergySurfacePtr& surf,const CSmallString& name)
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

bool CGHSEnergyIntegrate::Integrate(void)
{

    CGHSIntegratorGPR0A   integrator;

    integrator.SetGDerProxy(GDerProxy);
    integrator.SetHDerProxy(HDerProxy);
    integrator.SetSDerProxy(SDerProxy);

    integrator.SetOutputFES(FES);
    integrator.SetOutputHES(HES);
    integrator.SetOutputSES(SES);

    if( Options.IsOptLoadHyprmsSet() ){
        integrator.LoadGPRHyprms(Options.GetOptLoadHyprms());
    } else {
        integrator.SetSigmaF2(Options.GetOptSigmaF2());
        integrator.SetCoVar(Options.GetOptCoVar());
        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetSigmaN2(Options.GetOptSigmaN2());
    }

    integrator.EnableConstraints(Options.GetOptImposeConstraints());
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

    if( Options.IsOptMFInfoSet() ){
       // integrator.PrepForMFInfo();
    }

    if(integrator.Integrate(vout) == false) {
        ES_ERROR("unable to integrate ABF accumulator");
        return(false);
    }

    if( Options.IsOptMFInfoSet() ){
       // if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero

void CGHSEnergyIntegrate::PrepareAccumulatorI(void)
{
    for(int ibin=0; ibin < Accumulator->GetNumOfBins(); ibin++) {
        // erase datapoints not properly sampled, preserve glueing
        if( (Accumulator->GetNumOfSamples(ibin) >= 0) && (Accumulator->GetNumOfSamples(ibin) <= Options.GetOptLimit()) ) {
            Accumulator->SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CGHSEnergyIntegrate::SyncFESWithAccu(void)
{
    for(int ibin=0; ibin < FES->GetNumOfBins(); ibin++) {
        FES->SetNumOfSamples(ibin,0);
    }

    for(int ibin=0; ibin < Accumulator->GetNumOfBins(); ibin++) {
        int osam = FES->GetNumOfSamples(ibin);
        int nsam = Accumulator->GetNumOfSamples(ibin);
        FES->SetNumOfSamples(ibin,osam+nsam);
    }
}

//------------------------------------------------------------------------------

void CGHSEnergyIntegrate::SyncAccuWithFES(void)
{
    for(int ibin=0; ibin < Accumulator->GetNumOfBins(); ibin++) {
        if( FES->GetNumOfSamples(ibin) <= 0 ) {
            Accumulator->SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CGHSEnergyIntegrate::PrintAccuStat(void)
{
    // calculate sampled area
    double maxbins = Accumulator->GetNumOfBins();
    int    sampled = 0;
    for(int ibin=0; ibin < Accumulator->GetNumOfBins(); ibin++) {
        if( Accumulator->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
    }
    if( maxbins > 0 ){
        vout << "   -- Sampled area:               "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%";
    }
    vout << endl;
}

//------------------------------------------------------------------------------

void CGHSEnergyIntegrate::PrintSampledStat(void)
{
    SyncFESWithAccu();

    // calculate sampled area
    double maxbins = FES->GetNumOfBins();
    int    sampled = 0;
    int    holes = 0;
    int    glued = 0;
    for(int ibin=0; ibin < FES->GetNumOfBins(); ibin++) {
        if( FES->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( FES->GetNumOfSamples(ibin) < 0 ) {
            glued++;
        }
        if( FES->GetNumOfSamples(ibin) == -1 ) {
            holes++;
        }
    }
    if( (maxbins > 0) && (glued != 0) ){
        vout << "   -- Sampled area:               "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
    if( glued > 0 ){
        vout << "   -- All inter/extrapolated area:"
             << setw(6) << glued << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << glued/maxbins*100 <<"%" << endl;
    }
    if( holes > 0 ){
        vout << "   -- Interpolated area:       "
             << setw(6) << holes << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << holes/maxbins*100 <<"%" << endl;
    }
    if( (glued-holes) > 0 ){
        vout << "   -- Extrapolated area:       "
             << setw(6) << (glued-holes) << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << (glued-holes)/maxbins*100 <<"%" << endl;
    }
    if( glued+sampled > 0 ){
        vout << "   -- Total area:                 "
             << setw(6) << (glued+sampled) << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << (glued+sampled)/maxbins*100 <<"%" << endl;
    }
}

//------------------------------------------------------------------------------

void CGHSEnergyIntegrate::FloodFillTest(void)
{
    vout << "   Searching for discontinuous regions ..." << endl;
    int seedid = 1;

    SyncFESWithAccu();

    FFSeeds.CreateVector(FES->GetNumOfBins());
    FFSeeds.SetZero();
    IPos.CreateVector(FES->GetNumOfCVs());
    TPos.CreateVector(FES->GetNumOfCVs());

    double maxbins = FES->GetNumOfBins();
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

    for(int ibin=0; ibin < FES->GetNumOfBins(); ibin++) {
        if( FFSeeds[ibin] != maxseedid ) {
            FES->SetNumOfSamples(ibin,0);
        }
    }

    SyncAccuWithFES();
}

//------------------------------------------------------------------------------

bool CGHSEnergyIntegrate::InstallNewSeed(int seedid,bool unsampled)
{
    for(int ibin=0; ibin < FES->GetNumOfBins(); ibin++) {
        if( unsampled ){
            if( (FFSeeds[ibin] == 0) && ( FES->GetNumOfSamples(ibin) == 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        } else {
            if( (FFSeeds[ibin] == 0) && ( FES->GetNumOfSamples(ibin) != 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

int CGHSEnergyIntegrate::FillSeed(int seedid,bool unsampled)
{
    int newsamples = 0;
    int ndir = 1;
    for(int j=0; j < FES->GetNumOfCVs(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < FES->GetNumOfBins(); ibin++) {
        if( unsampled ){
            if( FES->GetNumOfSamples(ibin) > 0 ) continue; // skip sampled regions
        } else {
            if( FES->GetNumOfSamples(ibin) == 0 ) continue; // skip unsampled regions
        }
        if( FFSeeds[ibin] != seedid ) continue; // skip different regions

        // convert to ipont
        FES->GetIPoint(ibin,IPos);

        // in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = FES->GetGlobalIndex(TPos);
            if( tbin >= 0 ){
                if( FFSeeds[tbin] == 0 ){
                    if( unsampled ){
                        if( FES->GetNumOfSamples(tbin) == 0 ){
                            FFSeeds[tbin] = seedid;
                            newsamples++;
                        }
                    } else {
                        if( FES->GetNumOfSamples(tbin) != 0 ){
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

void CGHSEnergyIntegrate::GetTPoint(CSimpleVector<int>& ipos,int d,CSimpleVector<int>& tpos)
{
    for(int k=FES->GetNumOfCVs()-1; k >= 0; k--) {
        int ibin = d % 3 - 1;
        tpos[k] = ibin + ipos[k];
        d = d / 3;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSEnergyIntegrate::Finalize(void)
{
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

