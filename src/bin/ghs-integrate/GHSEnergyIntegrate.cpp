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
// -------------
#include <GHSIntegratorGPR0C.hpp>
#include <GHSIntegratorGPRcC.hpp>
// -------------
#include <ABFProxy_dG.hpp>
#include <ABFProxy_mTdS.hpp>
#include <PMFProxy_dH.hpp>

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
    Accu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        Accu->Load(name);
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

// -------
    vout << format("   ** FES [from ABF dG(x)/dx]") << endl;
    GDerProxy = CABFProxy_dG_Ptr(new CABFProxy_dG());
    GDerProxy->Init(Accu);

    FES = CEnergySurfacePtr(new CEnergySurface);
    FES->Allocate(Accu);
    FES->SetSLevel(Options.GetOptSLevel());

    if( Options.IsOptGlobalMinSet() ){
        FES->SetGlobalMin(Options.GetOptGlobalMin());
    }

// -------
    vout << format("   ** HES [from dH(x)]") << endl;
    HEneProxy = CPMFProxy_dH_Ptr(new CPMFProxy_dH());
    HEneProxy->Init(Accu);

    HES = CEnergySurfacePtr(new CEnergySurface);
    HES->Allocate(Accu);
    HES->SetSLevel(Options.GetOptSLevel());

// -------
    vout << format("   ** SES [from ABF -TdS(x)/dx]") << endl;
    SDerProxy = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS());
    SDerProxy->Init(Accu);

    SES = CEnergySurfacePtr(new CEnergySurface);
    SES->Allocate(Accu);
    SES->SetSLevel(Options.GetOptSLevel());
    vout << "   Done." << endl;

    vout << endl;
    vout << format("%02d:Statistics of the input PMF accumulator")%State << endl;
    State++;
    PrintAccuStat();
    vout << "   Done." << endl;

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:PMF accumulator integration")%State << endl;
    State++;
    if( Options.GetOptRealm() == "GHS_dH" ) {
        if( Integrate0C() == false ) return(false);
    } else if( Options.GetOptRealm() == "cGHS_dH" ) {
        if( IntegratecC() == false ) return(false);
    } else {
        RUNTIME_ERROR("unsupported realm");
    }

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

bool CGHSEnergyIntegrate::Integrate0C(void)
{

    CGHSIntegratorGPR0C   integrator;

    integrator.SetAccumulator(Accu);

    integrator.SetGDerProxy(GDerProxy);
    integrator.SetHEneProxy(HEneProxy);
    integrator.SetSDerProxy(SDerProxy);

    integrator.SetOutputFES(FES);
    integrator.SetOutputHES(HES);
    integrator.SetOutputSES(SES);

    if( Options.IsOptLoadHyprmsSet() ){
        integrator.LoadGPRHyprms(Options.GetOptLoadHyprms());
    } else {
        integrator.SetSigmaF2(Options.GetOptSigmaF2());
        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetSigmaN2(Options.GetOptSigmaN2());
    }

    integrator.SetIncludeError(Options.GetOptWithError());
    integrator.SetNoEnergy(Options.GetOptNoEnergy());
    integrator.SetUseNumDiff(Options.GetOptGPRNumDiff());

    integrator.SetRCond(Options.GetOptRCond());
    integrator.SetLAMethod(Options.GetOptLAMethod());
    integrator.SetUseInv(Options.GetOptGPRUseInv());
    integrator.SetKernel(Options.GetOptGPRKernel());
    integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());

    if( Options.IsOptMFInfoSet() ){
       integrator.PrepForMFInfo();
    }

    if(integrator.Integrate(vout) == false) {
        ES_ERROR("unable to integrate ABF accumulator");
        return(false);
    }
    vout << "   Done." << endl;

    if( Options.IsOptMFInfoSet() ){
    vout << endl;
    vout << format("%02d:MF Info file: %s")%State%string(Options.GetOptMFInfo()) << endl;
    State++;
        CSmallString mfinfo;
        mfinfo = Options.GetOptMFInfo();
        mfinfo << ".dG_dx";
    vout << format("   ** dG(x)/dx") << endl;
        if( integrator.WriteMFInfo(mfinfo,0) == false ) return(false);
        mfinfo = Options.GetOptMFInfo();
        mfinfo << ".dH";
    vout << format("   ** dH(x)") << endl;
        if( integrator.WriteMFInfo(mfinfo,1) == false ) return(false);
        mfinfo = Options.GetOptMFInfo();
        mfinfo << ".mTdS_dx";
    vout << format("   ** -TdS(x)/dx") << endl;
        if( integrator.WriteMFInfo(mfinfo,2) == false ) return(false);
    }
    vout << "   Done." << endl;

    return(true);
}

//------------------------------------------------------------------------------

bool CGHSEnergyIntegrate::IntegratecC(void)
{

    CGHSIntegratorGPRcC   integrator;

    integrator.SetAccumulator(Accu);

    integrator.SetGDerProxy(GDerProxy);
    integrator.SetHEneProxy(HEneProxy);
    integrator.SetSDerProxy(SDerProxy);

    integrator.SetOutputFES(FES);
    integrator.SetOutputHES(HES);
    integrator.SetOutputSES(SES);

    if( Options.IsOptLoadHyprmsSet() ){
        integrator.LoadGPRHyprms(Options.GetOptLoadHyprms());
    } else {
        integrator.SetSigmaF2(Options.GetOptSigmaF2());
        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetSigmaN2(Options.GetOptSigmaN2());
    }

    integrator.SetIncludeError(Options.GetOptWithError());
    integrator.SetNoEnergy(Options.GetOptNoEnergy());
    integrator.SetUseNumDiff(Options.GetOptGPRNumDiff());

    integrator.SetRCond(Options.GetOptRCond());
    integrator.SetLAMethod(Options.GetOptLAMethod());
    integrator.SetUseInv(Options.GetOptGPRUseInv());
    integrator.SetKernel(Options.GetOptGPRKernel());
    integrator.SetCalcLogPL(Options.GetOptGPRCalcLogPL());

    if( Options.IsOptMFInfoSet() ){
       integrator.PrepForMFInfo();
    }

    if(integrator.Integrate(vout) == false) {
        ES_ERROR("unable to integrate ABF accumulator");
        return(false);
    }
    vout << "   Done." << endl;

    if( Options.IsOptMFInfoSet() ){
    vout << endl;
    vout << format("%02d:MF Info file: %s")%State%string(Options.GetOptMFInfo()) << endl;
    State++;
        CSmallString mfinfo;
        mfinfo = Options.GetOptMFInfo();
        mfinfo << ".dG_dx";
    vout << format("   ** dG(x)/dx") << endl;
        if( integrator.WriteMFInfo(mfinfo,0) == false ) return(false);
        mfinfo = Options.GetOptMFInfo();
        mfinfo << ".dH";
    vout << format("   ** dH(x)") << endl;
        if( integrator.WriteMFInfo(mfinfo,1) == false ) return(false);
        mfinfo = Options.GetOptMFInfo();
        mfinfo << ".mTdS_dx";
    vout << format("   ** -TdS(x)/dx") << endl;
        if( integrator.WriteMFInfo(mfinfo,2) == false ) return(false);
    }
    vout << "   Done." << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSEnergyIntegrate::PrintAccuStat(void)
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

