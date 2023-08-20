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
#include "Enthalpy.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <ESPrinter.hpp>
#include <iomanip>
#include <EnergyProxyInit.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CEnthalpy)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnthalpy::CEnthalpy(void)
{
    State = 1;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CEnthalpy::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CIntOpts
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

    HEOutputName = Options.GetProgArg(Options.GetNumberOfProgArgs()-1);

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-enthalpy (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(HEOutputName != "-") {
        vout << "# Enthalpy file (out)       : " << HEOutputName << endl;
    } else {
        vout << "# Enthalpy file (out)       : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
        vout << "# Processed realm           : " << Options.GetOptRealm() << endl;
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptMethod() == "raw" ) {
        vout << "# Method                    :     RAW (raw data from  accumulator)" << endl;
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
    if( OutputFile.Open(HEOutputName,"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CEnthalpy::Run(void)
{
// load accumulator
    State = 1;

    vout << endl;
    vout << format("%02d:Loading PMF accumulators ...")%State << endl;
    State++;
    for(int i=0; i < Options.GetNumberOfProgArgs()-1; i++){
        CSmallString name = Options.GetProgArg(i);
        vout << format("** PMF Accumulator #%05d: %s")%(i+1)%string(name) << endl;
        CPMFAccumulatorPtr p_accu(new CPMFAccumulator);
        try {
            p_accu->Load(name);
        } catch(...) {
            CSmallString error;
            error << "unable to load the input PMF accumulator file '" << name << "'";
            ES_ERROR(error);
            return(false);
        }
        Accumulators.push_back(p_accu);
    }
    vout << "   Done" << endl;

    if( Accumulators.size() == 0 ){
        CSmallString error;
        error << "no PMF accumulator was loaded";
        ES_ERROR(error);
        return(false);
    }

// realms
    vout << endl;
    vout << format("%02d:Initializing %s realm ...")%State%Options.GetOptRealm() << endl;
    vout << format(  "   Number of loaded PMF accumulators = %d")%Accumulators.size() << endl;
    State++;
    for(size_t i=0; i < Accumulators.size(); i++){
        CPMFAccumulatorPtr  accu  = Accumulators[i];
        CEnergyProxyPtr     proxy = CEnergyProxyInit::InitProxy(Options.GetOptRealm(),accu);
        proxy->Init(accu);
        EnergyProxies.push_back(proxy);
    }

    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accu->SetNCorr(Options.GetOptNCorr());

// -----------------------------------------------------------------------------
    vout << endl;
    vout << format("%02d:Statistics of input accumulators")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done." << endl;

// -----------------------------------------------------------------------------

    HES = CEnergySurfacePtr(new CEnergySurface);
    HES->Allocate(Accumulators[0]);
    HES->SetSLevel(Options.GetOptSLevel());

    if( Options.GetOptMethod() == "raw" ){
        if( Accumulators.size() > 1 ){
            CSmallString error;
            error << "more than one PMF accumulator was loaded";
            ES_ERROR(error);
            return(false);
        }
        Accumulators[0]->SetNCorr(Options.GetOptNCorr());
        vout << endl;
        vout << format("%02d:Raw absolute enthalpy")%State << endl;
        GetRawEnthalpy();

        if( ! Options.GetOptAbsolute() ){
            AdjustGlobalMin();
        }

        vout << "      SigmaF2   = " << setprecision(5) << HES->GetSigmaF2() << endl;
        vout << "      SigmaF    = " << setprecision(5) << HES->GetSigmaF() << endl;
        if( Options.GetOptWithError() ){
        vout << "      RMSError  = " << setprecision(5) << HES->GetRMSError() << endl;
        vout << "      MaxError  = " << setprecision(5) << HES->GetMaxError() << endl;
        }
        State++;
        vout << "   Done." << endl;
        Accumulators[0]->SetNCorr(1.0);    // to prevent possible errors if used later

    } else if ( Options.GetOptMethod() == "gpr" ){
        vout << endl;
        vout << format("%02d:GPR interpolated enthalpy")%State << endl;
        CSmootherGPR   entgpr;

        entgpr.SetOutputES(HES);
        for(size_t i=0; i < EnergyProxies.size(); i++){
            entgpr.AddInputEnergyProxy(EnergyProxies[i]);
        }

        if( Options.IsOptLoadHyprmsSet() ){
            entgpr.LoadGPRHyprms(Options.GetOptLoadHyprms());
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
        entgpr.SetCalcLogPL(Options.GetOptGPRCalcLogPL());
        entgpr.SetKernel(Options.GetOptGPRKernel());

        if(entgpr.Interpolate(vout) == false) {
            ES_ERROR("unable to interpolate enthalpy");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( entgpr.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

        State++;
        vout << "   Done." << endl;
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    if( Options.GetOptAbsolute() == false ){
        if( ! Options.IsOptGlobalMinSet() ){
            HES->ApplyOffset(Options.GetOptOffset() - HES->GetGlobalMinimumValue());
        } else {
            HES->ApplyOffset(Options.GetOptOffset());
        }
    }

    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            HES->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            HES->AdaptUnsampledToMaxEnergy();
        }
    }

// -----------------------------------------------------------------------------
// print energy surface

    if(PrintHES() == false) {
        ES_ERROR("unable to print enthalpy");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CEnthalpy::AdjustGlobalMin(void)
{
    CSimpleVector<double>  GPos;
    GPos.CreateVector(HES->GetNumOfCVs());

// adjust global minimum
    if( Options.IsOptGlobalMinSet()  ){
        string sspec(Options.GetOptGlobalMin());

        // remove "x" from the string
        replace (sspec.begin(), sspec.end(), 'x' , ' ');

        // parse values of CVs
        stringstream str(sspec);
        for(int i=0; i < HES->GetNumOfCVs(); i++){
            double val;
            str >> val;
            if( ! str ){
                CSmallString error;
                error << "unable to decode CV value for position: " << i+1;
                RUNTIME_ERROR(error);
            }
            GPos[i] = HES->GetCV(i)->GetIntValue(val);
        }

        // GPos.CreateVector(NCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << HES->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < HES->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << HES->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << endl;

        vout << "      Closest bin found at: ";
        // find the closest bin
        CSimpleVector<double>   pos;
        pos.CreateVector(HES->GetNumOfCVs());
        double minv = 0.0;
        int    glb_bin = 0;
        for(int ibin=0; ibin < HES->GetNumOfBins(); ibin++){
            HES->GetPoint(ibin,pos);
            double dist2 = 0.0;
            for(int cv=0; cv < HES->GetNumOfCVs(); cv++){
                dist2 = dist2 + (pos[cv]-GPos[cv])*(pos[cv]-GPos[cv]);
            }
            if( ibin == 0 ){
                minv = dist2;
                glb_bin = 0;
            }
            if( dist2 < minv ){
                minv = dist2;
                glb_bin = ibin;
            }
        }

        HES->GetPoint(glb_bin,GPos);

        vout << setprecision(5) << HES->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < HES->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << HES->GetCV(i)->GetRealValue(GPos[i]);
        }

        double glb_min = HES->GetEnergy(glb_bin);
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(int ibin=0; ibin < HES->GetNumOfBins(); ibin++){
            int samples = HES->GetNumOfSamples(ibin);
            if( samples != 0 ){
                double ene = HES->GetEnergy(ibin);
                HES->SetEnergy(ibin,ene-glb_min);
            }
        }

    } else {
        // search for global minimum
        double glb_min = 0.0;
        for(int ibin=0; ibin < HES->GetNumOfBins(); ibin++){
            int samples = HES->GetNumOfSamples(ibin);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = HES->GetEnergy(ibin);
            if( (ibin == 0) || (glb_min > value) ){
                glb_min = value;
                HES->GetPoint(ibin,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum found at: ";
        vout << setprecision(5) << HES->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < HES->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << HES->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(int ibin=0; ibin < HES->GetNumOfBins(); ibin++){
            int samples = HES->GetNumOfSamples(ibin);
            if( samples != 0 ){
                double ene = HES->GetEnergy(ibin);
                HES->SetEnergy(ibin,ene-glb_min);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CEnthalpy::GetRawEnthalpy(void)
{
    for(int ibin=0; ibin < Accumulators[0]->GetNumOfBins(); ibin++){
        int    nsamples = EnergyProxies[0]->GetNumOfSamples(ibin);
        double ent = EnergyProxies[0]->GetValue(ibin,E_PROXY_VALUE);
        double error = EnergyProxies[0]->GetValue(ibin,E_PROXY_ERROR);
        HES->SetNumOfSamples(ibin,nsamples);
        HES->SetEnergy(ibin,ent);
        HES->SetError(ibin,error);
    }
}

//------------------------------------------------------------------------------

bool CEnthalpy::PrintHES(void)
{
    vout << endl;
    vout << format("%02d:Writing results to file: %s")%State%string(HEOutputName) << endl;

    if( OutputFile.Open(HEOutputName,"w") == false ){
        ES_ERROR("unable to open output file");
        return(false);
    }

    State++;
    CESPrinter printer;

    if(Options.GetOptPrintAll()) {
        printer.SetSampleLimit(0);
    } else {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());

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

    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetPrintedES(HES);

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

void CEnthalpy::WriteHeader(void)
{
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(OutputFile);
        Accumulators[0]->PrintInfo(OutputFile);
    }
}

//------------------------------------------------------------------------------

void CEnthalpy::PrintSampledStat(void)
{
    for(size_t i=0; i < Accumulators.size(); i++){
        CPMFAccumulatorPtr accu = Accumulators[i];
        vout << format("** PMF Accumulator #%05d ...")%(i+1) << endl;
        // calculate sampled area
        double maxbins = accu->GetNumOfBins();
        int    sampled = 0;
        int    limit = 0;
        for(int ibin=0; ibin < accu->GetNumOfBins(); ibin++) {
            if( accu->GetNumOfSamples(ibin) > 0 ) {
                sampled++;
            }
            if( accu->GetNumOfSamples(ibin) > Options.GetOptLimit() ) {
                limit++;
            } else {
                accu->SetNumOfSamples(ibin,0);
            }
        }
        if( maxbins > 0 ){
            vout << " Sampled area: "
                 << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
            vout << " Within limit: "
                 << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%" << endl;
        }
        vout << endl;
        if( accu->CheckCVSInfo(Accumulators[0]) == false ){
            CSmallString error;
            error << "inconsistent dimensions of two PMF accumulators";
            RUNTIME_ERROR(error);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnthalpy::Finalize(void)
{
    // close files if they are own by program
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-enthalpy terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

