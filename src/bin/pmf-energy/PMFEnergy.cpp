// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "PMFEnergy.hpp"
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

MAIN_ENTRY(CPMFEnergy)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFEnergy::CPMFEnergy(void)
{
    State = 1;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFEnergy::Init(int argc,char* argv[])
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

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-energy (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

        vout << "# PMF accumulator (in)      : " << Options.GetArgAccuFile() << endl;

    if( Options.GetArgENEFile() != "-") {
        vout << "# Energy file (out)         : " << Options.GetArgENEFile() << endl;
    } else {
        vout << "# Energy file (out)         : - (standard output)" << endl;
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
    if( OutputFile.Open(Options.GetArgENEFile(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFEnergy::Run(void)
{
// load accumulator
    State = 1;

    vout << endl;
    vout << format("%02d:Loading PMF accumulator ...")%State << endl;
    State++;
    vout << format("   ** Name: %s")%string(Options.GetArgAccuFile()) << endl;
    Accu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        Accu->Load(Options.GetArgAccuFile());
    } catch(...) {
        CSmallString error;
        error << "unable to load the input PMF accumulator file '" << Options.GetArgAccuFile() << "'";
        ES_ERROR(error);
        return(false);
    }
    vout << "   Done" << endl;

// realms
    vout << endl;
    vout << format("%02d:Initializing %s realm ...")%State%Options.GetOptRealm() << endl;
    State++;
    EneProxy = CEnergyProxyInit::InitProxy(Options.GetOptRealm(),Accu);
    EneProxy->Init(Accu);
    vout << format(  "   %s")%EneProxy->GetFullDescription() << endl;

    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accu->SetNCorr(Options.GetOptNCorr());

// -----------------------------------------------------------------------------
    vout << endl;
    vout << format("%02d:Statistics of input accumulator")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done." << endl;

// -----------------------------------------------------------------------------

    ENE = CEnergySurfacePtr(new CEnergySurface);
    ENE->Allocate(Accu);
    ENE->SetSLevel(Options.GetOptSLevel());

    if( Options.IsOptGlobalMinSet() ){
        ENE->SetGlobalMin(Options.GetOptGlobalMin());
    }

    if( Options.GetOptMethod() == "raw" ){
        vout << endl;
        vout << format("%02d:Raw absolute energy")%State << endl;
        GetRawEnthalpy();

        if( ! Options.GetOptAbsolute() ){
            AdjustGlobalMin();
        }

        vout << "      SigmaF2   = " << setprecision(5) << ENE->GetSigmaF2() << endl;
        vout << "      SigmaF    = " << setprecision(5) << ENE->GetSigmaF() << endl;
        if( Options.GetOptWithError() ){
        vout << "      RMSError  = " << setprecision(5) << ENE->GetRMSError() << endl;
        vout << "      MaxError  = " << setprecision(5) << ENE->GetMaxError() << endl;
        }
        State++;
        vout << "   Done." << endl;

    } else if ( Options.GetOptMethod() == "gpr" ){
        vout << endl;
        vout << format("%02d:GPR interpolated energy")%State << endl;
        CSmootherGPR   entgpr;

        entgpr.SetOutputES(ENE);
        entgpr.SetInputEnergyProxy(EneProxy);

        // these two lines must be here - it can be overwritten in LoadGPRHyprms
        entgpr.SetKernel(Options.GetOptGPRKernel());
        entgpr.UseFirstKernelDerivatives(Options.GetOptUseFDKernel());

        if( Options.IsOptLoadHyprmsSet() ){
            entgpr.LoadGPRHyprms(Options.GetOptLoadHyprms());
        } else {
            entgpr.SetSigmaF2(Options.GetOptSigmaF2());
            entgpr.SetWFac(Options.GetOptWFac());
            entgpr.SetSigmaN2(Options.GetOptSigmaN2());
        }

        entgpr.SetIncludeError(Options.GetOptWithError());

        entgpr.SetRCond(Options.GetOptRCond());
        entgpr.SetLAMethod(Options.GetOptLAMethod());
        entgpr.SetCalcLogPL(Options.GetOptGPRCalcLogPL());

        if( Options.IsOptMFInfoSet() ){
            entgpr.PrepForMFInfo();
        }

        if(entgpr.Interpolate(vout) == false) {
            ES_ERROR("unable to interpolate energy");
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
            ENE->ApplyOffset(Options.GetOptOffset() - ENE->GetGlobalMinimumValue());
        } else {
            ENE->ApplyOffset(Options.GetOptOffset());
        }
    }

    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            ENE->AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            ENE->AdaptUnsampledToMaxEnergy();
        }
    }

// -----------------------------------------------------------------------------
// print energy surface

    if(PrintENE() == false) {
        ES_ERROR("unable to print energy");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CPMFEnergy::AdjustGlobalMin(void)
{
    CSimpleVector<double>  GPos;
    GPos.CreateVector(ENE->GetNumOfCVs());

// adjust global minimum
    if( Options.IsOptGlobalMinSet()  ){
        string sspec(Options.GetOptGlobalMin());

        // remove "x" from the string
        replace (sspec.begin(), sspec.end(), 'x' , ' ');

        // parse values of CVs
        stringstream str(sspec);
        for(int i=0; i < ENE->GetNumOfCVs(); i++){
            double val;
            str >> val;
            if( ! str ){
                CSmallString error;
                error << "unable to decode CV value for position: " << i+1;
                RUNTIME_ERROR(error);
            }
            GPos[i] = ENE->GetCV(i)->GetIntValue(val);
        }

        // GPos.CreateVector(NCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << ENE->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < ENE->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << ENE->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << endl;

        vout << "      Closest bin found at: ";
        // find the closest bin
        CSimpleVector<double>   pos;
        pos.CreateVector(ENE->GetNumOfCVs());
        double minv = 0.0;
        int    glb_bin = 0;
        for(int ibin=0; ibin < ENE->GetNumOfBins(); ibin++){
            ENE->GetPoint(ibin,pos);
            double dist2 = 0.0;
            for(int cv=0; cv < ENE->GetNumOfCVs(); cv++){
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

        ENE->GetPoint(glb_bin,GPos);

        vout << setprecision(5) << ENE->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < ENE->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << ENE->GetCV(i)->GetRealValue(GPos[i]);
        }

        double glb_min = ENE->GetEnergy(glb_bin);
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(int ibin=0; ibin < ENE->GetNumOfBins(); ibin++){
            int samples = ENE->GetNumOfSamples(ibin);
            if( samples != 0 ){
                double ene = ENE->GetEnergy(ibin);
                ENE->SetEnergy(ibin,ene-glb_min);
            }
        }

    } else {
        // search for global minimum
        double glb_min = 0.0;
        for(int ibin=0; ibin < ENE->GetNumOfBins(); ibin++){
            int samples = ENE->GetNumOfSamples(ibin);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = ENE->GetEnergy(ibin);
            if( (ibin == 0) || (glb_min > value) ){
                glb_min = value;
                ENE->GetPoint(ibin,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum found at: ";
        vout << setprecision(5) << ENE->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < ENE->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << ENE->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(int ibin=0; ibin < ENE->GetNumOfBins(); ibin++){
            int samples = ENE->GetNumOfSamples(ibin);
            if( samples != 0 ){
                double ene = ENE->GetEnergy(ibin);
                ENE->SetEnergy(ibin,ene-glb_min);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CPMFEnergy::GetRawEnthalpy(void)
{
    for(int ibin=0; ibin < EneProxy->GetNumOfBins(); ibin++){
        int    nsamples = EneProxy->GetNumOfSamples(ibin);
        double ent = EneProxy->GetValue(ibin,E_PROXY_MEAN);
        double error = EneProxy->GetValue(ibin,E_PROXY_SEM);
        ENE->SetNumOfSamples(ibin,nsamples);
        ENE->SetEnergy(ibin,ent);
        ENE->SetError(ibin,error);
    }
}

//------------------------------------------------------------------------------

bool CPMFEnergy::PrintENE(void)
{
    vout << endl;
    vout << format("%02d:Writing results to file ...")%State << endl;
    vout << format("   ** Name: %s")%string(Options.GetArgENEFile()) << endl;

    if( OutputFile.Open(Options.GetArgENEFile(),"w") == false ){
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
    printer.SetPrintedES(ENE);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output energy file");
        return(false);
    }
    vout << "   Done." << endl;

    return(true);
}

//------------------------------------------------------------------------------

void CPMFEnergy::WriteHeader(void)
{
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        Options.PrintOptions(OutputFile);
        Accu->PrintInfo(OutputFile);
    }
}

//------------------------------------------------------------------------------

void CPMFEnergy::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = EneProxy->GetNumOfBins();
    int    sampled = 0;
    int    limit = 0;
    for(int ibin=0; ibin < EneProxy->GetNumOfBins(); ibin++) {
        if( EneProxy->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( EneProxy->GetNumOfSamples(ibin) > Options.GetOptLimit() ) {
            limit++;
        } else {
            EneProxy->SetNumOfSamples(ibin,0);
        }
    }
    if( maxbins > 0 ){
        vout << "   Sampled area: "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" ;
        vout << " ... Within limit: "
             << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%";
    }
    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFEnergy::Finalize(void)
{
    // close files if they are own by program
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-energy terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

