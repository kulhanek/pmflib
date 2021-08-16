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
#include <PMFProxy_dH.hpp>

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

    AccuName = Options.GetProgArg(0);
    HEOutputName = Options.GetProgArg(1);

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-enthalpy (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(AccuName != "-") {
        vout << "#  accumulator file (in) : " << AccuName << endl;
    } else {
        vout << "#  accumulator file (in) : - (standard input)" << endl;
    }
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
    if( InputFile.Open(Options.GetArgAccuName(),"r") == false ){
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

bool CEnthalpy::Run(void)
{
// load accumulator
    State = 1;

// -----------------------------------------------------------------------------
// setup accu, energy proxy, and output FES
    Accu        = CPMFAccumulatorPtr(new CPMFAccumulator);
    HES         = CEnergySurfacePtr(new CEnergySurface);
    EneProxy    = CPMFProxy_dH_Ptr();

    if( Options.GetOptRealm() == "<Etot>" ){
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_ETOT);
        EneProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "<Epot>" ) {
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EPOT);
        EneProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "<Ekin>" ) {
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EKIN);
        EneProxy = proxy;
// -----------------------------------------------
    } else if ( Options.GetOptRealm() == "<Erst>" ) {
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_ERST);
        EneProxy = proxy;
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << Options.GetOptRealm() ;
        RUNTIME_ERROR(error);
    }

    vout << endl;
    vout << format("%02d:Loading  accumulator: %s")%State%string(AccuName) << endl;
    State++;
    try {
        Accu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input  accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    // print Accu info
    Accu->PrintInfo(vout);
    EneProxy->Init(Accu);

    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accu->SetNCorr(Options.GetOptNCorr());

// -----------------------------------------------------------------------------
    HES->Allocate(Accu);
    HES->SetSLevel(Options.GetOptSLevel());

    vout << endl;
    vout << format("%02d:Statistics of input  accumulator")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done." << endl;

// sampling limit -------------------------------
    vout << endl;
    vout << format("%02d:Preparing  accumulator for analysis (sampling limit)")%State << endl;
    State++;
    PrepareAccumulatorI();
    PrintSampledStat();
    vout << "   Done." << endl;

// -----------------------------------------------------------------------------

    if( Options.GetOptMethod() == "raw" ){
        Accu->SetNCorr(Options.GetOptNCorr());
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
        Accu->SetNCorr(1.0);    // to prevent possible errors if used later

    } else if ( Options.GetOptMethod() == "gpr" ){
        vout << endl;
        vout << format("%02d:GPR interpolated enthalpy")%State << endl;
        CSmootherGPR   entgpr;

        entgpr.SetInputEnergyProxy(EneProxy);
        entgpr.SetOutputES(HES);

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
    string sspec(Options.GetOptGlobalMin());

    // remove "x" from the string
    replace (sspec.begin(), sspec.end(), 'x' , ' ');

    // parse values of CVs
    CSimpleVector<double>  GPos;
    GPos.CreateVector(Accu->GetNumOfCVs());
    stringstream str(sspec);
    for(int i=0; i < Accu->GetNumOfCVs(); i++){
        double val;
        str >> val;
        if( ! str ){
            CSmallString error;
            error << "unable to decode CV value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        GPos[i] = Accu->GetCV(i)->GetIntValue(val);
    }

// adjust global minimum
    if( Options.IsOptGlobalMinSet()  ){
        // GPos.CreateVector(NCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << Accu->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < Accu->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << Accu->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << endl;

        vout << "      Closest global minimum found at: ";
        // find the closest bin
        CSimpleVector<double>   pos;
        pos.CreateVector(Accu->GetNumOfCVs());
        double minv = 0.0;
        int    glb_bin = 0;
        for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
            Accu->GetPoint(ibin,pos);
            double dist2 = 0.0;
            for(int cv=0; cv < Accu->GetNumOfCVs(); cv++){
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

        Accu->GetPoint(glb_bin,GPos);

        vout << setprecision(5) << Accu->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < Accu->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << Accu->GetCV(i)->GetRealValue(GPos[i]);
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
        GPos.CreateVector(Accu->GetNumOfCVs());
        double glb_min = 0.0;
        for(int ibin=0; ibin < HES->GetNumOfBins(); ibin++){
            int samples = HES->GetNumOfSamples(ibin);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = HES->GetEnergy(ibin);
            if( (ibin == 0) || (glb_min > value) ){
                glb_min = value;
                Accu->GetPoint(ibin,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum found at: ";
        vout << setprecision(5) << Accu->GetCV(0)->GetRealValue(GPos[0]);
        for(int i=1; i < Accu->GetNumOfCVs(); i++){
            vout << "x" << setprecision(5) << Accu->GetCV(i)->GetRealValue(GPos[i]);
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

void CEnthalpy::LoadGPRHyprms(CSmootherGPR& gpr)
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

void CEnthalpy::GetRawEnthalpy(void)
{
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
        int    nsamples = EneProxy->GetNumOfSamples(ibin);
        double ent = EneProxy->GetValue(ibin,E_PROXY_VALUE);
        double error = EneProxy->GetValue(ibin,E_PROXY_ERROR);
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
        Accu->PrintInfo(OutputFile);
    }
}

//------------------------------------------------------------------------------

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero

void CEnthalpy::PrepareAccumulatorI(void)
{
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        // erase datapoints not properly sampled, preserve glueing
        if( (EneProxy->GetNumOfSamples(ibin) >= 0) && (EneProxy->GetNumOfSamples(ibin) <= Options.GetOptLimit()) ) {
            EneProxy->SetNumOfSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

void CEnthalpy::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = Accu->GetNumOfBins();
    int    sampled = 0;
    int    holes = 0;
    int    glued = 0;
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( EneProxy->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( EneProxy->GetNumOfSamples(ibin) < 0 ) {
            glued++;
        }
        if( EneProxy->GetNumOfSamples(ibin) == -1 ) {
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

void CEnthalpy::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
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

