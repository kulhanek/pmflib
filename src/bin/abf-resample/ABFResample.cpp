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

#include "ABFResample.hpp"
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

MAIN_ENTRY(CABFResample)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFResample::CABFResample(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFResample::Init(int argc,char* argv[])
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
    vout << "# abf-resample (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# Input accumulator  : " << Options.GetArgInAccuName() << endl;
    vout << "# New number of bins : " << Options.GetArgNBins() << endl;
    vout << "# Output accumulator : " << Options.GetArgOutAccuName() << endl;
    vout << "# ------------------------------------------------" << endl;

// ---------------------------------------
    vout << "# Linear algebra     : " << Options.GetOptLAMethod() << endl;
    if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
    vout << "# SVD rcond          : " << setprecision(3) << Options.GetOptRCond() << endl;
    }
// ---------------------------------------
    vout << "# ------------------------------------------------" << endl;
    if( Options.IsOptLoadHyprmsSet() ){
    vout << "# GPR hyperprms file : " << Options.GetOptLoadHyprms() << endl;
        // actual values are printed in detailed output from integrator
    } else {
    vout << "# SigmaF2            : " << setprecision(3) << Options.GetOptSigmaF2() << endl;
    vout << "# NCorr              : " << setprecision(3) << Options.GetOptNCorr() << endl;
    vout << "# Width factor wfac  : " << Options.GetOptWFac() << endl;
    vout << "# SigmaN2            : " << Options.GetOptSigmaN2() << endl;
    }

    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
    vout << "# Sampling limit     : all bins will be taken into account" << endl;
    } else {
    vout << "# Sampling limit     : " << Options.GetOptLimit() << endl;
    }

    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFResample::Run(void)
{
    State = 1;

// load accumulator
    vout << endl;
    vout << format("%02d:Loading the input PMF accumulator ...")%State << endl;
    State++;
    vout << format("   ** PMF Accumulator: %s")%string(Options.GetArgInAccuName()) << endl;
    InAccu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        InAccu->Load(Options.GetArgInAccuName());
    } catch(...) {
        CSmallString error;
        error << "unable to load the input PMF accumulator file '" << Options.GetArgInAccuName() << "'";
        ES_ERROR(error);
        return(false);
    }
    vout << "   Done" << endl;

// realms
    vout << endl;
    vout << format("%02d:Initializing the dG (ABF) realm ...")%State << endl;
    State++;

    if( CABFProxy_dG::IsCompatible(InAccu) ){
        DerProxy = CABFProxy_dG_Ptr(new CABFProxy_dG);
    } else {
        CSmallString error;
        error << "incompatible method: " << InAccu->GetMethod() << " with requested realm: dG (ABF)";
        RUNTIME_ERROR(error);
    }
    DerProxy->Init(InAccu);
    vout << "   Done" << endl;

    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accu->SetNCorr(Options.GetOptNCorr());
    FES = CEnergySurfacePtr(new CEnergySurface);
    FES->Allocate(InAccu);
    FES->SetSLevel(1.0);

    vout << endl;
    vout << format("%02d:Statistics of input PMF accumulator")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done." << endl;

// resample accu
    vout << endl;
    vout << format("%02d:Preparing the resampled PMF accumulator ...")%State << endl;
    State++;
    vout <<        "   ** Original discretization(s) = ";
    vout << InAccu->GetCV(0)->GetNumOfBins();
    for(int i=1; i < InAccu->GetNumOfCVs(); i++){
        vout << "x" << InAccu->GetCV(i)->GetNumOfBins();
    }
    vout << endl;
    vout <<        "   ** New discretization(s)      = " << Options.GetArgNBins() << endl;
    DecodeIList(Options.GetArgNBins(),NBins,"nbins");
    PrepResampledAccu();
    vout << "   Done." << endl;

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:ABF accumulator integration (GPR)")%State << endl;
    State++;
    if( Integrate() == false ) return(false);
    vout << "   Done." << endl;

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d:Resampling ABF forces ...")%State << endl;
    State++;
    ResampleAccu();
    vout << "   Done." << endl;

// print result ---------------------------------
    vout << endl;
    vout << format("%02d:Writing results to PMF accumulator: %s")%State%string(Options.GetArgOutAccuName()) << endl;
    State++;

    State++;
    try {
        OutAccu->Save(Options.GetArgOutAccuName());
    } catch(...) {
        ES_ERROR("unable to save the PMF accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    return(true);
}

//------------------------------------------------------------------------------

void CABFResample::PrepResampledAccu(void)
{
    OutAccu = InAccu->DuplicateHeader();
    OutAccu->ResampleCVs(InAccu,NBins);

//    call pmf_accu_write_rbuf_B(abfaccu%PMFAccuType,iounit,'NSAMPLES',   'AD',abfaccu%nsamples)
//    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'MICF',       'WA',abfaccu%micf,  'NSAMPLES')
//    call pmf_accu_write_rbuf_M(abfaccu%PMFAccuType,iounit,'M2ICF',      'M2',abfaccu%m2icf, 'NSAMPLES','MICF')

    NSAMPLES = OutAccu->CreateSectionData("NSAMPLES",   "AD","R","B");
    MICF     = OutAccu->CreateSectionData("MICF",       "WA","R","M","NSAMPLES");
    M2ICF    = OutAccu->CreateSectionData("M2ICF",      "M2","R","M","NSAMPLES","MICF");
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFResample::Integrate(void)
{
    Integrator.SetOutputES(FES);
    Integrator.AddInputEnergyDerProxy(DerProxy);

    if( Options.IsOptLoadHyprmsSet() ){
        LoadGPRHyprms(Integrator);
    } else {
        Integrator.SetSigmaF2(Options.GetOptSigmaF2());
        Integrator.SetNCorr(Options.GetOptNCorr());
        Integrator.SetWFac(Options.GetOptWFac());
        Integrator.SetSigmaN2(Options.GetOptSigmaN2());
    }

    Integrator.SetFastError(true);
    Integrator.SetIncludeError(false);
    Integrator.SetNoEnergy(false);
    Integrator.SetUseNumDiff(false);
    Integrator.IncludeGluedAreas(false);

    Integrator.SetRCond(Options.GetOptRCond());
    Integrator.SetLAMethod(Options.GetOptLAMethod());
    Integrator.SetUseInv(false);
    Integrator.SetKernel(Options.GetOptGPRKernel());
    Integrator.SetCalcLogPL(false);
    Integrator.SetUseZeroPoint(false);

    if(Integrator.Integrate(vout) == false) {
        ES_ERROR("unable to integrate ABF accumulator");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFResample::ResampleAccu(void)
{
    CSimpleVector<double> pos;
    pos.CreateVector(OutAccu->GetNumOfCVs());

    for(int ibin=0; ibin < OutAccu->GetNumOfBins(); ibin++){
        OutAccu->GetPoint(ibin,pos);
        for(int icv=0; icv < OutAccu->GetNumOfCVs(); icv++){
            double mf = Integrator.GetMeanForce(pos,icv);
            NSAMPLES->SetData(ibin,1.0);
            MICF->SetData(ibin,icv,mf);
            M2ICF->SetData(ibin,icv,0.0);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFResample::LoadGPRHyprms(CIntegratorGPR& gpr)
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
        } else if( key.find("SigmaN2#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode sigman2 key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            gpr.SetSigmaN2(cvind,value);
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

void CABFResample::PrintSampledStat(void)
{
    CPMFAccumulatorPtr accu = InAccu;

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
//            accu->SetNumOfSamples(ibin,1);
//            for(int icv=0; icv < accu->GetNumOfCVs(); icv++){
//                accu->SetData("MICF",ibin,icv,0.0);
//            }
        }
    }
    if( maxbins > 0 ){
        vout << "   ** Sampled area: "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" ;
        vout << " ... Within limit: "
             << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%";
    }
    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFResample::DecodeIList(const CSmallString& spec, std::vector<int>& ilist,const CSmallString& optionname)
{
    int ncvs = FES->GetNumOfCVs();

    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( (int)slist.size() > ncvs ){
        CSmallString error;
        error << "too many values (" << slist.size() << ") for " << optionname << " than required (" << ncvs << ")";
        RUNTIME_ERROR(error);
    }

    ilist.resize(ncvs);

    // parse values
    int  last_value = 0;
    for(int i=0; i < (int)slist.size(); i++){
        stringstream str(slist[i]);
        int  value = 0;
        str >> value;
        if( ! str ){
            CSmallString error;
            error << "unable to decode value for " << optionname << " at position: " << i+1;
            RUNTIME_ERROR(error);
        }
        ilist[i] = value;
        last_value = value;
    }

    // pad the rest with the last value
    for(int i=slist.size(); i < ncvs; i++){
        ilist[i] = last_value;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFResample::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    CSmallTime dur;
    dur = dt - StartTime;

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-resample terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

