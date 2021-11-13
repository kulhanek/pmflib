// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <SmallTimeAndDate.hpp>
#include "PMFAccuInfo.hpp"
#include <iomanip>
#include <boost/format.hpp>
#include <ABFProxy_dG.hpp>
#include <ABFProxy_mTdS.hpp>
#include <CSTProxy_dG.hpp>
#include <CSTProxy_mTdS.hpp>
#include <CSTProxy_MTC.hpp>
#include <PMFProxy_dH.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;


//------------------------------------------------------------------------------

MAIN_ENTRY(CPMFAccuInfo)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccuInfo::CPMFAccuInfo(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFAccuInfo::Init(int argc,char* argv[])
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

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << debug;
    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-info (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetProgArg(0) != "-") {
        vout << "# PMF accumulator (in)  : " << Options.GetProgArg(0) << endl;
    } else {
        vout << "# PMF accumulator (in)  : - (standard input)" << endl;
    }
    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetProgArg(0),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CPMFAccuInfo::Run(void)
{
    Accu = CPMFAccumulatorPtr(new CPMFAccumulator);

// load accumulator
    vout << "#" << endl;
    vout << "# 1) Loading the PMF accumulator: " << Options.GetProgArg(0) << endl;
    try {
        Accu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input PMF accumulator file");
        return(false);
    }
    vout << "#   Done" << endl;

// perform action
    if( Options.GetNumberOfProgArgs() == 1 ) {
        PrintInfo();
// -----------------------------------------------
    } else  if( Options.GetProgArg(1)  == "info" ){
        PrintInfo();
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "list-sections" ){
        ListSections();
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-section" ){
        GetSection(Options.GetProgArg(2));
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-derivative" ){
        GetDerivative(Options.GetProgArg(2));
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-energy" ){
        GetEnergy(Options.GetProgArg(2));
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-mtc" ){
        GetMTC();
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-mean" ){
        GetMean(Options.GetProgArg(2));
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "get-tseries" ){
        GetTSeries(Options.GetProgArg(2));
// -----------------------------------------------
    } else if( (Options.GetProgArg(1) == "NSAMPLES") || (Options.GetProgArg(1) == "nsamples") ){
        GetSection("NSAMPLES");
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unrecognized action: " << Options.GetProgArg(1);
        ES_ERROR(error);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuInfo::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << "#" << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-info terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        vout << high;
        ErrorSystem.PrintErrors(vout);
        vout << endl;
        vout << debug;
    } else {
        vout << "#" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuInfo::PrintInfo(void)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# 2) Info about the PMF accumulator."<< endl;

    vout << high;
    // print CVS info
    Accu->PrintAccuInfo(vout);
    Accu->PrintCVSInfo(vout);
    vout << endl;

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::ListSections(void)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# 2) List of data sections stored in the PMF accumulator."<< endl;

    vout << high;
    Accu->ListSections(vout);
    vout << endl;

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetSection(const CSmallString& name)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# Data section: " << name << endl;

    vout << high;

    CPMFAccuDataPtr data = Accu->GetSectionData(name);

// -----------------------------------------------
    if( (data->GetMode() == "M") || (data->GetMode() == "B") ){
        if( Options.GetOptNoHeader() == false ){
            PrintHeader(name,false);
        }
        // get data
        Values.CreateVector(Accu->GetNumOfBins());

        if( data->GetMode() == "M" ){
            for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
                Values[ibin] = data->GetData(ibin,Options.GetOptCV()-1);
            }
        }

        if( data->GetMode() == "B" ){
            for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
                Values[ibin] = data->GetData(ibin);
            }
        }

        // print data
        PrintData(false);
// -----------------------------------------------
    }  else  if (data->GetMode() == "C") {

        // get data
        Values.CreateVector(Accu->GetNumOfCVs());
        for(int icv=0; icv < Accu->GetNumOfCVs(); icv++){
            Values[icv] = data->GetData(icv);
        }

        PrintDataPerCV(name);
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported data section mode (" << data->GetMode() << ")";
        RUNTIME_ERROR(error);
    }

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetDerivative(const CSmallString& name)
{
    CEnergyDerProxyPtr      der_proxy;

// init energy-der proxy
    if( name == "dG" ){
        if( CABFProxy_dG::IsCompatible(Accu) ){
            der_proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        } else if (CCSTProxy_dG::IsCompatible(Accu) ) {
            der_proxy    = CCSTProxy_dG_Ptr(new CCSTProxy_dG);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm: " <<  name;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( name == "dG_p" ) {
        CABFProxy_dG_Ptr proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        proxy->SetType(ABF_MICF_POT);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "dG_k" ) {
        CABFProxy_dG_Ptr proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        proxy->SetType(ABF_MICF_KIN);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS" ) {
        if( CABFProxy_mTdS::IsCompatible(Accu) ){
            der_proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        } else if (CCSTProxy_mTdS::IsCompatible(Accu) ) {
            der_proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm: " <<  name;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( name == "-TdS_PP" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_PP);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_PK" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_PK);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_PR" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_PR);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_KP" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_KP);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_KK" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_KK);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_KR" ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_KR);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_HP" ) {
        CCSTProxy_mTdS_Ptr proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        proxy->SetType(CST_C11HP);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_HK" ) {
        CCSTProxy_mTdS_Ptr proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        proxy->SetType(CST_C11HK);
        der_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "-TdS_HR" ) {
        CCSTProxy_mTdS_Ptr proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        proxy->SetType(CST_C11HR);
        der_proxy = proxy;
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << name;
        RUNTIME_ERROR(error);
    }

    der_proxy->Init(Accu);

    vout << high;

    Values.CreateVector(Accu->GetNumOfBins());
    Sigmas.CreateVector(Accu->GetNumOfBins());
    Errors.CreateVector(Accu->GetNumOfBins());

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
        Values[ibin] = der_proxy->GetValue(ibin,Options.GetOptCV()-1,E_PROXY_VALUE);
        Sigmas[ibin] = der_proxy->GetValue(ibin,Options.GetOptCV()-1,E_PROXY_SIGMA);
        Errors[ibin] = der_proxy->GetValue(ibin,Options.GetOptCV()-1,E_PROXY_ERROR);
    }

    if( Options.GetOptNoHeader() == false ){
        CSmallString title;
        title << "derivative for: " << name;
        PrintHeader(title,true);
    }

    // print data
    PrintData(true);

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetEnergy(const CSmallString& name)
{
    CPMFProxy_dH_Ptr    ene_proxy;

    if( name == "<Etot>" ){
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_ETOT);
        ene_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "<Epot>" ) {
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EPOT);
        ene_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "<Ekin>" ) {
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_EKIN);
        ene_proxy = proxy;
// -----------------------------------------------
    } else if ( name == "<Erst>" ) {
        CPMFProxy_dH_Ptr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy->SetType(PMF_ERST);
        ene_proxy = proxy;
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << name ;
        RUNTIME_ERROR(error);
    }

    ene_proxy->Init(Accu);

    vout << high;

    Values.CreateVector(Accu->GetNumOfBins());
    Sigmas.CreateVector(Accu->GetNumOfBins());
    Errors.CreateVector(Accu->GetNumOfBins());

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
        Values[ibin] = ene_proxy->GetValue(ibin,E_PROXY_VALUE);
        Sigmas[ibin] = ene_proxy->GetValue(ibin,E_PROXY_SIGMA);
        Errors[ibin] = ene_proxy->GetValue(ibin,E_PROXY_ERROR);
    }

    if( Options.GetOptNoHeader() == false ){
        CSmallString title;
        title << "derivative for: " << name;
        PrintHeader(title,true);
    }

    // print data
    PrintData(true);

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetMTC(void)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# Data section: MTC (from CST)" << endl;

    vout << high;

    CCSTProxy_MTC_Ptr mtc_proxy   = CCSTProxy_MTC_Ptr(new CCSTProxy_MTC);
    mtc_proxy->Init(Accu);

    Values.CreateVector(Accu->GetNumOfBins());
    Sigmas.CreateVector(Accu->GetNumOfBins());
    Errors.CreateVector(Accu->GetNumOfBins());

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
        Values[ibin] = mtc_proxy->GetValue(ibin,E_PROXY_VALUE);
        Sigmas[ibin] = mtc_proxy->GetValue(ibin,E_PROXY_SIGMA);
        Errors[ibin] = mtc_proxy->GetValue(ibin,E_PROXY_ERROR);
    }

    if( Options.GetOptNoHeader() == false ){
        PrintHeader("MTC (from CST)",true);
    }

    // print data
    PrintData(true);

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetMean(const CSmallString& name)
{
    CSmallString mean;
    CSmallString m2;

    mean << "M" << name;
    m2   << "M2" << name;

    vout << "#" << endl;
    vout << "# mean for : " << mean << " (" << m2 << ")" << endl;

    vout << high;

    CPMFAccuDataPtr mval_sd  = Accu->GetSectionData(mean);
    CPMFAccuDataPtr m2val_sd = Accu->GetSectionData(m2);

    if( mval_sd->GetMSName() != m2val_sd->GetMSName() ){
        RUNTIME_ERROR("inconsistent MSName");
    }

    CPMFAccuDataPtr nsamples_sd = Accu->GetSectionData(mval_sd->GetMSName());

    Values.CreateVector(Accu->GetNumOfBins());
    Sigmas.CreateVector(Accu->GetNumOfBins());
    Errors.CreateVector(Accu->GetNumOfBins());

    double  ncorr       = Accu->GetNCorr();

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){

        double  nsamples    = nsamples_sd->GetData(ibin);
        double  mene        = mval_sd->GetData(ibin,Options.GetOptCV()-1);
        double  m2ene       = m2val_sd->GetData(ibin,Options.GetOptCV()-1);

        if( nsamples > 0 ) {
            Values[ibin] = mene;
            Sigmas[ibin] = sqrt(m2ene / nsamples);
            Errors[ibin] = sqrt(m2ene * ncorr) / nsamples;
        } else {
            Values[ibin] = 0.0;
            Sigmas[ibin] = 0.0;
            Errors[ibin] = 0.0;
        }
    }

    if( Options.GetOptNoHeader() == false ){
        CSmallString title;
        title << "mean for: " << name;
        PrintHeader(title,true);
    }

    // print data
    PrintData(true);

    vout << debug;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::GetTSeries(const CSmallString& name)
{
    vout << "#" << endl;
    vout << "# time series for : " << name << endl;

    vout << high;

    CPMFAccuDataPtr tser  = Accu->GetSectionData(name);
    double          dt    = Accu->GetTimeStep();

    double time = 0.0;
    for(int tstep=0; tstep < Accu->GetNSTLimit(); tstep++){
        double  value   = tser->GetData(tstep,Options.GetOptCV()-1);
        vout << format("%15d %15.2f ")%tstep%time << format(Options.GetOptOSFormat())%value << endl;;
        time += dt;
    }

    vout << debug;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuInfo::PrintHeader(const CSmallString& sec_name, bool print_errors)
{
    // print header
    vout << "#" << endl;
    int num = Accu->GetNumOfCVs() + 1;
    if( print_errors ) num += 2;

// -----------------------------------------------
    vout << "#";
    for(int i=0; i < num; i++) {
        vout << format("%15d ")%(i+1);
    }
    vout << endl;
// -----------------------------------------------
    vout << "#";
    for(int i=0; i < num ; i++) {
        vout << "--------------- ";
    }
    vout << endl;
// -----------------------------------------------
    vout << "#";
    for(int i=0; i < Accu->GetNumOfCVs(); i++) {
        vout << format("%15s ")%Accu->GetCV(i)->GetName();
    }

    if( print_errors ) {
        vout << format("%47s")%sec_name;
    } else {
        vout << format("%15s")%sec_name;
    }
    vout << endl;
// -----------------------------------------------
    vout << "#";
    for(int i=0; i < Accu->GetNumOfCVs(); i++) {
        vout << "--------------- ";
    }
    if( print_errors ) {
        vout << "-----------------------------------------------";
    } else {
        vout << "---------------";
    }
    vout << endl;
// -----------------------------------------------
    vout << "#";
    for(int i=0; i < Accu->GetNumOfCVs(); i++) {
        CSmallString unit;
        unit << "[" << Accu->GetCV(i)->GetUnit() << "]";
        vout << format("%15s ")%unit;
    }

    if( print_errors ) {
        vout << "     <<X> [i.u.]     s(X) [i.u.]   s(<<X>) [i.u.]" << endl;
    } else {
        vout << "         [i.u.]" << endl;
    }

// -----------------------------------------------
    vout << "#";
    for(int i=0; i < num; i++) {
        vout << "--------------- ";
    }
    vout << endl;
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::PrintData(bool print_errors)
{
    CSimpleVector<double>   pos;
    CSimpleVector<int>      ipos;

    pos.CreateVector(Accu->GetNumOfCVs());
    ipos.CreateVector(Accu->GetNumOfCVs());

    int last_cv = -1;

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){

        // do we have enough samples?
        double nsamples = Accu->GetNumOfSamples(ibin);
        if( Options.IsOptLimitSet() ){
            if( nsamples < Options.GetOptLimit() ) continue;
        }

    // write block delimiter - required by GNUPlot
        if( Options.GetOptNoGNUPlot() == false ) {
            int ncvs = Accu->GetNumOfCVs();
            Accu->GetIPoint(ibin,ipos);

            if( (last_cv >= 0) && (ipos[ncvs-1] != last_cv + 1) ){
                vout << endl;
            }
            last_cv = ipos[ncvs-1];
        }

        Accu->GetPoint(ibin,pos);

        CSmallString xform;
        xform = " " + Options.GetOptIXFormat();

        // print point position
        for(int i=0; i < Accu->GetNumOfCVs(); i++) {
            double xvalue = Accu->GetCV(i)->GetRealValue(pos[i]);
            vout << format(xform)%xvalue;
        }

        CSmallString yform;
        yform = " " + Options.GetOptOSFormat();
        // and value

        double value = Values[ibin];
        vout << format(yform)%value;

        if( print_errors ){
            double sigma = Sigmas[ibin];
            vout << format(yform)%sigma;

            double error = Errors[ibin];
            vout << format(yform)%error;
        }

        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CPMFAccuInfo::PrintDataPerCV(const CSmallString& sec_name)
{
// print header
    vout << "#" << endl;
    vout << "#              1               2" << endl;
    vout << "#--------------- ---------------" << endl;
    vout << "#             CV " << format("%15s")%sec_name << endl;
    vout << "#--------------- ---------------" << endl;
    vout << "#                          [i.u.]" << endl;
    vout << "#--------------- ---------------" << endl;

// print data
    CSmallString form;
    form = " %15d " + Options.GetOptOSFormat();
    for(int icv=0; icv < Accu->GetNumOfCVs(); icv++) {
        double value = Values[icv];
        vout << format(form)%(icv+1)%value << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

