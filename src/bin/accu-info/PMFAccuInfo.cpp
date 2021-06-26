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
    } else if( (Options.GetProgArg(1) == "NSAMPLES") || (Options.GetProgArg(1) == "nsamples") ){
        GetSection("NSAMPLES");
// -----------------------------------------------
    } else if( Options.GetProgArg(1) == "micf" ){
        GetMICF();
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

void CPMFAccuInfo::GetMICF(void)
{
// prepare accumulator --------------------------
    vout << "#" << endl;
    vout << "# Data section: MICF (from ABF)" << endl;

    vout << high;

    CABFProxy_dG der_proxy;
    der_proxy.Init(Accu);

    Values.CreateVector(Accu->GetNumOfBins());
    Sigmas.CreateVector(Accu->GetNumOfBins());
    Errors.CreateVector(Accu->GetNumOfBins());

    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++){
        Values[ibin] = der_proxy.GetValue(ibin,Options.GetOptCV()-1,E_PROXY_VALUE);
        Sigmas[ibin] = der_proxy.GetValue(ibin,Options.GetOptCV()-1,E_PROXY_SIGMA);
        Errors[ibin] = der_proxy.GetValue(ibin,Options.GetOptCV()-1,E_PROXY_ERROR);
    }

    if( Options.GetOptNoHeader() == false ){
        PrintHeader("MICF (from ABF)",true);
    }

    // print data
    PrintData(true);

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

