// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <GHSIntegratorGPR0A.hpp>
#include <ErrorSystem.hpp>
#include <FortranMatrix.hpp>
#include <Vector.hpp>
#include <algorithm>
#include <SciLapack.hpp>
#include <iomanip>
#include <SciBlas.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <set>

// OpenMP support
#if defined(_OPENMP)
#include <omp.h>
#endif

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGHSIntegratorGPR0A::CGHSIntegratorGPR0A(void)
{
    GDerProxy           = NULL;
    HDerProxy           = NULL;
    SDerProxy           = NULL;
    GSurface            = NULL;
    HSurface            = NULL;
    SSurface            = NULL;

    GPRSize             = 0;
    NumOfUsedBins       = 0;
    NumOfCVs            = 0;
    NumOfBins           = 0;

    ConstrainedTK       = false;
    IncludeError        = false;
    NoEnergy            = false;
    GlobalMinSet        = false;
    GPosSet             = false;
    GPosBin             = 0;

    UseNumDiff          = false;
    Method              = EGPRLA_LU;
    Kernel              = EGPRK_ARDSE;

    NumOfThreads        = 1;

    UseInv              = false;
    NeedInv             = false;
    FastErrors          = true;

    KSInverted          = false;
}

//------------------------------------------------------------------------------

CGHSIntegratorGPR0A::~CGHSIntegratorGPR0A(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPR0A::SetGDerProxy(CEnergyDerProxyPtr p_proxy)
{
    if( p_proxy == NULL ) return;                 // no-proxy
    if( p_proxy->GetAccu() == NULL ) return;      // no PMFAccu

    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)p_proxy->GetAccu()->GetNumOfCVs();
        NumOfBins = (size_t)p_proxy->GetAccu()->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)p_proxy->GetAccu()->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_proxy->GetAccu()->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    GDerProxy = p_proxy;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetHDerProxy(CEnergyDerProxyPtr p_proxy)
{
    if( p_proxy == NULL ) return;                 // no-proxy
    if( p_proxy->GetAccu() == NULL ) return;      // no PMFAccu

    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)p_proxy->GetAccu()->GetNumOfCVs();
        NumOfBins = (size_t)p_proxy->GetAccu()->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)p_proxy->GetAccu()->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_proxy->GetAccu()->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    HDerProxy = p_proxy;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSDerProxy(CEnergyDerProxyPtr p_proxy)
{
    if( p_proxy == NULL ) return;                 // no-proxy
    if( p_proxy->GetAccu() == NULL ) return;      // no PMFAccu

    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)p_proxy->GetAccu()->GetNumOfCVs();
        NumOfBins = (size_t)p_proxy->GetAccu()->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)p_proxy->GetAccu()->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_proxy->GetAccu()->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    SDerProxy = p_proxy;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetOutputFES(CEnergySurfacePtr p_surf)
{
    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)p_surf->GetNumOfCVs();
        NumOfBins = (size_t)p_surf->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)p_surf->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_surf->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    GSurface = p_surf;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetOutputHES(CEnergySurfacePtr p_surf)
{
    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)p_surf->GetNumOfCVs();
        NumOfBins = (size_t)p_surf->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)p_surf->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_surf->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    HSurface = p_surf;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetOutputSES(CEnergySurfacePtr p_surf)
{
    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)p_surf->GetNumOfCVs();
        NumOfBins = (size_t)p_surf->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)p_surf->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_surf->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    SSurface = p_surf;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSigmaF2(const CSmallString& spec)
{
    string          sspec(spec);
    vector<string>  ssigmaf2;

    split(ssigmaf2,sspec,is_any_of("x"),token_compress_on);

    if( ssigmaf2.size() > 3){
        CSmallString error;
        error << "too many sigmaf2 (" << ssigmaf2.size() << ") than required (" << 3 << ")";
        RUNTIME_ERROR(error);
    }

    SigmaF2.CreateVector(3);

    // parse values of sigmaf2
    double last_sigmaf2 = 0.1;
    for(size_t i=0; i < ssigmaf2.size(); i++){
        stringstream str(ssigmaf2[i]);
        str >> last_sigmaf2;
        if( ! str ){
            CSmallString error;
            error << "unable to decode sigmaf2 value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        SigmaF2[i] = last_sigmaf2;
    }

    // pad the rest with the last value
    for(size_t i=ssigmaf2.size(); i < 3; i++){
        SigmaF2[i] = last_sigmaf2;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSigmaF2(CSimpleVector<double>& sigmaf2)
{
    if( sigmaf2.GetLength() != 3 ){
        RUNTIME_ERROR("dimmension inconsistent in the source and target");
    }
    SigmaF2 = sigmaf2;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSigmaF2(size_t cvind, double value)
{
    if( cvind >= 3 ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is SigmaF2 initialized?
    if( SigmaF2.GetLength() == 0 ){
        SigmaF2.CreateVector(3);
    }
    SigmaF2[cvind] = value;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetCoVar(const CSmallString& spec)
{
    string          sspec(spec);
    vector<string>  scovar;

    split(scovar,sspec,is_any_of("x"),token_compress_on);

    if( scovar.size() > 3){
        CSmallString error;
        error << "too many covariances (" << scovar.size() << ") than required (" << 3 << ")";
        RUNTIME_ERROR(error);
    }

    CoVar.CreateVector(3);

    // parse values of scovar
    double last_covar = 0.1;
    for(size_t i=0; i < scovar.size(); i++){
        stringstream str(scovar[i]);
        str >> last_covar;
        if( ! str ){
            CSmallString error;
            error << "unable to decode scovar value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        CoVar[i] = last_covar;
    }

    // pad the rest with the last value
    for(size_t i=scovar.size(); i < 3; i++){
        CoVar[i] = last_covar;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetCoVar(CSimpleVector<double>& covar)
{
    if( covar.GetLength() != 3 ){
        RUNTIME_ERROR("dimmension inconsistent in the source and target");
    }
    CoVar = covar;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetCoVar(size_t cvind, double value)
{
    if( cvind >= 3 ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is CoVar initialized?
    if( CoVar.GetLength() == 0 ){
        CoVar.CreateVector(3);
    }
    CoVar[cvind] = value;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSigmaN2(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetSigmaN2");
    }

    string          sspec(spec);
    vector<string>  ssigman2;

    split(ssigman2,sspec,is_any_of("x"),token_compress_on);

    if( ssigman2.size() > 3*NumOfCVs ){
        CSmallString error;
        error << "too many sigman2 (" << ssigman2.size() << ") than required (" << NumOfCVs << ")";
        RUNTIME_ERROR(error);
    }

    SigmaN2.CreateVector(3*NumOfCVs);

    // parse values of sigman2
    double last_sigman2 = 0.1;
    for(size_t i=0; i < ssigman2.size(); i++){
        stringstream str(ssigman2[i]);
        str >> last_sigman2;
        if( ! str ){
            CSmallString error;
            error << "unable to decode sigman2 value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        SigmaN2[i] = last_sigman2;
    }

    // pad the rest with the last value
    for(size_t i=ssigman2.size(); i < 3*NumOfCVs; i++){
        SigmaN2[i] = last_sigman2;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSigmaN2(CSimpleVector<double>& sigman2)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("FES is not set for SetSigmaN2");
    }
    if( sigman2.GetLength() != 3*NumOfCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    SigmaN2 = sigman2;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetSigmaN2(size_t cvind, double value)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("FES is not set for SetSigmaN2");
    }
    if( cvind >= 3*NumOfCVs ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is SigmaN2 initialized?
    if( SigmaN2.GetLength() == 0 ){
        SigmaN2.CreateVector(3*NumOfCVs);
    }
    SigmaN2[cvind] = value;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetWFac(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }

    string          sspec(spec);
    vector<string>  swfacs;

    split(swfacs,sspec,is_any_of("x"),token_compress_on);

    if( swfacs.size() > NumOfCVs ){
        CSmallString error;
        error << "too many wfacs (" << swfacs.size() << ") than required (" << NumOfCVs << ")";
        RUNTIME_ERROR(error);
    }

    WFac.CreateVector(NumOfCVs);

    // parse values of wfac
    double last_wfac = 3.0;
    for(size_t i=0; i < swfacs.size(); i++){
        stringstream str(swfacs[i]);
        str >> last_wfac;
        if( ! str ){
            CSmallString error;
            error << "unable to decode wfac value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        WFac[i] = last_wfac;
    }

    // pad the rest with the last value
    for(size_t i=swfacs.size(); i < NumOfCVs; i++){
        WFac[i] = last_wfac;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetWFac(CSimpleVector<double>& wfac)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("FES is not set for SetWFac");
    }
    if( wfac.GetLength() != NumOfCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    WFac = wfac;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetWFac(size_t cvind, double value)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("FES is not set for SetWFac");
    }
    if( cvind >= NumOfCVs ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is wfac initialized?
    if( WFac.GetLength() == 0 ){
        WFac.CreateVector(NumOfCVs);
    }
    WFac[cvind] = value;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::LoadGPRHyprms(const CSmallString& name)
{
    ifstream fin;
    fin.open(name);
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << name;
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
        if( key.find("SigmaF2#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode sigmaf2 key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            SetSigmaF2(cvind,value);
        } else         if( key.find("CoVar#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode covar key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            SetCoVar(cvind,value);
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
            SetWFac(cvind,value);
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
            SetSigmaN2(cvind,value);
        } else {
            CSmallString error;
            error << "GPR hyperparameters file, unrecognized key: " << key.c_str();
            RUNTIME_ERROR(error);
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::EnableConstraints(bool set)
{
    ConstrainedTK = set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetIncludeError(bool set)
{
    IncludeError = set;
    if( FastErrors == false ){
        NeedInv |= set;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetUseNumDiff(bool set)
{
    UseNumDiff = set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetNoEnergy(bool set)
{
    NoEnergy = set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetLAMethod(EGPRLAMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetLAMethod(const CSmallString& method)
{
    if( method == "svd" ){
        SetLAMethod(EGPRLA_SVD);
    } else if( method == "svd2" ){
        SetLAMethod(EGPRLA_SVD2);
    } else if( method == "lu" ) {
        SetLAMethod(EGPRLA_LU);
    } else if( method == "ll" ) {
        SetLAMethod(EGPRLA_LL);
    } else if( method == "default" ) {
        SetLAMethod(EGPRLA_LU);
    } else {
        CSmallString error;
        error << "Specified method '" << method << "' for linear algebra is not supported. "
                 "Supported methods are: svd (simple driver), svd2 (conquer and divide driver), lu, ll (Choleskff decomposition), default (=lu)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetKernel(const CSmallString& kernel)
{
    if( kernel == "ardse" ){
        Kernel = EGPRK_ARDSE;
    } else if( kernel == "ardmc52" ) {
        Kernel = EGPRK_ARDMC52;
    } else if( kernel == "default" ) {
        Kernel = EGPRK_ARDSE;
    } else {
        CSmallString error;
        error << "Specified kernel '" << kernel << "' is not supported. "
                 "Supported kernels are: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), default(=ardse)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetGlobalMin(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetGlobalMin");
    }
    if( GSurface == NULL ) {
        RUNTIME_ERROR("FES is not set");
    }

    GlobalMinSet = true;
    string sspec(spec);

    // remove "x" from the string
    replace (sspec.begin(), sspec.end(), 'x' , ' ');

    // parse values of CVs
    GPos.CreateVector(NumOfCVs);
    stringstream str(sspec);
    for(size_t i=0; i < NumOfCVs; i++){
        double val;
        str >> val;
        if( ! str ){
            CSmallString error;
            error << "unable to decode CV value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        GPos[i] = GSurface->GetCV(i)->GetIntValue(val);
    }

    GPosSet = true;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetGlobalMin(const CSimpleVector<double>& pos)
{
    GlobalMinSet = true;
    GPos = pos;
    GPosSet = true;
}

//------------------------------------------------------------------------------

CSimpleVector<double> CGHSIntegratorGPR0A::GetGlobalMin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPos);
}

//------------------------------------------------------------------------------

int CGHSIntegratorGPR0A::GetGlobalMinBin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPosBin);
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetUseInv(bool set)
{
    UseInv = set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::PrepForHyprmsGrd(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetCalcLogPL(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::SetFastError(bool set)
{
    FastErrors = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CGHSIntegratorGPR0A::Integrate(CVerboseStr& vout,bool nostat)
{
    PrintExecInfo(vout);

    if( GDerProxy == NULL ){
        RUNTIME_ERROR("no input data - FES");
    }
    if( HDerProxy == NULL ){
        RUNTIME_ERROR("no input data - HES");
    }
    if( SDerProxy == NULL ){
        RUNTIME_ERROR("no input data - SES");
    }
    if( GSurface == NULL ) {
        RUNTIME_ERROR("GES output is not set");
    }
    if( HSurface == NULL ) {
        RUNTIME_ERROR("HES output is not set");
    }
    if( SSurface == NULL ) {
        RUNTIME_ERROR("SES output is not set");
    }
    if( NumOfCVs == 0 ) {
        RUNTIME_ERROR("number of CVs is zero");
    }
    if( NumOfBins == 0 ) {
        RUNTIME_ERROR("number of bins is zero");
    }
    if( SigmaF2.GetLength() != 3 ){
        RUNTIME_ERROR("sigmaf2 is not set");
    }
    if( CoVar.GetLength() != 3 ){
        RUNTIME_ERROR("covar is not set");
    }
    if( WFac.GetLength() != NumOfCVs ){
        RUNTIME_ERROR("wfac is not set");
    }
    if( SigmaN2.GetLength() != 3*NumOfCVs ){
        RUNTIME_ERROR("sigman2 is not set");
    }

    // GPR setup
    CVLengths2.CreateVector(NumOfCVs);
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*GSurface->GetCV(i)->GetRange()/GSurface->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }

    // number of data points, use dH source as it can have smaller number of samples (NTDS vs NSAMPLES)
    NumOfUsedBins = 0;
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( HDerProxy->GetNumOfSamples(ibin) > 0 ) NumOfUsedBins++;
    }
    GPRSize = 3 * NumOfUsedBins * NumOfCVs;

    // create sampled map
    SampledMap.resize(NumOfUsedBins);
    size_t ind = 0;
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( HDerProxy->GetNumOfSamples(ibin) <= 0 ) continue;
        SampledMap[ind] = ibin;
        ind++;
    }

    // init GPR arrays
    GPRModel.CreateVector(GPRSize);
    Y.CreateVector(GPRSize);
    KS.CreateMatrix(GPRSize,GPRSize);
    TK.CreateMatrix(3,3);

    // print hyperparameters
        vout        << "   Hyperparameters ..." << endl;
    for(size_t k=0; k < 3*NumOfCVs; k++ ){
        vout << format("      SigmaF2#%-2d= %10.4f")%(k+1)%SigmaF2[k] << endl;
    }
    for(size_t k=0; k < 3*NumOfCVs; k++ ){
        vout << format("      CoVar#%-2d  = %10.4f")%(k+1)%CoVar[k] << endl;
    }

    for(size_t k=0; k < NumOfCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(size_t k=0; k < 3*NumOfCVs; k++ ){
        vout << format("      SigmaN2#%-2d= %10.4e")%(k+1)%SigmaN2[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    if( ! nostat ){
        // and finally some statistics
//        for(size_t k=0; k < NumOfCVs; k++ ){
//            vout << "      RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
//        }

        // and log of marginal likelihood
        vout << "      logML     = " << setprecision(5) << GetLogML() << endl;
        if( NeedInv || UseInv ){
            // and log of pseudo-likelihood
            vout << "      logPL     = " << setprecision(5) << GetLogPL() << endl;
        }
    }

    // finalize EneSurface if requested
    if( ! NoEnergy ){
        CalculateEnergy(vout);
        if( IncludeError ){
            if( FastErrors ){
   //             CalculateCovs(vout);
   //             CalculateErrorsFromCov(vout);
            } else {
   //             CalculateErrors(GPos,vout);
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::PrintExecInfo(CVerboseStr& vout)
{
    NumOfThreads = 1;

#if defined(_OPENMP)
    {
        NumOfThreads = omp_get_max_threads();
        vout << "   OpenMP - number of threads: " << NumOfThreads << endl;
    }
#else
    vout << "   No OpenMP - sequential mode." << endl;
#endif
    RunBlasLapackPar();
    CSciLapack::PrintExecInfo(vout);
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CGHSIntegratorGPR0A::TrainGP(CVerboseStr& vout)
{
    if( UseNumDiff ) {
        vout << "   Creating K+Sigma and Y (numeric differentation) ..." << endl;
    } else {
        vout << "   Creating K+Sigma and Y ..." << endl;
    }
        vout << "      Kernel    = " << GetKernelName() << endl;
        vout << "      Dim       = " << GPRSize << " x " << GPRSize << endl;

        if( ConstrainedTK ){
        vout << "      Task Covs = constrained" << endl;
        } else {
        vout << "      Task Covs = unconstrained" << endl;
        }

// construct Y
    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        for(size_t ii=0; ii < NumOfCVs; ii++){
            double mf0 = GDerProxy->GetValue(ibin,ii,E_PROXY_VALUE);
            Y[offset0 + indi*NumOfCVs+ii] = mf0;

            double mf1 = HDerProxy->GetValue(ibin,ii,E_PROXY_VALUE);
            Y[offset1 + indi*NumOfCVs+ii] = mf1;

            double mf2 = SDerProxy->GetValue(ibin,ii,E_PROXY_VALUE);
            Y[offset2 + indi*NumOfCVs+ii] = mf2;
        }
    }

// construct KS
    CreateTK();
    CreateKS();

    RunBlasLapackPar();

// inverting the K+Sigma
    int result = 0;
    switch(Method){
        case(EGPRLA_LU):
            if( UseInv ){
                vout << "   Inverting K+Sigma by LU ..." << endl;
                result = CSciLapack::invLU(KS,logdetK);
                if( result != 0 ) return(false);
                KSInverted = true;
                // calculate weights
                vout << "   Calculating weights B ..." << endl;
                CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
            } else {
                GPRModel = Y;
                if( NeedInv ){
                    vout << "   Training GPR + K+Sigma inversion by LU ..." << endl;
                    result = CSciLapack::solvleLUInv(KS,GPRModel,logdetK);
                    if( result != 0 ) return(false);
                    KSInverted = true;
                } else {
                    vout << "   Training GPR by LU ..." << endl;
                    result = CSciLapack::solvleLU(KS,GPRModel,logdetK);
                    if( result != 0 ) return(false);
                }
            }
            break;
        case(EGPRLA_LL):
            if( UseInv ){
                vout << "   Inverting K+Sigma by LL ..." << endl;
                result = CSciLapack::invLL(KS,logdetK);
                if( result != 0 ) return(false);
                KSInverted = true;
                // calculate weights
                vout << "   Calculating weights B ..." << endl;
                CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
            } else {
                GPRModel = Y;
                if( NeedInv ){
                    vout << "   Training GPR + K+Sigma inversion by LL ..." << endl;
                    result = CSciLapack::solvleLLInv(KS,GPRModel,logdetK);
                    if( result != 0 ) return(false);
                    KSInverted = true;
                } else {
                    vout << "   Training GPR by LL ..." << endl;
                    result = CSciLapack::solvleLL(KS,GPRModel,logdetK);
                    if( result != 0 ) return(false);
                }
            }
            break;
        case(EGPRLA_SVD):{
            vout << "   Inverting K+Sigma by SVD (divide and conquer driver) ..." << endl;
            int rank = 0;
            double realRCond = 0;
            result = CSciLapack::invSVD2(KS,logdetK,RCond,rank,realRCond);
            vout << "      Rank = " << rank << "; Info = " << result << "; Real rcond = " << scientific << realRCond << fixed << endl;
            if( result != 0 ) return(false);
            // calculate weights
            vout << "   Calculating weights B ..." << endl;
            CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
            KSInverted = true;
            }
            break;
        case(EGPRLA_SVD2):{
            vout << "   Inverting K+Sigma by SVD2 (simple driver) ..." << endl;
            int rank = 0;
            double realRCond = 0;
            result = CSciLapack::invSVD1(KS,logdetK,RCond,rank,realRCond);
            vout << "      Rank = " << rank << "; Info = " << result << "; Real rcond = " << scientific << realRCond << fixed << endl;
            if( result != 0 ) return(false);
            // calculate weights
            vout << "   Calculating weights B ..." << endl;
            CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
            KSInverted = true;
            }
            break;
    default:
        INVALID_ARGUMENT("unsupported method");
    }

    return(true);
}

//------------------------------------------------------------------------------

const CSmallString CGHSIntegratorGPR0A::GetKernelName(void)
{
    switch(Kernel){
    case(EGPRK_ARDSE):
        return("ARD squared exponential");
    case(EGPRK_ARDMC52):
        return("ARD Matern class 5/2");
    default:
        RUNTIME_ERROR("not implemented");
    }
}

////------------------------------------------------------------------------------

// TK is multitask covariance matrix
void CGHSIntegratorGPR0A::CreateTK(void)
{
    TK[0][0] = SigmaF2[0];
    TK[1][1] = SigmaF2[1];
    TK[2][2] = SigmaF2[2];

    TK[0][1] = CoVar[0];
    TK[1][0] = CoVar[0];

    TK[0][2] = CoVar[1];
    TK[2][0] = CoVar[1];

    TK[1][2] = CoVar[2];
    TK[2][1] = CoVar[2];

    if( ConstrainedTK ){
        // impose dG/dx - dH/dx - (-TdS/dx) = 0 constraint
        CFortranMatrix A;
        A.CreateMatrix(3,3);
        CFortranMatrix B;
        B.CreateMatrix(3,3);

        A[0][0] =  0.0;
        A[0][1] = -1.0;
        A[0][2] =  1.0;

        A[1][0] =  1.0;
        A[1][1] =  0.0;
        A[1][2] =  1.0;

        A[2][0] = -1.0;
        A[2][1] = -1.0;
        A[2][2] =  0.0;

        CSciBlas::gemm(1.0,A,TK,0.0,B);
        CSciBlas::gemm(1.0,'N',B,'T',A,0.0,TK);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CreateKS(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    // generate main kernel block
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];

        GSurface->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];

            GSurface->GetPoint(jbin,jpos);

            // calc Kblock
            if( UseNumDiff ){
                GetKernelDer2Num(ipos,jpos,kblock);
            } else {
                GetKernelDer2Ana(ipos,jpos,kblock);
            }

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    KS[offset0 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TK[0][0]*kblock[ii][jj];
                    KS[offset0 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TK[0][1]*kblock[ii][jj];
                    KS[offset0 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TK[0][2]*kblock[ii][jj];

                    KS[offset1 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TK[1][0]*kblock[ii][jj];
                    KS[offset1 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TK[1][1]*kblock[ii][jj];
                    KS[offset1 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TK[1][2]*kblock[ii][jj];

                    KS[offset2 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TK[2][0]*kblock[ii][jj];
                    KS[offset2 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TK[2][1]*kblock[ii][jj];
                    KS[offset2 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TK[2][2]*kblock[ii][jj];
                }
            }
        }
    }

// error of data points;
    int offsetn0 = 0*NumOfCVs;
    int offsetn1 = 1*NumOfCVs;
    int offsetn2 = 2*NumOfCVs;

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        for(size_t ii=0; ii < NumOfCVs; ii++){
            KS[offset0 + indi*NumOfCVs+ii][offset0 + indi*NumOfCVs+ii] += SigmaN2[offsetn0 + ii];
            KS[offset1 + indi*NumOfCVs+ii][offset1 + indi*NumOfCVs+ii] += SigmaN2[offsetn1 + ii];
            KS[offset2 + indi*NumOfCVs+ii][offset2 + indi*NumOfCVs+ii] += SigmaN2[offsetn2 + ii];
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& kff,int task)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CSimpleVector<double> kder;
    kder.CreateVector(NumOfCVs);

    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    // main kernel matrix
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];

        GSurface->GetPoint(jbin,jpos);

        // calc Kder
        if( UseNumDiff ){
            GetKernelDerNumJ(ip,jpos,kder);
        } else {
            GetKernelDerAnaJ(ip,jpos,kder);
        }

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff[offset0 + indj*NumOfCVs+jj] = TK[0][task]*kder[jj];
            kff[offset1 + indj*NumOfCVs+jj] = TK[1][task]*kder[jj];
            kff[offset2 + indj*NumOfCVs+jj] = TK[2][task]*kder[jj];
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& kff2)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];

        GSurface->GetPoint(jbin,jpos);

        // calc Kder
        if( UseNumDiff ){
            GetKernelDer2Num(ip,jpos,kblock);
        } else {
            GetKernelDer2Ana(ip,jpos,kblock);
        }

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff2[indj*NumOfCVs+jj] = kblock[icoord][jj];
        }
    }

    int offseti;

    offseti = 1*NumOfUsedBins*NumOfCVs;
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff2[offseti + indj*NumOfCVs+jj] = kff2[indj*NumOfCVs+jj];
        }
    }

    offseti = 2*NumOfUsedBins*NumOfCVs;
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff2[offseti + indj*NumOfCVs+jj] = kff2[indj*NumOfCVs+jj];
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDerAnaI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDerNumI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NumOfCVs);

    double  dh = 1e-3;
    double  v1,v2;

    // off diagonal elements
    for(size_t ii=0; ii < NumOfCVs; ii++) {

        tip = jp;
        tip[ii] = ip[ii] - dh;
        v1 = GetKernelValue(tip,jp);

        tip = jp;
        tip[ii] = ip[ii] + dh;
        v2 = GetKernelValue(tip,jp);

        kder[ii] = (v2 - v1)/(2.0*dh);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDerAnaJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDerNumJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    CSimpleVector<double> tjp;
    tjp.CreateVector(NumOfCVs);

    double  dh = 1e-3;
    double  v1,v2;

    // off diagonal elements
    for(size_t ii=0; ii < NumOfCVs; ii++) {

        tjp = jp;
        tjp[ii] = jp[ii] - dh;
        v1 = GetKernelValue(ip,tjp);

        tjp = jp;
        tjp[ii] = jp[ii] + dh;
        v2 = GetKernelValue(ip,tjp);

        kder[ii] = (v2 - v1)/(2.0*dh);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDer2Ana(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    kblock.SetZero();

    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double pre = exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                GSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] -= pre*du/dd;
                    if( ii == jj ){
                        kblock[ii][ii] += pre/CVLengths2[ii];
                    }
                }
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pr = -(5.0/3.0)*exp(-sqrt(5.0)*scdist);
            double d1 = pr*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                GSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] += 5.0*pr*du/dd;
                    if( (ii == jj) ){
                        kblock[ii][jj] -= d1/(CVLengths2[ii]);
                    }
                }
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDer2Num(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NumOfCVs);
    CSimpleVector<double> tjp;
    tjp.CreateVector(NumOfCVs);

    double  dh = 1e-4;
    double  v1,v2,v3,v4;

    // off diagonal elements
    for(size_t ii=0; ii < NumOfCVs; ii++) {
        for(size_t jj=0; jj < NumOfCVs; jj++) {

            tip = ip;
            tip[ii] = ip[ii] + dh;
            tjp = jp;
            tjp[jj] = jp[jj] + dh;
            v1 = GetKernelValue(tip,tjp);

            tip = ip;
            tip[ii] = ip[ii] + dh;
            tjp = jp;
            tjp[jj] = jp[jj] - dh;
            v2 = GetKernelValue(tip,tjp);

            tip = ip;
            tip[ii] = ip[ii] - dh;
            tjp = jp;
            tjp[jj] = jp[jj] + dh;
            v3 = GetKernelValue(tip,tjp);

            tip = ip;
            tip[ii] = ip[ii] - dh;
            tjp = jp;
            tjp[jj] = jp[jj] - dh;
            v4 = GetKernelValue(tip,tjp);

            kblock[ii][jj] = (v1 - v2 - v3 + v4)/(4.0*dh*dh);
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetKernelDer2AnaWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    kblock.SetZero();

    double wf = WFac[cv];
    double wd3 = 1.0/(CVLengths2[cv]*wf);
    double wd5 = wd3/CVLengths2[cv];
    double dc = GSurface->GetCV(cv)->GetDifference(ip[cv],jp[cv]);

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double arg = exp(-0.5*scdist2);
            double argd = arg*dc*dc*wd3;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                GSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] -= argd*du/dd;
                    if( (cv == ii) && (cv != jj) ) {
                        kblock[ii][jj] += 2.0*arg*du*wd3/(CVLengths2[jj]);
                    } else if( (cv == jj) && (cv != ii) ) {
                        kblock[ii][jj] += 2.0*arg*du*wd3/(CVLengths2[ii]);
                    } else if( (cv == ii) && (cv == jj) ) {
                        kblock[ii][jj] += 4.0*arg*du*wd5;
                    }
                    if( (ii == jj) ){
                        kblock[ii][ii] += argd/CVLengths2[ii];
                        if( ii == cv ){
                            kblock[ii][ii] -= 2.0*arg*wd3;
                        }
                    }
                }
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pr = -(5.0/3.0)*exp(-sqrt(5.0)*scdist);
            double prd = 0.0;
            if( scdist > 0 ){
               prd = -(5.0/3.0)*sqrt(5.0)*exp(-sqrt(5.0)*scdist)*(1.0/scdist)*wd3*dc*dc;
            }
            double d1 = pr*(sqrt(5.0)*scdist+1.0);
            double d1d = -(25.0/3.0)*exp(-sqrt(5.0)*scdist)*wd3*dc*dc;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                GSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] += 5.0*prd*du/dd;
                    if( (cv == ii) && (cv != jj) ) {
                        kblock[ii][jj] -= 10.0*pr*du*wd3/(CVLengths2[jj]);
                    } else if( (cv == jj) && (cv != ii) ) {
                        kblock[ii][jj] -= 10.0*pr*du*wd3/(CVLengths2[ii]);
                    } else if( (cv == ii) && (cv == jj) ) {
                        kblock[ii][jj] -= 20.0*pr*du*wd5;
                    }
                    if( (ii == jj) ){
                        kblock[ii][ii] -= d1d/(CVLengths2[ii]);
                        if( ii == cv ){
                            kblock[ii][ii] += 2.0*d1*wd3;
                        }

                    }
                }
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

// for internal debugging
void CGHSIntegratorGPR0A::GetKernelDer2NumWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    CFortranMatrix kblock1,kblock2;

    kblock1.CreateMatrix(NumOfCVs,NumOfCVs);
    kblock2.CreateMatrix(NumOfCVs,NumOfCVs);

    double  dh = 1e-3;

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*GSurface->GetCV(i)->GetRange()/GSurface->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*GSurface->GetCV(i)->GetRange()/GSurface->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDer2Ana(ip,jp,kblock1);

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*GSurface->GetCV(i)->GetRange()/GSurface->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*GSurface->GetCV(i)->GetRange()/GSurface->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDer2Ana(ip,jp,kblock2);

    for(size_t ii=0; ii < NumOfCVs; ii++){
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kblock[ii][jj] = (kblock2[ii][jj]-kblock1[ii][jj])/(2.0*dh);
        }
    }

    // restore original CVLengths2
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*GSurface->GetCV(i)->GetRange()/GSurface->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPR0A::GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
{
    // get kernel value
    switch(Kernel){
        case(EGPRK_ARDSE): {
            // calculate scaled distance
            double scdist2 = 0.0;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                scdist2 += du*du/dd;
            }
            return(exp(-0.5*scdist2));
        }
        break;
        case(EGPRK_ARDMC52):{
            // calculate scaled distance
            double scdist2 = 0.0;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = GSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                scdist2 += du*du/dd;
            }
            double scdist = sqrt(scdist2);
            return((1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2)*exp(-sqrt(5.0)*scdist));
        }
        break;
        default:
            RUNTIME_ERROR("not implemented");
        break;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPR0A::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating dG, dH, -TdS surfaces ..." << endl;
    vout << "   >>>>>>>>" << endl;

// calculate energies
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

// basic EneSurface update
    #pragma omp parallel for
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        int nsamples = HDerProxy->GetNumOfSamples(ibin);
        GSurface->SetNumOfSamples(ibin,nsamples);
        GSurface->SetEnergy(ibin,0.0);
        GSurface->SetError(ibin,0.0);
        HSurface->SetNumOfSamples(ibin,nsamples);
        HSurface->SetEnergy(ibin,0.0);
        HSurface->SetError(ibin,0.0);
        SSurface->SetNumOfSamples(ibin,nsamples);
        SSurface->SetEnergy(ibin,0.0);
        SSurface->SetError(ibin,0.0);
    }

// update FES
    #pragma omp parallel for firstprivate(jpos)
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        int jbin = SampledMap[indj];
        GSurface->GetPoint(jbin,jpos);
        double dg,dh,mtds;
        GetValues(jpos,dg,dh,mtds);
        GSurface->SetEnergy(jbin,dg);
        HSurface->SetEnergy(jbin,dh);
        SSurface->SetEnergy(jbin,mtds);
    }

// update FES
    if( GlobalMinSet ){
        // GPos.CreateVector(NumOfCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "      dG(x) Global minimum provided at: ";
        vout << setprecision(5) << GSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << GSurface->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << endl;

        // find the closest bin
        CSimpleVector<double>   pos;
        pos.CreateVector(NumOfCVs);
        double minv = 0.0;
        GPosBin = 0;
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            GSurface->GetPoint(ibin,pos);
            double dist2 = 0.0;
            for(size_t cv=0; cv < NumOfCVs; cv++){
                dist2 = dist2 + (pos[cv]-GPos[cv])*(pos[cv]-GPos[cv]);
            }
            if( ibin == 0 ){
                minv = dist2;
                GPosBin = 0;
            }
            if( dist2 < minv ){
                minv = dist2;
                GPosBin = ibin;
            }
        }

        GSurface->GetPoint(GPosBin,GPos);
        GPosSet = true;

        vout << "      dG(x) Closest bin found at: ";
        vout << setprecision(5) << GSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << GSurface->GetCV(i)->GetRealValue(GPos[i]);
        }

        double glb_min = GSurface->GetEnergy(GPosBin);
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            int jbin = SampledMap[indj];
            GSurface->SetEnergy(jbin,GSurface->GetEnergy(jbin)-glb_min);
        }
    } else {
        // search for global minimum
        GPos.CreateVector(NumOfCVs);
        bool   first = true;
        double glb_min = 0.0;
        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            int jbin = SampledMap[indj];
            int samples = GSurface->GetNumOfSamples(jbin);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = GSurface->GetEnergy(jbin);
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                GPosBin = jbin;
                GSurface->GetPoint(jbin,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      dG(x) Global minimum found at: ";
        vout << setprecision(5) << GSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << GSurface->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            int jbin = SampledMap[indj];
            GSurface->SetEnergy(jbin,GSurface->GetEnergy(jbin)-glb_min);
        }
    }

    GPosSet = true;

    vout << "      dG(x) SigmaF2   = " << setprecision(5) << GSurface->GetSigmaF2() << endl;
    vout << "      dG(x) SigmaF    = " << setprecision(5) << GSurface->GetSigmaF() << endl;

    vout << "   >>>>>>>>" << endl;
    double glb_min = HSurface->GetEnergy(GPosBin);
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        int jbin = SampledMap[indj];
        HSurface->SetEnergy(jbin,HSurface->GetEnergy(jbin)-glb_min);
    }
    vout << "      dH(x) SigmaF2   = " << setprecision(5) << HSurface->GetSigmaF2() << endl;
    vout << "      dH(x) SigmaF    = " << setprecision(5) << HSurface->GetSigmaF() << endl;

    vout << "   >>>>>>>>" << endl;
    glb_min = SSurface->GetEnergy(GPosBin);
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        int jbin = SampledMap[indj];
        SSurface->SetEnergy(jbin,SSurface->GetEnergy(jbin)-glb_min);
    }
    vout << "      -TdS(x) SigmaF2 = " << setprecision(5) << SSurface->GetSigmaF2() << endl;
    vout << "      -TdS(x) SigmaF  = " << setprecision(5) << SSurface->GetSigmaF() << endl;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetValues(const CSimpleVector<double>& position,double& dg, double& dh, double& mtds)
{
    CSimpleVector<double>   kff;
    kff.CreateVector(GPRSize);

    RunBlasLapackSeq();

// FIXME
    CreateKff(position,kff,0);
    dg   = CSciBlas::dot(kff,GPRModel);

    CreateKff(position,kff,1);
    dh   = CSciBlas::dot(kff,GPRModel);

    CreateKff(position,kff,2);
    mtds = CSciBlas::dot(kff,GPRModel);

    //cout << dg << " " << dh << " " << mtds << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGHSIntegratorGPR0A::GetLogML(void)
{
    double ml = 0.0;

    // http://www.gaussianprocess.org/gpml/chapters/RW5.pdf
    // page 113

    RunBlasLapackPar();

    ml -= CSciBlas::dot(Y,GPRModel);
    ml -= logdetK;
    ml -= GPRSize * log(2*M_PI);
    ml *= 0.5;

    return(ml);
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPR0A::GetLogPL(void)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("logPL requires K+Sigma inverted matrix");
    }

    double loo = 0.0;

    // http://www.gaussianprocess.org/gpml/chapters/RW5.pdf
    // page 116-117

    for(size_t i=0; i < GPRSize; i++){
        loo -= log(1.0/KS[i][i]);
        loo -= GPRModel[i]*GPRModel[i]/KS[i][i];
    }
    loo -= GPRSize*log(2*M_PI);

    loo *= 0.5;

    return(loo);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


