// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <IntegratorGPR.hpp>
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

CIntegratorGPR::CIntegratorGPR(void)
{
    GPRSize             = 0;
    NumOfUsedBins       = 0;
    NumOfValues         = 0;
    NumOfCVs            = 0;
    NumOfBins           = 0;

    SigmaF2             = 15.0;
    NCorr               = 1.0;

    IncludeError        = false;
    IncludeGluedBins    = false;
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
    UseZeroPoint        = false;
    FastErrors          = true;

    KSInverted          = false;
}

//------------------------------------------------------------------------------

CIntegratorGPR::~CIntegratorGPR(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CIntegratorGPR::AddInputEnergyDerProxy(CEnergyDerProxyPtr p_proxy)
{
    if( p_proxy == NULL ) return;                 // no-proxy
    if( p_proxy->GetAccu() == NULL ) return;      // no PMFAccu

    if( NumOfCVs != (size_t)p_proxy->GetAccu()->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_proxy->GetAccu()->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    DerProxyItems.push_back(p_proxy);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::ClearInputEnergyDerProxies(void)
{
    DerProxyItems.clear();
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetOutputES(CEnergySurfacePtr p_surf)
{
    NumOfCVs = 0;
    NumOfBins = 0;
    EneSurface = p_surf;

    if( EneSurface != NULL ){
        NumOfCVs = EneSurface->GetNumOfCVs();
        NumOfBins = EneSurface->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetSigmaF2(double sigf2)
{
    SigmaF2 = sigf2;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetNCorr(double value)
{
    NCorr = value;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetSigmaN2(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetSigmaN2");
    }

    string          sspec(spec);
    vector<string>  ssigman2;

    split(ssigman2,sspec,is_any_of("x"),token_compress_on);

    if( ssigman2.size() > NumOfCVs ){
        CSmallString error;
        error << "too many sigman2 (" << ssigman2.size() << ") than required (" << NumOfCVs << ")";
        RUNTIME_ERROR(error);
    }

    SigmaN2.CreateVector(NumOfCVs);

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
    for(size_t i=ssigman2.size(); i < NumOfCVs; i++){
        SigmaN2[i] = last_sigman2;
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetSigmaN2(CSimpleVector<double>& sigman2)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("FES is not set for SetSigmaN2");
    }
    if( sigman2.GetLength() != NumOfCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    SigmaN2 = sigman2;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetSigmaN2(size_t cvind, double value)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("FES is not set for SetSigmaN2");
    }
    if( cvind >= NumOfCVs ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is SigmaN2 initialized?
    if( SigmaN2.GetLength() == 0 ){
        SigmaN2.CreateVector(NumOfCVs);
    }
    SigmaN2[cvind] = value;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetWFac(const CSmallString& spec)
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

void CIntegratorGPR::SetWFac(CSimpleVector<double>& wfac)
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

void CIntegratorGPR::SetWFac(size_t cvind, double value)
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

void CIntegratorGPR::LoadGPRHyprms(const CSmallString& name)
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
        if( (key == "SigmaF2") ||(key == "SigmaF2#1") ){
            SetSigmaF2(value);
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
        } else if( (key == "NCorr") || (key == "NCorr#1") ){
            SetNCorr(value);
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

void CIntegratorGPR::SetIncludeError(bool set)
{
    IncludeError = set;
    if( FastErrors == false ){
        NeedInv |= set;
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetUseNumDiff(bool set)
{
    UseNumDiff = set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetNoEnergy(bool set)
{
    NoEnergy = set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetLAMethod(EGPRLAMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetLAMethod(const CSmallString& method)
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

void CIntegratorGPR::SetKernel(const CSmallString& kernel)
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

void CIntegratorGPR::IncludeGluedAreas(bool set)
{
    IncludeGluedBins = set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetGlobalMin(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetGlobalMin");
    }
    if( EneSurface == NULL ) {
        RUNTIME_ERROR("ES is not set");
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
        GPos[i] = EneSurface->GetCV(i)->GetIntValue(val);
    }

    GPosSet = true;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetGlobalMin(const CSimpleVector<double>& pos)
{
    GlobalMinSet = true;
    GPos = pos;
    GPosSet = true;
}

//------------------------------------------------------------------------------

CSimpleVector<double> CIntegratorGPR::GetGlobalMin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPos);
}

//------------------------------------------------------------------------------

int CIntegratorGPR::GetGlobalMinBin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPosBin);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetUseInv(bool set)
{
    UseInv = set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::PrepForHyprmsGrd(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetCalcLogPL(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetUseZeroPoint(bool set)
{
    if( set == true ){
        if( GlobalMinSet == false ){
            RUNTIME_ERROR("position of global minimum was not set");
        }
    }
    UseZeroPoint = set;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetFastError(bool set)
{
    FastErrors = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CIntegratorGPR::Integrate(CVerboseStr& vout,bool nostat)
{
    PrintExecInfo(vout);

    if( DerProxyItems.size() == 0 ){
        RUNTIME_ERROR("no input data");
    }
    if( EneSurface == NULL ) {
        RUNTIME_ERROR("ES is not set");
    }
    if( NumOfCVs == 0 ) {
        RUNTIME_ERROR("number of CVs is zero");
    }
    if( NumOfBins == 0 ) {
        RUNTIME_ERROR("number of bins is zero");
    }
    if( WFac.GetLength() == 0 ){
        RUNTIME_ERROR("wfac is not set");
    }
    if( SigmaN2.GetLength() == 0 ){
        RUNTIME_ERROR("sigman2 is not set");
    }

    // GPR setup
    CVLengths2.CreateVector(NumOfCVs);
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }

    // number of data points
    NumOfUsedBins = 0;
    for(size_t i=0; i < DerProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( DerProxyItems[i]->GetNumOfSamples(ibin) > 0 ) NumOfUsedBins++;
        }
    }
    GPRSize = NumOfUsedBins * NumOfCVs;
    if( UseZeroPoint ) GPRSize++;

    // create sampled map
    SampledMap.resize(NumOfUsedBins);
    DerProxyMap.resize(NumOfUsedBins);
    size_t ind = 0;
    for(size_t i=0; i < DerProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( DerProxyItems[i]->GetNumOfSamples(ibin) <= 0 ) continue;
            SampledMap[ind] = ibin;
            DerProxyMap[ind] = i;
            ind++;
        }
    }

    // init GPR arrays
    GPRModel.CreateVector(GPRSize);
    Y.CreateVector(GPRSize);
    KS.CreateMatrix(GPRSize,GPRSize);

    // print hyperparameters
        vout        << "   Hyperparameters ..." << endl;
        vout << format("      SigmaF2   = %10.4f")%SigmaF2 << endl;
        vout << format("      NCorr     = %10.4f")%NCorr << endl;

    for(size_t k=0; k < NumOfCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(size_t k=0; k < NumOfCVs; k++ ){
        vout << format("      SigmaN2#%-2d= %10.4e")%(k+1)%SigmaN2[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    if( ! nostat ){
        // and finally some statistics
        for(size_t k=0; k < NumOfCVs; k++ ){
            vout << "      RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
        }

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
                CalculateCovs(vout);
                CalculateErrorsFromCov(vout);
            } else {
                CalculateErrors(GPos,vout);
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::PrintExecInfo(CVerboseStr& vout)
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

void CIntegratorGPR::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CIntegratorGPR::TrainGP(CVerboseStr& vout)
{
    if( UseZeroPoint ){
        vout << "   Zero-point included in GPR ..." << endl;
    }

    if( UseNumDiff ) {
        vout << "   Creating K+Sigma and Y (numeric differentation) ..." << endl;
    } else {
        vout << "   Creating K+Sigma and Y ..." << endl;
    }
        vout << "      Kernel    = " << GetKernelName() << endl;
        vout << "      Dim       = " << GPRSize << " x " << GPRSize << endl;

// construct Y
    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];
        for(size_t ii=0; ii < NumOfCVs; ii++){
            double mf = item->GetValue(ibin,ii,E_PROXY_VALUE);
            Y[indi*NumOfCVs+ii] = mf;
        }
    }
    if( UseZeroPoint ){
        Y[GPRSize-1] = 0.0;
    }

// construct KS
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

const CSmallString CIntegratorGPR::GetKernelName(void)
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

//------------------------------------------------------------------------------

void CIntegratorGPR::CreateKS(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    // main kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];

            EneSurface->GetPoint(jbin,jpos);

            // calc Kblock
            if( UseNumDiff ){
                GetKernelDer2Num(ipos,jpos,kblock);
            } else {
                GetKernelDer2Ana(ipos,jpos,kblock);
            }

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    KS[indi*NumOfCVs+ii][indj*NumOfCVs+jj] = kblock[ii][jj];
                }
            }
        }
    }
    if( UseZeroPoint ){
        CSimpleVector<double> kder;
        kder.CreateVector(NumOfCVs);

        #pragma omp parallel for firstprivate(ipos,kder)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t ibin = SampledMap[indi];

            EneSurface->GetPoint(ibin,ipos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumI(ipos,GPos,kder);
            } else {
                GetKernelDerAnaI(ipos,GPos,kder);
            }

            // distribute to matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                KS[indi*NumOfCVs+ii][GPRSize-1] = kder[ii];
            }
        }

        #pragma omp parallel for firstprivate(jpos,kder)
        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];

            EneSurface->GetPoint(jbin,jpos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumJ(GPos,jpos,kder);
            } else {
                GetKernelDerAnaJ(GPos,jpos,kder);
            }

            // distribute to matrix
            for(size_t jj=0; jj < NumOfCVs; jj++){
                KS[GPRSize-1][indj*NumOfCVs+jj] = kder[jj];
            }
        }

        KS[GPRSize-1][GPRSize-1] = GetKernelValue(GPos,GPos);
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];
        for(size_t ii=0; ii < NumOfCVs; ii++){
            double er = item->GetValue(ibin,ii,E_PROXY_ERROR);
            KS[indi*NumOfCVs+ii][indi*NumOfCVs+ii] += er*er*NCorr + SigmaN2[ii];
        }
    }
    // zero point has no uncertainty
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& kff)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CSimpleVector<double> kder;
    kder.CreateVector(NumOfCVs);

    // main kernel matrix
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];

        EneSurface->GetPoint(jbin,jpos);

        // calc Kder
        if( UseNumDiff ){
            GetKernelDerNumJ(ip,jpos,kder);
        } else {
            GetKernelDerAnaJ(ip,jpos,kder);
        }

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff[indj*NumOfCVs+jj] = kder[jj];
        }
    }

    if( UseZeroPoint ){
        kff[GPRSize-1] = GetKernelValue(ip,GPos);
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& kff2)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];

        EneSurface->GetPoint(jbin,jpos);

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

    if( UseZeroPoint ){
        CSimpleVector<double> kder;
        kder.CreateVector(NumOfCVs);

        // calc Kder
        if( UseNumDiff ){
            GetKernelDerNumI(ip,GPos,kder);
        } else {
            GetKernelDerAnaI(ip,GPos,kder);
        }

        kff2[GPRSize-1] = kder[icoord];
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::GetKernelDerAnaI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = SigmaF2*exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*SigmaF2*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
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

void CIntegratorGPR::GetKernelDerNumI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
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

void CIntegratorGPR::GetKernelDerAnaJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = SigmaF2*exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*SigmaF2*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
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

void CIntegratorGPR::GetKernelDerNumJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
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

void CIntegratorGPR::GetKernelDer2Ana(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    kblock.SetZero();

    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double pre = SigmaF2*exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                EneSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
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
            double pr = -SigmaF2*(5.0/3.0)*exp(-sqrt(5.0)*scdist);
            double d1 = pr*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                EneSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
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

void CIntegratorGPR::GetKernelDer2Num(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
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

void CIntegratorGPR::GetKernelDer2AnaWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    kblock.SetZero();

    double wf = WFac[cv];
    double wd3 = 1.0/(CVLengths2[cv]*wf);
    double wd5 = wd3/CVLengths2[cv];
    double dc = EneSurface->GetCV(cv)->GetDifference(ip[cv],jp[cv]);

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double arg = SigmaF2*exp(-0.5*scdist2);
            double argd = arg*dc*dc*wd3;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                EneSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
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
            double pr = -SigmaF2*(5.0/3.0)*exp(-sqrt(5.0)*scdist);
            double prd = 0.0;
            if( scdist > 0 ){
               prd = -SigmaF2*(5.0/3.0)*sqrt(5.0)*exp(-sqrt(5.0)*scdist)*(1.0/scdist)*wd3*dc*dc;
            }
            double d1 = pr*(sqrt(5.0)*scdist+1.0);
            double d1d = -(25.0/3.0)*SigmaF2*exp(-sqrt(5.0)*scdist)*wd3*dc*dc;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                EneSurface->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
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
void CIntegratorGPR::GetKernelDer2NumWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    CFortranMatrix kblock1,kblock2;

    kblock1.CreateMatrix(NumOfCVs,NumOfCVs);
    kblock2.CreateMatrix(NumOfCVs,NumOfCVs);

    double  dh = 1e-3;

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDer2Ana(ip,jp,kblock1);

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
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
        double l = WFac[i]*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }
}

//------------------------------------------------------------------------------

double CIntegratorGPR::GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
{
    // get kernel value
    switch(Kernel){
        case(EGPRK_ARDSE): {
            // calculate scaled distance
            double scdist2 = 0.0;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                scdist2 += du*du/dd;
            }
            return(SigmaF2*exp(-0.5*scdist2));
        }
        break;
        case(EGPRK_ARDMC52):{
            // calculate scaled distance
            double scdist2 = 0.0;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                scdist2 += du*du/dd;
            }
            double scdist = sqrt(scdist2);
            return(SigmaF2*(1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2)*exp(-sqrt(5.0)*scdist));
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

void CIntegratorGPR::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating EneSurface ..." << endl;

// create map for bins with calculated energy and error
    std::set<size_t>    vset;
    for(size_t i=0; i < DerProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            int samples = DerProxyItems[i]->GetNumOfSamples(ibin);
            if( IncludeGluedBins ){
                if( samples == 0 ) continue;
            } else {
                if( samples <= 0 ) continue;
            }
            vset.insert(ibin);
        }
    }
    NumOfValues = vset.size();
    ValueMap.resize(NumOfValues);

    size_t indi = 0;
    std::set<size_t>::iterator  it = vset.begin();
    std::set<size_t>::iterator  ie = vset.end();

    while( it != ie ){
        ValueMap[indi] = *it;
        it++;
        indi++;
    }

// calculate energies
    CSimpleVector<double> jpos;
    CSimpleVector<double> values;

    jpos.CreateVector(NumOfCVs);
    values.CreateVector(NumOfValues);

    #pragma omp parallel for firstprivate(jpos)
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        EneSurface->GetPoint(j,jpos);
        values[indj] = GetValue(jpos);
    }

// basic EneSurface update
    for(size_t i=0; i < DerProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            int nsamples = DerProxyItems[i]->GetNumOfSamples(ibin);
            int osamples = EneSurface->GetNumOfSamples(ibin);
            EneSurface->SetNumOfSamples(ibin,nsamples+osamples);
            EneSurface->SetEnergy(ibin,0.0);
            EneSurface->SetError(ibin,0.0);
        }
    }

// update FES
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        EneSurface->SetEnergy(j,values[indj]);
    }

// update FES
    if( GlobalMinSet ){
        // GPos.CreateVector(NumOfCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << EneSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << EneSurface->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << endl;

        // find the closest bin
        CSimpleVector<double>   pos;
        pos.CreateVector(NumOfCVs);
        double minv = 0.0;
        GPosBin = 0;
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            EneSurface->GetPoint(ibin,pos);
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

        EneSurface->GetPoint(GPosBin,GPos);
        GPosSet = true;

        vout << "      Closest bin found at: ";
        vout << setprecision(5) << EneSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << EneSurface->GetCV(i)->GetRealValue(GPos[i]);
        }

        double glb_min = EneSurface->GetEnergy(GPosBin);
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->SetEnergy(j,EneSurface->GetEnergy(j)-glb_min);
        }
    } else {
        // search for global minimum
        GPos.CreateVector(NumOfCVs);
        bool   first = true;
        double glb_min = 0.0;
        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            int samples = EneSurface->GetNumOfSamples(j);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = values[indj];
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                GPosBin = j;
                EneSurface->GetPoint(j,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum found at: ";
        vout << setprecision(5) << EneSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << EneSurface->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->SetEnergy(j,values[indj]-glb_min);
        }
    }

    GPosSet = true;

        vout << "      SigmaF2   = " << setprecision(5) << EneSurface->GetSigmaF2() << endl;
    if( IncludeGluedBins ){
        vout << "      SigmaF2 (including glued bins) = " << setprecision(5) << EneSurface->GetSigmaF2(true) << endl;
    }
        vout << "      SigmaF    = " << setprecision(5) << EneSurface->GetSigmaF() << endl;
}

//------------------------------------------------------------------------------

double CIntegratorGPR::GetValue(const CSimpleVector<double>& position)
{
    CSimpleVector<double>   kff;
    kff.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(position,kff);
    double energy = CSciBlas::dot(kff,GPRModel);
    return(energy);
}

//------------------------------------------------------------------------------

double CIntegratorGPR::GetMeanForce(const CSimpleVector<double>& position,size_t icoord)
{
    CSimpleVector<double>   kff2;
    kff2.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff2(position,icoord,kff2);
    double mf = CSciBlas::dot(kff2,GPRModel);
    return(mf);
}

//------------------------------------------------------------------------------

double CIntegratorGPR::GetMeanForceVar(const CSimpleVector<double>& position,size_t icoord)
{
    if( KSInverted != true ) {
        RUNTIME_ERROR("KS must be inverted!");
    }

    CSimpleVector<double>   kff2;
    CSimpleVector<double>   ik;

    kff2.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff2(position,icoord,kff2);
    CSciBlas::gemv(1.0,KS,kff2,0.0,ik);     // KS must be inverted

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    GetKernelDer2Ana(position,position,kblock);

    // FIXME - valiadate
    double var = kblock[icoord][icoord] - CSciBlas::dot(kff2,ik);
    return(var);
}

//------------------------------------------------------------------------------

double CIntegratorGPR::GetRMSR(size_t cv)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

    double rmsr = 0.0;

    #pragma omp parallel for firstprivate(ipos) reduction(+:rmsr)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);

        double mfi = item->GetValue(ibin,cv,E_PROXY_VALUE);
        double mfp = GetMeanForce(ipos,cv);
        double diff = mfi - mfp;
        rmsr += diff*diff;
    }

    double nsamples = NumOfUsedBins;

    if( nsamples > 0 ){
        rmsr /= nsamples;
    }
    if( rmsr > 0.0 ){
        rmsr = sqrt(rmsr);
    }

    return(rmsr);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::PrepForMFInfo(void)
{
    NeedInv = true;
}

//------------------------------------------------------------------------------

bool CIntegratorGPR::WriteMFInfo(const CSmallString& name)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> ipos;
    CSimpleVector<double> mfi;
    CSimpleVector<double> mfie;
    CSimpleVector<double> mfp;
    CSimpleVector<double> mfpe;

    ipos.CreateVector(NumOfCVs);
    mfi.CreateVector(GPRSize);
    mfie.CreateVector(GPRSize);
    mfp.CreateVector(GPRSize);
    mfpe.CreateVector(GPRSize);

    // calculate
    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);
        for(size_t k=0; k < NumOfCVs; k++){
            mfi[indi*NumOfCVs+k] = item->GetValue(ibin,k,E_PROXY_VALUE);
            double mfe = item->GetValue(ibin,k,E_PROXY_ERROR);
            mfie[indi*NumOfCVs+k] = mfe;            // this is a sigma

            mfp[indi*NumOfCVs+k] = GetMeanForce(ipos,k);
            double mfv = GetMeanForceVar(ipos,k);   // this is a variance, sigma^2
            mfpe[indi*NumOfCVs+k] = sqrt(mfv);      // get a sigma
        }
    }

    ofstream ofs(name);
    if( ! ofs ){
        CSmallString error;
        error << "unable to open file '" << name << "' for derivatives";
        ES_ERROR(error);
        return(false);
    }

    // print
    size_t accu_prev = DerProxyMap[0];
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        size_t accu = DerProxyMap[indi];

        if( accu != accu_prev ){
            ofs << endl;
            accu_prev = accu;
        }

        EneSurface->GetPoint(ibin,ipos);

        for(size_t c=0; c < NumOfCVs; c++){
            ofs << format("%20.16f ")%ipos[c];
        }

        for(size_t k=0; k < NumOfCVs; k++){
            ofs << format(" %20.16f %20.16f %20.16f %20.16f")%mfi[indi*NumOfCVs+k]%mfie[indi*NumOfCVs+k]%mfp[indi*NumOfCVs+k]%mfpe[indi*NumOfCVs+k];
        }

        ofs << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::FilterByMFZScore(double zscore,CVerboseStr& vout)
{
    vout << high;
    vout << "   Running Z-test filter on MF errors ..." << endl;

    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return;
    }

    // we work with variances
    zscore *= zscore;

    CSimpleVector<int>      flags;
    CSimpleVector<double>   sig2;
    CSimpleVector<int>      maxi;
    CSimpleVector<double>   maxzscore;
    CSimpleVector<double>   mferror2;

    flags.CreateVector(NumOfBins);
    sig2.CreateVector(NumOfCVs);
    maxi.CreateVector(NumOfCVs);
    maxzscore.CreateVector(NumOfCVs);
    mferror2.CreateVector(GPRSize);

    flags.Set(1);

    vout << high;
    vout << "      Pre-calculating MF errors ..." << endl;

    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

    // precalculate values
    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);

        for(size_t k=0; k < NumOfCVs; k++){
            double mf = item->GetValue(ibin,k,E_PROXY_VALUE);
            double diff2 = mf - GetMeanForce(ipos,k);
            diff2 *= diff2;
            mferror2[indi*NumOfCVs+k] = diff2;
        }
    }

    vout << "      Searching for MF outliers ..." << endl;
    vout << debug;

    bool testme = true;
    while( testme ){

        testme = false;

        double  count = 0;
        sig2.SetZero();

        // calc variances - we assume zero mean on errors
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            if( flags[i] != 0 ) {
                for(size_t k=0; k < NumOfCVs; k++){
                    sig2[k] += mferror2[indi*NumOfCVs+k];
                }
                count++;
            }
        }

        if( count == 0 ) break;
        for(size_t k=0; k < NumOfCVs; k++){
            sig2[k] /= count;
        }

        // filter
        bool first = true;
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            if( flags[i] != 0 ){

                for(size_t k=0; k < NumOfCVs; k++){
                    double diff2 = mferror2[indi*NumOfCVs+k];
                    double zscore2 = diff2 / sig2[k];
                    if( first == true ){
                        maxzscore[k] = zscore2;
                        maxi[k] = i;
                    } else {
                        if( maxzscore[k] < zscore2 ) {
                            maxi[k] = i;
                            maxzscore[k] = zscore2;

                        }
                    }
                }

                first = false;
            }
        }

        for(size_t k=0; k < NumOfCVs; k++){
            if( maxzscore[k] > zscore ){
                flags[maxi[k]] = 0;
                testme = true;
                vout << "   outlier found at " << maxi[k] << " for cv " << k << " with z-score " << sqrt(maxzscore[k]) << endl;
            }
        }
    }

    // apply limits
    size_t outliers = 0;
    for(size_t i=0; i < NumOfBins; i++){
        if( flags[i] == 0 ){
            for(size_t j=0; j < DerProxyItems.size(); j++){
                DerProxyItems[j]->GetAccu()->SetNumOfSamples(i,0);
            }
            outliers++;
        }
    }

    vout << high;
    vout << "      Number of outliers = " << outliers << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CIntegratorGPR::GetLogML(void)
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

double CIntegratorGPR::GetLogPL(void)
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

//------------------------------------------------------------------------------

void CIntegratorGPR::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("GetLogMLDerivatives requires K+Sigma inverted matrix");
    }

    if( EneSurface == NULL ) {
        RUNTIME_ERROR("ES is not set");
    }
    if( NumOfCVs <= 0 ){
        RUNTIME_ERROR("NumOfCVs <= NULL");
    }
    if( GPRSize <= 0 ){
        RUNTIME_ERROR("GPRSize <= NULL");
    }

    Kder.CreateMatrix(GPRSize,GPRSize);

    CFortranMatrix          ata;            // alphaT*alpha
    ata.CreateMatrix(GPRSize,GPRSize);

    // calc ATA matrix
    #pragma omp parallel for
    for(size_t i=0; i < GPRSize; i++){
        for(size_t j=0; j < GPRSize; j++){
            ata[i][j] = GPRModel[i]*GPRModel[j];
        }
    }

    size_t ind = 0;

    for(size_t prm=0; prm < flags.size(); prm++){
        // shall we calc der?
        if( flags[prm] == false ) {
            continue;
        }

        // calc Kder
        // NumOfCVs = 1
        // 0; 1; 2; 3
        // 0; 1<1+NumOfCVs; NumOfCVs+1; 2+NumOfCVs<2+2*NumOfCVs+1
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( (prm >= 1) && (prm < 1+NumOfCVs) ){
            // wfac
            size_t cv = prm - 1;
            CalcKderWRTWFac(cv);
        } else if( prm == NumOfCVs+1 ){
            // ncorr
            CalcKderWRTNCorr();
        } else if( (prm >= 2+NumOfCVs) && (prm < 3+2*NumOfCVs) ){
            size_t cv = prm - (2+NumOfCVs);
            CalcKderWRTSigmaN2(cv);
        } else {
            RUNTIME_ERROR("prm out-of-range");
        }

        // calc trace
        double tr = 0.0;
        #pragma omp parallel for reduction(+:tr)
        for(size_t i=0; i < GPRSize; i++){
            for(size_t j=0; j < GPRSize; j++){
                tr += (ata[i][j]-KS[i][j])*Kder[j][i];
            }
        }

        // finalize derivative
        der[ind] += 0.5*tr;
        ind++;
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("GetLogPLDerivatives requires K+Sigma inverted matrix");
    }

    if( EneSurface == NULL ) {
        RUNTIME_ERROR("ES is not set");
    }
    if( NumOfCVs <= 0 ){
        RUNTIME_ERROR("NumOfCVs <= NULL");
    }
    if( GPRSize <= 0 ){
        RUNTIME_ERROR("GPRSize <= NULL");
    }

    Kder.CreateMatrix(GPRSize,GPRSize);

    CFortranMatrix zj;
    zj.CreateMatrix(GPRSize,GPRSize);

    CSimpleVector<double> za;
    za.CreateVector(GPRSize);

    size_t ind = 0;

    for(size_t prm=0; prm < flags.size(); prm++){
        // shall we calc der?
        if( flags[prm] == false ) {
            continue;
        }

        // calc Kder
        // NumOfCVs = 1
        // 0; 1; 2; 3
        // 0; 1<1+NumOfCVs; NumOfCVs+1; 2+NumOfCVs<2+2*NumOfCVs+1
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( (prm >= 1) && (prm < 1+NumOfCVs) ){
            // wfac
            size_t cv = prm - 1;
            CalcKderWRTWFac(cv);
        } else if( prm == NumOfCVs+1 ){
            // ncorr
            CalcKderWRTNCorr();
        } else if( (prm >= 2+NumOfCVs) && (prm < 3+2*NumOfCVs) ){
            size_t cv = prm - (2+NumOfCVs);
            CalcKderWRTSigmaN2(cv);
        } else {
            RUNTIME_ERROR("prm out-of-range");
        }

        RunBlasLapackPar();

        // calc Zj
        CSciBlas::gemm(1.0,KS,Kder,0.0,zj);

        // calc zj * alpha
        CSciBlas::gemv(1.0,zj,GPRModel,0.0,za);

        // derivative
        double loo = 0.0;
        #pragma omp parallel for reduction(+:loo)
        for(size_t i=0; i < GPRSize; i++){
            double zk = 0.0;
            for(size_t j=0; j < GPRSize; j++){
                zk += zj[i][j]*KS[j][i];
            }
            double top;
            top  = GPRModel[i]*za[i];
            top -= 0.5*(1.0 + GPRModel[i]*GPRModel[i]/KS[i][i])*zk;
            loo += top/KS[i][i];
        }

        // finalize derivative
        der[ind] += loo;
        ind++;
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CalcKderWRTSigmaF2(void)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];

            EneSurface->GetPoint(jbin,jpos);

            // calc Kblock
            if( UseNumDiff ){
                GetKernelDer2Num(ipos,jpos,kblock);
            } else {
                GetKernelDer2Ana(ipos,jpos,kblock);
            }

            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    // SigmaF2 cannot be zero
                    Kder[indi*NumOfCVs+ii][indj*NumOfCVs+jj] = kblock[ii][jj]/SigmaF2;
                }
            }
        }
    }

    if( UseZeroPoint ){
        CSimpleVector<double> kder;
        kder.CreateVector(NumOfCVs);

        #pragma omp parallel for firstprivate(ipos,kder)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t ibin = SampledMap[indi];

            EneSurface->GetPoint(ibin,ipos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumI(ipos,GPos,kder);
            } else {
                GetKernelDerAnaI(ipos,GPos,kder);
            }

            // distribute to matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                Kder[indi*NumOfCVs+ii][GPRSize-1] = kder[ii]/SigmaF2;
            }
        }

        #pragma omp parallel for firstprivate(jpos,kder)
        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];

            EneSurface->GetPoint(jbin,jpos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumJ(GPos,jpos,kder);
            } else {
                GetKernelDerAnaJ(GPos,jpos,kder);
            }

            // distribute to matrix
            for(size_t jj=0; jj < NumOfCVs; jj++){
                Kder[GPRSize-1][indj*NumOfCVs+jj] = kder[jj]/SigmaF2;
            }
        }

        Kder[GPRSize-1][GPRSize-1] = GetKernelValue(GPos,GPos)/SigmaF2;
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CalcKderWRTNCorr(void)
{
    Kder.SetZero();

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];

        for(size_t ii=0; ii < NumOfCVs; ii++){
            double er = item->GetValue(ibin,ii,E_PROXY_ERROR);
            Kder[indi*NumOfCVs+ii][indi*NumOfCVs+ii] = er*er;
        }
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CalcKderWRTWFac(size_t cv)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];

            EneSurface->GetPoint(jbin,jpos);

            GetKernelDer2AnaWFac(ipos,jpos,cv,kblock);

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    Kder[indi*NumOfCVs+ii][indj*NumOfCVs+jj] = kblock[ii][jj];
                }
            }
        }
    }

    if( UseZeroPoint ){
        RUNTIME_ERROR("not implemented");
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CalcKderWRTSigmaN2(size_t cv)
{
    Kder.SetZero();

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        for(size_t ii=0; ii < NumOfCVs; ii++){
            if( ii == cv ){
                Kder[indi*NumOfCVs+ii][indi*NumOfCVs+ii] = 1.0;
            }
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CIntegratorGPR::GetVar(CSimpleVector<double>& lpos)
{
    if( KSInverted != true ) {
        RUNTIME_ERROR("KS must be inverted!");
    }

    CSimpleVector<double>   kff;
    CSimpleVector<double>   ik;

    kff.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(lpos,kff);
    CSciBlas::gemv(1.0,KS,kff,0.0,ik);
    double cov = GetKernelValue(lpos,lpos) - CSciBlas::dot(kff,ik);
    return(cov);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& lrcov,double& rrvar)
{
    if( KSInverted != true ) {
        RUNTIME_ERROR("KS must be inverted!");
    }

    CSimpleVector<double>   kffr;
    CSimpleVector<double>   kffl;
    CSimpleVector<double>   ik;

    kffl.CreateVector(GPRSize);
    kffr.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(rpos,kffr);
    CreateKff(lpos,kffl);

    CSciBlas::gemv(1.0,KS,kffr,0.0,ik);

    lrcov = GetKernelValue(lpos,rpos) - CSciBlas::dot(kffl,ik);
    rrvar = GetKernelValue(rpos,rpos) - CSciBlas::dot(kffr,ik);
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }

    vout << "   Calculating FES error ..." << endl;
    CSmallTime st;
    st.GetActualTime();

    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    double  vargp = GetVar(gpos);
    int     nbatches = 0;

    #pragma omp parallel for firstprivate(jpos)
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        EneSurface->GetPoint(j,jpos);

        double varfc,covfg;
        GetCovVar(gpos,jpos,covfg,varfc);

        double error = varfc + vargp - 2.0*covfg;
        if( error > 0 ){
            error = sqrt(error);
        } else {
            error = 0.0;
        }
        EneSurface->SetError(j,error);

        #pragma omp atomic
        nbatches++;

#if defined(_OPENMP)
        int tnum = omp_get_thread_num();
#else
        int tnum = 0;
#endif
        if( tnum == 0){
            CSmallTime ct;
            ct.GetActualTime();
            if( (ct - st).GetSecondsFromBeginning() > 5*60 ){
                int comp = nbatches*100 / NumOfValues;
                vout << format("      completed %6d/%6d - %2d%%")%nbatches%NumOfValues%comp << endl;
                st = ct;
            }
        }
    }

    vout << "      RMSError  = " << setprecision(5) << EneSurface->GetRMSError() << endl;
    vout << "      MaxError  = " << setprecision(5) << EneSurface->GetMaxError() << endl;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::CalculateErrorsFromCov(CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }

    vout << "   Calculating FES error ..." << endl;

    // find global minimum
    size_t iglb = 0;

    if( GlobalMinSet ){
        iglb = NumOfValues; // the last item in Cov(i,i)
    } else {
        bool   first = true;
        double glb_min = 0.0;
        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            double value = EneSurface->GetEnergy(j);
            if( first || (glb_min > value) ){
                iglb = indj;
                first = false;
                glb_min = value;
            }
        }
    }

    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];

        double varfc = Cov[indj][indj];
        double covfg = Cov[indj][iglb];
        double vargp = Cov[iglb][iglb];

        double error = varfc + vargp - 2.0*covfg;
        if( error > 0 ){
            error = sqrt(error);
        } else {
            error = 0.0;
        }
        EneSurface->SetError(j,error);
    }

    vout << "      RMSError  = " << setprecision(5) << EneSurface->GetRMSError() << endl;
    vout << "      MaxError  = " << setprecision(5) << EneSurface->GetMaxError() << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CIntegratorGPR::CalculateCovs(CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }
    if( GPRSize == 0 ){
        RUNTIME_ERROR("GPRSize == 0");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("NumOfCVs == 0");
    }

// ------------------------------------
        vout << "   Calculating covariances ..." << endl;
    if( UseNumDiff ) {
        vout << "      Creating K+Sigma (numeric differentation) ..." << endl;
    } else {
        vout << "      Creating K+Sigma ..." << endl;
    }
        vout << "         Dim    = " << GPRSize << " x " << GPRSize << endl;

    CFortranMatrix KSInv;   // backup KS, which is inverted
    KSInv = KS;

    CreateKS();

    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

// ------------------------------------
    CSimpleVector<double> kff;
    kff.CreateVector(GPRSize);

    size_t nvals = NumOfValues;
    if( GlobalMinSet ){
        nvals++; // include explicitly set global minimum
    }

    CFortranMatrix  Kr;
    Kr.CreateMatrix(GPRSize,nvals);

        vout << "      Constructing kff ..." << endl;
        vout << "         Dim    = " << GPRSize << " x " << nvals << endl;

    #pragma omp parallel for firstprivate(ipos,kff)
    for(size_t indi=0; indi < NumOfValues; indi++){
        size_t i = ValueMap[indi];
        EneSurface->GetPoint(i,ipos);
        CreateKff(ipos,kff);
        for(size_t k=0; k < GPRSize; k++){
            Kr[k][indi] = kff[k];
        }
    }
    if( GlobalMinSet ){
        CreateKff(GPos,kff);
        for(size_t k=0; k < GPRSize; k++){
            Kr[k][nvals-1] = kff[k];
        }
    }

// ------------------------------------
    CFortranMatrix  Kl;
    Kl = Kr;

    RunBlasLapackPar();
    int result = 0;
    switch(Method){
        case(EGPRLA_LU):
            vout << format("      Solving (K+Sigma)^(-1)*kff by LU ...") << endl;
            result = CSciLapack::solvleLU(KS,Kr);
            if( result != 0 ) return;
            break;
        case(EGPRLA_LL):
            vout << format("      Solving (K+Sigma)^(-1)*kff by LL ...") << endl;
            result = CSciLapack::solvleLL(KS,Kr);
            if( result != 0 ) return;
            break;
            break;
    default:
        INVALID_ARGUMENT("unsupported method");
    }

        vout << "      Calculating Cov ..." << endl;
        vout << "         Dim    = " << nvals << " x " << nvals << endl;

// ------------------------------------
    Cov.CreateMatrix(nvals,nvals);

    #pragma omp parallel for firstprivate(ipos,jpos)
    for(size_t indi=0; indi < NumOfValues; indi++){
        size_t i = ValueMap[indi];
        EneSurface->GetPoint(i,ipos);
        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->GetPoint(j,jpos);
            Cov[indi][indj] = GetKernelValue(ipos,jpos);

        }
    }
    if( GlobalMinSet ){
        #pragma omp parallel for firstprivate(ipos)
        for(size_t indi=0; indi < NumOfValues; indi++){
            size_t i = ValueMap[indi];
            EneSurface->GetPoint(i,ipos);
            Cov[indi][nvals-1] = GetKernelValue(ipos,GPos);
            Cov[nvals-1][indi] = Cov[indi][nvals-1];
        }
        Cov[nvals-1][nvals-1] = GetKernelValue(GPos,GPos);
    }

    RunBlasLapackPar();
    CSciBlas::gemm(-1.0,'T',Kl,'N',Kr,1.0,Cov);

    // restore inverted KS
    KS = KSInv;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergySurfacePtr CIntegratorGPR::ReduceFES(const std::vector<bool>& keepcvs)
{
    if( keepcvs.size() != NumOfCVs ){
        RUNTIME_ERROR("keepcvs.size() != NumOfCVs");
    }

    bool only_variances = false;

    CEnergySurfacePtr p_rsurf = CEnergySurfacePtr(new CEnergySurface());
    double            temp    = EneSurface->GetTemperature();

    p_rsurf->Allocate(EneSurface,keepcvs);
    p_rsurf->Clear();

    CSimpleVector<int>    midx;
    midx.CreateVector(NumOfCVs);

    CSimpleVector<int>    ridx;
    ridx.CreateVector(p_rsurf->GetNumOfCVs());

    const double R = 1.98720425864083e-3;

    CSimpleVector<size_t>  IdxMap;
    IdxMap.CreateVector(NumOfValues);

// calculate weights
    for(size_t indi=0; indi < NumOfValues; indi++){
        size_t mbin = ValueMap[indi];
        double ene = EneSurface->GetEnergy(mbin);
        EneSurface->GetIPoint(mbin,midx);
        EneSurface->ReduceIPoint(keepcvs,midx,ridx);
        size_t rbin = p_rsurf->IPoint2Bin(ridx);
        IdxMap[indi] = rbin;
        double w = exp(-ene/(R*temp));
        p_rsurf->SetEnergy(rbin,p_rsurf->GetEnergy(rbin) + w);
        p_rsurf->SetNumOfSamples(rbin,1);
    }

// calculate errors
    for(size_t rbin = 0; rbin < (size_t)p_rsurf->GetNumOfBins(); rbin++){
        double err = 0.0;
        // err is now variance
        for(size_t indi=0; indi < NumOfValues; indi++){
            if( IdxMap[indi] != rbin ) continue;
            size_t mbini = ValueMap[indi];
            double enei = EneSurface->GetEnergy(mbini);
            double wi = exp(-enei/(R*temp));
            for(size_t indj=0; indj < NumOfValues; indj++){
                if( IdxMap[indj] != rbin ) continue;
                if( only_variances ){
                    if( indi != indj ) continue;
                }
                size_t mbinj = ValueMap[indj];
                double enej = EneSurface->GetEnergy(mbinj);
                double wj = exp(-enej/(R*temp));
                err = err + wi*wj*Cov[indi][indj];
            }
        }

        // p_rsurf->GetEnergy(rbin) contains sum of all weights
        err = err / (p_rsurf->GetEnergy(rbin)*p_rsurf->GetEnergy(rbin));

        // switch to error
        if( err > 0 ){
            err = sqrt(err);
        }
        p_rsurf->SetError(rbin,err);
    }

// transform back to FE
    for(size_t rbin = 0; rbin < (size_t)p_rsurf->GetNumOfBins(); rbin++){
        double w = p_rsurf->GetEnergy(rbin);
        p_rsurf->SetEnergy(rbin,-R*temp*log(w));
    }

// move global minimum
    double gmin = p_rsurf->GetGlobalMinimumValue();
    p_rsurf->ApplyOffset(-gmin);

    return(p_rsurf);
}

//------------------------------------------------------------------------------
