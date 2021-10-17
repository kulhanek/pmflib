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

#include <SmootherGPR.hpp>
#include <EnergySurface.hpp>
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

CSmootherGPR::CSmootherGPR(void)
{
    GPRSize             = 0;
//    NumOfUsedBins       = 0;
    NumOfValues         = 0;
    NumOfCVs                = 0;
    NumOfBins           = 0;

    SigmaF2             = 15.0;
    NCorr               = 1.0;

    IncludeError        = false;
    GlobalMinSet        = false;
    GPosSet             = false;
    GPosBin             = 0;

    Method              = EGPRLA_LU;
    Kernel              = EGPRK_ARDSE;

    NumOfThreads        = 1;

    UseInv              = false;
    NeedInv             = false;

    GlbMinValue         = 0.0;
}

//------------------------------------------------------------------------------

CSmootherGPR::~CSmootherGPR(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSmootherGPR::AddInputEnergyProxy(CEnergyProxyPtr p_prx)
{
    if( p_prx == NULL ) return;                 // no-proxy
    if( p_prx->GetAccu() == NULL ) return;      // no PMFAccu

    if( NumOfCVs != (size_t)p_prx->GetAccu()->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)p_prx->GetAccu()->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    EneProxyItems.push_back(p_prx);
}

//------------------------------------------------------------------------------

void CSmootherGPR::ClearInputEnergyProxies(void)
{
    EneProxyItems.clear();
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetOutputES(CEnergySurfacePtr p_surf)
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

void CSmootherGPR::SetSigmaF2(double sigf2)
{
    SigmaF2 = sigf2;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetNCorr(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetNCorr");
    }

    string          sspec(spec);
    vector<string>  sncorrs;

    split(sncorrs,sspec,is_any_of("x"),token_compress_on);

    if( sncorrs.size() > NumOfCVs ){
        CSmallString error;
        error << "too many ncorrs (" << sncorrs.size() << ") than required (" << NumOfCVs << ")";
        RUNTIME_ERROR(error);
    }

    if( sncorrs.size() > 1 ){
        CSmallString error;
        error << "too many ncorrs (" << sncorrs.size() << ") but only one allowed";
        RUNTIME_ERROR(error);
    }

    // parse values of ncorr
    double last_ncorr = 1.0;
    for(size_t i=0; i < sncorrs.size(); i++){
        stringstream str(sncorrs[i]);
        str >> last_ncorr;
        if( ! str ){
            CSmallString error;
            error << "unable to decode ncorr value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        NCorr = last_ncorr;
    }

    // pad the rest with the last value
    NCorr = last_ncorr;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetNCorr(double value)
{
    NCorr = value;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetWFac(const CSmallString& spec)
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

void CSmootherGPR::SetWFac(CSimpleVector<double>& wfac)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }
    if( wfac.GetLength() != NumOfCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    WFac = wfac;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetWFac(size_t cvind, double value)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
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

void CSmootherGPR::SetIncludeError(bool set)
{
    IncludeError = set;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetLAMethod(EGPRLAMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetLAMethod(const CSmallString& method)
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
                 "Supported methods are: svd (simple driver), svd2 (conquer and divide driver), lu, ll (Cholesky decomposition), default (=lu)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetKernel(const CSmallString& kernel)
{
    if( kernel == "ardse" ){
        Kernel = EGPRK_ARDSE;
    } else if( kernel == "ardmc52" ) {
        Kernel = EGPRK_ARDMC52;
    } else if( kernel == "ardmc32" ) {
        Kernel = EGPRK_ARDMC32;
    } else if( kernel == "ardmc12" ) {
        Kernel = EGPRK_ARDMC12;
    } else if( kernel == "default" ) {
        Kernel = EGPRK_ARDSE;
    } else {
        CSmallString error;
        error << "Specified kernel '" << kernel << "' is not supported. "
                 "Supported kernels are: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), ardmc32 (ARD Matern class 3/2), ardmc12 (ARD Matern class 1/2), default(=ardse)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetGlobalMin(const CSmallString& spec)
{
    GlobalMinSet = true;
    string sspec(spec);
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetGlobalMin");
    }

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

void CSmootherGPR::SetGlobalMin(const CSimpleVector<double>& pos)
{
    GlobalMinSet = true;
    GPos = pos;
    GPosSet = true;
}

//------------------------------------------------------------------------------

CSimpleVector<double> CSmootherGPR::GetGlobalMin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("")
    }
    return(GPos);
}

//------------------------------------------------------------------------------

int CSmootherGPR::GetGlobalMinBin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPosBin);
}

//------------------------------------------------------------------------------

double CSmootherGPR::GetGlobalMinValue(void) const
{
    return(GlbMinValue);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSmootherGPR::Interpolate(CVerboseStr& vout,bool nostat)
{
    PrintExecInfo(vout);

    if( EneProxyItems.size() == 0 ){
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

    // GPR setup
    CVLengths2.CreateVector(NumOfCVs);
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*EneSurface->GetCV(i)->GetRange()/EneSurface->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }

    // number of data points
    GPRSize = 0;
    for(size_t i=0; i < EneProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( EneProxyItems[i]->GetNumOfSamples(ibin) > 0 ) GPRSize++;
        }
    }

    // create sampled map
    SampledMap.resize(GPRSize);
    EneProxyMap.resize(GPRSize);
    size_t ind = 0;
    for(size_t i=0; i < EneProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( EneProxyItems[i]->GetNumOfSamples(ibin) <= 0 ) continue;
            SampledMap[ind] = ibin;
            EneProxyMap[ind] = i;
            ind++;
        }
    }

    // init GPR arrays
    GPRModel.CreateVector(GPRSize);
    Y.CreateVector(GPRSize);
    KS.CreateMatrix(GPRSize,GPRSize);

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
        vout << format("      SigmaF2   = %10.4f")%SigmaF2 << endl;
        vout << format("      NCorr     = %10.4f")%NCorr << endl;

    for(size_t k=0; k < NumOfCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    // finalize HES
    CalculateEnergy(vout);

    if( ! nostat ){
        // and log of marginal likelihood
        vout << "      logML     = " << setprecision(5) << GetLogML() << endl;
        if( NeedInv || UseInv ){
            // and log of pseudo-likelihood
            vout << "      logPL     = " << setprecision(5) << GetLogPL() << endl;
        }
    }

    if( IncludeError ){
        CalculateCovs(vout);
        CalculateErrorsFromCov(vout);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CSmootherGPR::PrintExecInfo(CVerboseStr& vout)
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

void CSmootherGPR::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CSmootherGPR::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//------------------------------------------------------------------------------

bool CSmootherGPR::WriteMFInfo(const CSmallString& name)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> ipos;
    CSimpleVector<double> mfi;
    CSimpleVector<double> mfp;

    ipos.CreateVector(NumOfCVs);
    mfi.CreateVector(GPRSize);
    mfp.CreateVector(GPRSize);

    // calculate
    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < GPRSize; indi++){
        CEnergyProxyPtr  item = EneProxyItems[EneProxyMap[indi]];
        size_t           ibin = SampledMap[indi];

        mfi[indi] = item->GetValue(ibin,E_PROXY_VALUE);
        mfp[indi] = GetValue(ipos);
    }

    ofstream ofs(name);
    if( ! ofs ){
        CSmallString error;
        error << "unable to open file '" << name << "' for derivatives";
        ES_ERROR(error);
        return(false);
    }

    // print
    size_t accu_prev = EneProxyMap[0];
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t ibin = SampledMap[indi];
        size_t accu = EneProxyMap[indi];

        if( accu != accu_prev ){
            ofs << endl;
            accu_prev = accu;
        }

        EneSurface->GetPoint(ibin,ipos);

        for(size_t c=0; c < NumOfCVs; c++){
            ofs << format("%20.16f ")%ipos[c];
        }

        ofs << format(" %20.16f %20.16f")%mfi[indi]%mfp[indi];

        ofs << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSmootherGPR::TrainGP(CVerboseStr& vout)
{
    vout << "   Creating K+Sigma and Y ..." << endl;
    vout << "      Kernel    = " << GetKernelName() << endl;
    vout << "      Dim       = " << GPRSize << " x " << GPRSize << endl;

    Mean = 0.0;
    for(size_t indi=0; indi < GPRSize; indi++){
        CEnergyProxyPtr item = EneProxyItems[EneProxyMap[indi]];
        size_t          ibin = SampledMap[indi];
        Mean += item->GetValue(ibin,E_PROXY_VALUE);
    }
    Mean /= (double)GPRSize;

// construct Y
    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        CEnergyProxyPtr item = EneProxyItems[EneProxyMap[indi]];
        size_t          ibin = SampledMap[indi];
        double mf   = item->GetValue(ibin,E_PROXY_VALUE);
        Y[indi]     = mf - Mean;
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
                // calculate weights
                vout << "   Calculating weights B ..." << endl;
                CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
            } else {
                GPRModel = Y;
                if( NeedInv ){
                    vout << "   Training GPR + K+Sigma inversion by LU ..." << endl;
                    result = CSciLapack::solvleLUInv(KS,GPRModel,logdetK);
                    if( result != 0 ) return(false);
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
                // calculate weights
                vout << "   Calculating weights B ..." << endl;
                CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
            } else {
                GPRModel = Y;
                if( NeedInv ){
                    vout << "   Training GPR + K+Sigma inversion by LL ..." << endl;
                    result = CSciLapack::solvleLLInv(KS,GPRModel,logdetK);
                    if( result != 0 ) return(false);
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
            }
            break;
    default:
        INVALID_ARGUMENT("unsupported method");
    }

    return(true);
}

//------------------------------------------------------------------------------

const CSmallString CSmootherGPR::GetKernelName(void)
{
    switch(Kernel){
    case(EGPRK_ARDSE):
        return("ARD squared exponential");
    case(EGPRK_ARDMC52):
        return("ARD Matern class 5/2");
    case(EGPRK_ARDMC32):
        return("ARD Matern class 3/2");
    case(EGPRK_ARDMC12):
        return("ARD Matern class 1/2");
    default:
        RUNTIME_ERROR("not implemented");
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::CreateKS(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    // main kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos)
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < GPRSize; indj++){
            size_t jbin = SampledMap[indj];

            EneSurface->GetPoint(jbin,jpos);
            KS[indi][indj] = GetKernelValue(ipos,jpos);
        }
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t          ibin = SampledMap[indi];
        CEnergyProxyPtr item = EneProxyItems[EneProxyMap[indi]];
        double er = item->GetValue(ibin,E_PROXY_ERROR);
        KS[indi][indi] += er*er*NCorr;
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& kff)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    // main kernel matrix
    for(size_t indj=0; indj < GPRSize; indj++){
        size_t  jbin = SampledMap[indj];

        EneSurface->GetPoint(jbin,jpos);
        kff[indj] = GetKernelValue(ip,jpos);
    }

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSmootherGPR::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating energy ..." << endl;

// create map for bins with calculated energy and error
    std::set<size_t>    vset;
    for(size_t i=0; i < EneProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( EneProxyItems[i]->GetNumOfSamples(ibin) <= 0 ) continue;
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

// basic HES update
    for(size_t i=0; i < EneProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            int nsamples = EneProxyItems[i]->GetNumOfSamples(ibin);
            int osamples = EneSurface->GetNumOfSamples(ibin);
            EneSurface->SetNumOfSamples(ibin,nsamples+osamples);
            EneSurface->SetEnergy(ibin,0.0);
            EneSurface->SetError(ibin,0.0);
        }
    }

// update HES
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        EneSurface->SetEnergy(j,values[indj]);
    }

// update HES
    if( GlobalMinSet ){
        // GPos.CreateVector(NumOfCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating EneSurface ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << EneSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << EneSurface->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << endl;

// find the closest bin
        CSimpleVector<double>   pos;
        pos.CreateVector(NumOfBins);
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

        double GlbMinValue = EneSurface->GetEnergy(GPosBin);

        GlbMinValue = GetValue(GPos);
        vout << " (" << setprecision(5) << GlbMinValue << ")" << endl;
        vout << "      Offset    = " << setprecision(5) << GlbMinValue << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->SetEnergy(j,EneSurface->GetEnergy(j)-GlbMinValue);
        }
    } else {
        // search for global minimum
        GPos.CreateVector(NumOfCVs);
        bool   first = true;
        GlbMinValue = 0.0;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            int samples = EneSurface->GetNumOfSamples(j);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = values[indj];
            if( first || (GlbMinValue > value) ){
                GlbMinValue = value;
                first = false;
                EneSurface->GetPoint(j,GPos);
            }
        }

   //   vout << "   Calculating EneSurface ..." << endl;
        vout << "      Global minimum found at: ";
        vout << setprecision(5) << EneSurface->GetCV(0)->GetRealValue(GPos[0]);
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << EneSurface->GetCV(i)->GetRealValue(GPos[i]);
        }
        vout << " (" << setprecision(5) << GlbMinValue << ")" << endl;
        vout << "      Offset    = " << setprecision(5) << GlbMinValue << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->SetEnergy(j,values[indj]-GlbMinValue);
        }
    }

    GPosSet = true;

    vout << "      SigmaF2   = " << setprecision(5) << EneSurface->GetSigmaF2() << endl;
    vout << "      SigmaF    = " << setprecision(5) << EneSurface->GetSigmaF() << endl;
}

//------------------------------------------------------------------------------

double CSmootherGPR::GetValue(const CSimpleVector<double>& position)
{
    CSimpleVector<double>   kff;
    kff.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(position,kff);
    double energy = CSciBlas::dot(kff,GPRModel);
    return(energy + Mean);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSmootherGPR::CalculateErrorsFromCov(CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }

    vout << "   Calculating enthalpy error ..." << endl;

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

        // cout << varfc << " " << covfg << " " << vargp << endl;

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

void CSmootherGPR::CalculateCovs(CVerboseStr& vout)
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
    vout << "      Creating K+Sigma ..." << endl;
    vout << "         Dim    = " << GPRSize << " x " << GPRSize << endl;
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
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSmootherGPR::SetUseInv(bool set)
{
    UseInv = set;
}

//------------------------------------------------------------------------------

void CSmootherGPR::PrepForHyprmsGrd(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CSmootherGPR::SetCalcLogPL(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

double CSmootherGPR::GetLogML(void)
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

void CSmootherGPR::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
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
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( prm == 1 ){
            CalcKderWRTNCorr();
        } else {
            size_t cv = prm - 2;
            CalcKderWRTWFac(cv);
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

double CSmootherGPR::GetLogPL(void)
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

void CSmootherGPR::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
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
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( prm == 1 ){
            CalcKderWRTNCorr();
        } else {
            size_t cv = prm - 2;
            CalcKderWRTWFac(cv);
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

void CSmootherGPR::CalcKderWRTSigmaF2(void)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    #pragma omp parallel for firstprivate(ipos,jpos)
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t i = SampledMap[indi];
        EneSurface->GetPoint(i,ipos);

        for(size_t indj=0; indj < GPRSize; indj++){
            size_t j = SampledMap[indj];
            EneSurface->GetPoint(j,jpos);

            Kder[indi][indj] = GetKernelValue(ipos,jpos) / SigmaF2;
        }
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::CalcKderWRTNCorr(void)
{
    Kder.SetZero();

    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t          ibin = SampledMap[indi];
        CEnergyProxyPtr item = EneProxyItems[EneProxyMap[indi]];
        double er = item->GetValue(ibin,E_PROXY_ERROR);
        Kder[indi][indi] += er*er;
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::CalcKderWRTWFac(size_t cv)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    #pragma omp parallel for firstprivate(ipos,jpos)
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t i = SampledMap[indi];
        EneSurface->GetPoint(i,ipos);

        for(size_t indj=0; indj < GPRSize; indj++){
            size_t j = SampledMap[indj];
            EneSurface->GetPoint(j,jpos);

            if( indi != indj ){
                Kder[indi][indj] = GetKernelValueWFacDer(ipos,jpos,cv);
            } else {
                // FIXME - check validity - this avoids division by zero in GetKernelValueWFacDer
                Kder[indi][indj] = 0.0;
            }
        }
    }
}

//------------------------------------------------------------------------------

double CSmootherGPR::GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
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
        case(EGPRK_ARDSE): {
            return(SigmaF2*exp(-0.5*scdist2));
        }
        break;
    // -----------------------
        case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            return(SigmaF2*(1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2)*exp(-sqrt(5.0)*scdist));
        }
        break;
    // -----------------------
        case(EGPRK_ARDMC32):{
            double scdist = sqrt(scdist2);
            return(SigmaF2*(1.0+sqrt(3.0)*scdist)*exp(-sqrt(3.0)*scdist));
        }
        break;
    // -----------------------
        case(EGPRK_ARDMC12):{
            double scdist = sqrt(scdist2);
            return(SigmaF2*exp(-scdist));
        }
        break;
    // -----------------------
        default:
            RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

double CSmootherGPR::GetKernelValueWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = EneSurface->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    switch(Kernel){
        case(EGPRK_ARDSE): {
                double val = SigmaF2*exp(-0.5*scdist2);
                double du = EneSurface->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
                double dd = CVLengths2[cv];
                double wf = WFac[cv];
                double der = val * du*du / (dd * wf);
                return(der);
            }
            break;
    // -----------------------
        case(EGPRK_ARDMC52): {
            // possible division by zero is solved in CalcKderWRTWFac
            double scdist = sqrt(scdist2);
            double vex = SigmaF2*exp(-sqrt(5.0)*scdist);
            double val = vex*(1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2);
            double du = EneSurface->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double dexp =   val * sqrt(5.0) * du*du / (dd * wf) / scdist;       // exponential part
            double dpol = - vex * ( sqrt(5.0) * du*du / (dd * wf) / scdist + 2.0 * (5.0/3.0)* du*du / (dd * wf) );   // polynomial part
            return(dexp+dpol);
            }
            break;
    // -----------------------
        case(EGPRK_ARDMC32): {
            // possible division by zero is solved in CalcKderWRTWFac
            double scdist = sqrt(scdist2);
            double vex = SigmaF2*exp(-sqrt(3.0)*scdist);
            double val = vex*(1.0+sqrt(3.0)*scdist);
            double du = EneSurface->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double dexp =   val * sqrt(3.0) * du*du / (dd * wf) / scdist;   // exponential part
            double dpol = - vex * sqrt(3.0) * du*du / (dd * wf) / scdist;   // polynomial part
            return(dexp+dpol);
            }
            break;
    // -----------------------
        case(EGPRK_ARDMC12): {
            // possible division by zero is solved in CalcKderWRTWFac
            double scdist = sqrt(scdist2);
            double val = SigmaF2*exp(-scdist);
            double du = EneSurface->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double der = val * du*du / (dd * wf) / scdist;
            return(der);
            }
            break;
    // -----------------------
        default:
            RUNTIME_ERROR("not implemented");
            break;
    }
}

//------------------------------------------------------------------------------
