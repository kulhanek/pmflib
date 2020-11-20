// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <ABFIntegratorGPR.hpp>
#include <ABFAccumulator.hpp>
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

CABFIntegratorGPR::CABFIntegratorGPR(void)
{
    Accumulator         = NULL;
    FES                 = NULL;

    GPRSize             = 0;
    NumOfUsedBins       = 0;
    NumOfValues         = 0;
    NCVs                = 0;
    NumOfBins           = 0;

    SplitNCorr          = false;
    SigmaF2             = 15.0;

    IncludeError        = false;
    IncludeGluedBins    = false;
    NoEnergy            = false;
    GlobalMinSet        = false;
    GPosSet             = false;

    UseNumDiff          = false;
    Method              = EGPRLA_LU;
    Kernel              = EGPRK_ARDSE;

    NumOfThreads        = 1;

    UseInv              = false;
    NeedInv             = false;
    UseZeroPoint        = false;
    FastErrors          = true;
}

//------------------------------------------------------------------------------

CABFIntegratorGPR::~CABFIntegratorGPR(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorGPR::SetInputABFAccumulator(CABFAccumulator* p_accu)
{
    Accumulator = p_accu;

    if( Accumulator ){
        NCVs = Accumulator->GetNumberOfCoords();
        NumOfBins = Accumulator->GetNumberOfBins();
    } else {
        NCVs = 0;
        NumOfBins = 0;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetOutputFESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetSigmaF2(double sigf2)
{
    SigmaF2 = sigf2;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetNCorr(const CSmallString& spec)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetNCorr");
    }

    string          sspec(spec);
    vector<string>  sncorrs;

    split(sncorrs,sspec,is_any_of("x"),token_compress_on);

    if( sncorrs.size() > NCVs ){
        CSmallString error;
        error << "too many ncorrs (" << sncorrs.size() << ") than required (" << NCVs << ")";
        RUNTIME_ERROR(error);
    }

    if( (SplitNCorr == false) && (sncorrs.size() > 1) ){
        CSmallString error;
        error << "too many ncorrs (" << sncorrs.size() << ") but only one allowed without --splitncorr";
        RUNTIME_ERROR(error);
    }

    NCorr.CreateVector(NCVs);

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
        NCorr[i] = last_ncorr;
    }

    // pad the rest with the last value
    for(size_t i=sncorrs.size(); i < NCVs; i++){
        NCorr[i] = last_ncorr;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetSplitNCorr(bool set)
{
    SplitNCorr = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetNCorr(CSimpleVector<double>& ncorr)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetNCorr");
    }
    if( ncorr.GetLength() != NCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    NCorr = ncorr;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetNCorr(size_t cvind, double value)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetNCorr");
    }
    if( cvind >= NCVs ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is NCorr initialized?
    if( NCorr.GetLength() == 0 ){
        NCorr.CreateVector(NCVs);
    }
    NCorr[cvind] = value;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetWFac(const CSmallString& spec)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }

    string          sspec(spec);
    vector<string>  swfacs;

    split(swfacs,sspec,is_any_of("x"),token_compress_on);

    if( swfacs.size() > NCVs ){
        CSmallString error;
        error << "too many wfacs (" << swfacs.size() << ") than required (" << NCVs << ")";
        RUNTIME_ERROR(error);
    }

    WFac.CreateVector(NCVs);

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
    for(size_t i=swfacs.size(); i < NCVs; i++){
        WFac[i] = last_wfac;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetWFac(CSimpleVector<double>& wfac)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }
    if( wfac.GetLength() != NCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    WFac = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetWFac(size_t cvind, double value)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }
    if( cvind >= NCVs ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is wfac initialized?
    if( WFac.GetLength() == 0 ){
        WFac.CreateVector(NCVs);
    }
    WFac[cvind] = value;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetIncludeError(bool set)
{
    IncludeError = set;
    if( FastErrors == false ){
        NeedInv |= set;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetUseNumDiff(bool set)
{
    UseNumDiff = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetNoEnergy(bool set)
{
    NoEnergy = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetLAMethod(EGPRLAMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetLAMethod(const CSmallString& method)
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

void CABFIntegratorGPR::SetKernel(const CSmallString& kernel)
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

void CABFIntegratorGPR::IncludeGluedAreas(bool set)
{
    IncludeGluedBins = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetGlobalMin(const CSmallString& spec)
{
    GlobalMinSet = true;
    string sspec(spec);
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetGlobalMin");
    }

    // remove "x" from the string
    replace (sspec.begin(), sspec.end(), 'x' , ' ');

    // parse values of CVs
    GPos.CreateVector(NCVs);
    stringstream str(sspec);
    for(size_t i=0; i < NCVs; i++){
        str >> GPos[i];
        if( ! str ){
            CSmallString error;
            error << "unable to decode CV value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

    GPosSet = true;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetGlobalMin(const CSimpleVector<double>& pos)
{
    GlobalMinSet = true;
    GPos = pos;
    GPosSet = true;
}

//------------------------------------------------------------------------------

CSimpleVector<double> CABFIntegratorGPR::GetGlobalMin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("")
    }
    return(GPos);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetUseInv(bool set)
{
    UseInv = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::PrepForHyprmsGrd(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetCalcLogPL(bool set)
{
   NeedInv |= set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetUseZeroPoint(bool set)
{
    if( set == true ){
        if( GlobalMinSet == false ){
            RUNTIME_ERROR("position of global minimum was not set");
        }
    }
    UseZeroPoint = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetFastError(bool set)
{
    FastErrors = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::Integrate(CVerboseStr& vout,bool nostat)
{
    PrintExecInfo(vout);

    if( Accumulator == NULL ) {
        RUNTIME_ERROR("ABF accumulator is not set");
    }
    if( FES == NULL ) {
        RUNTIME_ERROR("FES is not set");
    }
    if( Accumulator->GetNumberOfCoords() == 0 ) {
        RUNTIME_ERROR("number of coordinates is zero");
    }
    if( Accumulator->GetNumberOfCoords() != FES->GetNumberOfCoords() ){
        RUNTIME_ERROR("inconsistent ABF and FES - CVs");
    }
    if( Accumulator->GetNumberOfBins() != FES->GetNumberOfPoints() ){
        RUNTIME_ERROR("inconsistent ABF and FES - points");
    }
    if( WFac.GetLength() == 0 ){
        RUNTIME_ERROR("wfac is not set");
    }
    if( NCorr.GetLength() == 0 ){
        RUNTIME_ERROR("ncorr is not set");
    }

    // GPR setup
    CVLengths2.CreateVector(NCVs);
    for(size_t i=0; i < NCVs; i++){
        double l = WFac[i]*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        CVLengths2[i] = l*l;
    }

    // number of data points
    NumOfUsedBins = 0;
    for(size_t i=0; i < NumOfBins; i++){
        if( Accumulator->GetNumberOfABFSamples(i) > 0 ) NumOfUsedBins++;
    }
    GPRSize = NumOfUsedBins*NCVs;
    if( UseZeroPoint ) GPRSize++;

    // create sampled map
    SampledMap.CreateVector(NumOfUsedBins);
    size_t ind = 0;
    for(size_t i=0; i < NumOfBins; i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;
        SampledMap[ind] = i;
        ind++;
    }

    // init GPR arrays
    GPRModel.CreateVector(GPRSize);
    Y.CreateVector(GPRSize);
    KS.CreateMatrix(GPRSize,GPRSize);

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
        vout << format("      SigmaF2   = %10.4f")%SigmaF2 << endl;
    if( SplitNCorr ){
        for(size_t k=0; k < NCVs; k++ ){
            vout << format("      NCorr#%-2d  = %10.4f")%(k+1)%NCorr[k] << endl;
        }
    } else {
        vout << format("      NCorr     = %10.4f")%NCorr[0] << endl;
    }

    for(size_t k=0; k < NCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    if( ! nostat ){
        // and finaly some statistics
        for(size_t k=0; k < NCVs; k++ ){
        vout << "      RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
        }

        // and lof of marginal likelihood
        vout << "      logML     = " << setprecision(5) << GetLogML() << endl;
        if( NeedInv || UseInv ){
            // and log of pseudo-likelihood
            vout << "      logPL     = " << setprecision(5) << GetLogPL() << endl;
        }
    }

    // finalize FES if requested
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

void CABFIntegratorGPR::PrintExecInfo(CVerboseStr& vout)
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

void CABFIntegratorGPR::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::TrainGP(CVerboseStr& vout)
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
        size_t i = SampledMap[indi];
        for(size_t ii=0; ii < NCVs; ii++){
            double mf = Accumulator->GetValue(ii,i,EABF_DG_VALUE);
            Y[indi*NCVs+ii] = mf;
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

const CSmallString CABFIntegratorGPR::GetKernelName(void)
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

void CABFIntegratorGPR::CreateKS(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NCVs,NCVs);

    // main kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            // calc Kblock
            if( UseNumDiff ){
                GetKernelDer2Num(ipos,jpos,kblock);
            } else {
                GetKernelDer2Ana(ipos,jpos,kblock);
            }

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    KS[indi*NCVs+ii][indj*NCVs+jj] = kblock[ii][jj];
                }
            }
        }
    }
    if( UseZeroPoint ){
        CSimpleVector<double> kder;
        kder.CreateVector(NCVs);

        #pragma omp parallel for firstprivate(ipos,kder)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            Accumulator->GetPoint(i,ipos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumI(ipos,GPos,kder);
            } else {
                GetKernelDerAnaI(ipos,GPos,kder);
            }

            // distribute to matrix
            for(size_t ii=0; ii < NCVs; ii++){
                KS[indi*NCVs+ii][GPRSize-1] = kder[ii];
            }
        }

        #pragma omp parallel for firstprivate(jpos,kder)
        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumJ(GPos,jpos,kder);
            } else {
                GetKernelDerAnaJ(GPos,jpos,kder);
            }

            // distribute to matrix
            for(size_t jj=0; jj < NCVs; jj++){
                KS[GPRSize-1][indj*NCVs+jj] = kder[jj];
            }
        }

        KS[GPRSize-1][GPRSize-1] = GetKernelValue(GPos,GPos);
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        for(size_t ii=0; ii < NCVs; ii++){
            double er = Accumulator->GetValue(ii,i,EABF_DG_ERROR);
            if( SplitNCorr ){
                KS[indi*NCVs+ii][indi*NCVs+ii] += er*er*NCorr[ii];
            } else {
                KS[indi*NCVs+ii][indi*NCVs+ii] += er*er*NCorr[0];
            }
        }
    }
    // zero point has no uncertainity
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& kff)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    CSimpleVector<double> kder;
    kder.CreateVector(NCVs);

    // main kernel matrix
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t j = SampledMap[indj];

        Accumulator->GetPoint(j,jpos);

        // calc Kder
        if( UseNumDiff ){
            GetKernelDerNumJ(ip,jpos,kder);
        } else {
            GetKernelDerAnaJ(ip,jpos,kder);
        }

        // distribute to vector
        for(size_t jj=0; jj < NCVs; jj++){
            kff[indj*NCVs+jj] = kder[jj];
        }
    }

    if( UseZeroPoint ){
        kff[GPRSize-1] = GetKernelValue(ip,GPos);
    }

}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& kff2)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NCVs,NCVs);

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t j = SampledMap[indj];

        Accumulator->GetPoint(j,jpos);

        // calc Kder
        if( UseNumDiff ){
            GetKernelDer2Num(ip,jpos,kblock);
        } else {
            GetKernelDer2Ana(ip,jpos,kblock);
        }

        // distribute to vector
        for(size_t jj=0; jj < NCVs; jj++){
            kff2[indj*NCVs+jj] = kblock[icoord][jj];
        }
    }

    if( UseZeroPoint ){
        CSimpleVector<double> kder;
        kder.CreateVector(NCVs);

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

void CABFIntegratorGPR::GetKernelDerAnaI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = SigmaF2*exp(-0.5*scdist2);
            for(size_t ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*SigmaF2*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
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

void CABFIntegratorGPR::GetKernelDerNumI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NCVs);

    double  dh = 1e-3;
    double  v1,v2;

    // off diagonal elements
    for(size_t ii=0; ii < NCVs; ii++) {

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

void CABFIntegratorGPR::GetKernelDerAnaJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = SigmaF2*exp(-0.5*scdist2);
            for(size_t ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*SigmaF2*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
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

void CABFIntegratorGPR::GetKernelDerNumJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    CSimpleVector<double> tjp;
    tjp.CreateVector(NCVs);

    double  dh = 1e-3;
    double  v1,v2;

    // off diagonal elements
    for(size_t ii=0; ii < NCVs; ii++) {

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

void CABFIntegratorGPR::GetKernelDer2Ana(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    kblock.SetZero();

    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double pre = SigmaF2*exp(-0.5*scdist2);
            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ip[jj],jp[jj]);
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
            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ip[jj],jp[jj]);
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

void CABFIntegratorGPR::GetKernelDer2Num(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NCVs);
    CSimpleVector<double> tjp;
    tjp.CreateVector(NCVs);

    double  dh = 1e-3;
    double  v1,v2,v3,v4;

    // off diagonal elements
    for(size_t ii=0; ii < NCVs; ii++) {
        for(size_t jj=0; jj < NCVs; jj++) {

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

void CABFIntegratorGPR::GetKernelDer2AnaWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    kblock.SetZero();

    double wf = WFac[cv];
    double wd3 = 1.0/(CVLengths2[cv]*wf);
    double wd5 = wd3/CVLengths2[cv];
    double dc = Accumulator->GetCoordinate(cv)->GetDifference(ip[cv],jp[cv]);

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double arg = SigmaF2*exp(-0.5*scdist2);
            double argd = arg*dc*dc*wd3;
            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ip[jj],jp[jj]);
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
            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ip[jj],jp[jj]);
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
void CABFIntegratorGPR::GetKernelDer2NumWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    CFortranMatrix kblock1,kblock2;

    kblock1.CreateMatrix(NCVs,NCVs);
    kblock2.CreateMatrix(NCVs,NCVs);

    double  dh = 1e-3;

    for(size_t i=0; i < NCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        } else {
            l = WFac[i]*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDer2Ana(ip,jp,kblock1);

    for(size_t i=0; i < NCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        } else {
            l = WFac[i]*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDer2Ana(ip,jp,kblock2);

    for(size_t ii=0; ii < NCVs; ii++){
        for(size_t jj=0; jj < NCVs; jj++){
            kblock[ii][jj] = (kblock2[ii][jj]-kblock1[ii][jj])/(2.0*dh);
        }
    }

    // restore original CVLengths2
    for(size_t i=0; i < NCVs; i++){
        double l = WFac[i]*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        CVLengths2[i] = l*l;
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
{
    // get kernel value
    switch(Kernel){
        case(EGPRK_ARDSE): {
            // calculate scaled distance
            double scdist2 = 0.0;
            for(size_t ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                scdist2 += du*du/dd;
            }
            return(SigmaF2*exp(-0.5*scdist2));
        }
        break;
        case(EGPRK_ARDMC52):{
            // calculate scaled distance
            double scdist2 = 0.0;
            for(size_t ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
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

void CABFIntegratorGPR::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating FES ..." << endl;

// create map for bins with calculated energy and error
    NumOfValues = 0;
    for(size_t i=0; i < NumOfBins; i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        if( IncludeGluedBins ){
            if( samples == 0 ) continue;
        } else {
            if( samples <= 0 ) continue;
        }
        NumOfValues++;
    }
    ValueMap.CreateVector(NumOfValues);

    size_t indi = 0;
    for(size_t i=0; i < NumOfBins; i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        if( IncludeGluedBins ){
            if( samples == 0 ) continue;
        } else {
            if( samples <= 0 ) continue;
        }
        ValueMap[indi]=i;
        indi++;
    }

// calculate energies
    CSimpleVector<double> jpos;
    CSimpleVector<double> values;

    jpos.CreateVector(NCVs);
    values.CreateVector(NumOfValues);

    #pragma omp parallel for firstprivate(jpos)
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        Accumulator->GetPoint(j,jpos);
        values[indj] = GetValue(jpos);
    }

// basic FES update
    for(size_t i=0; i < NumOfBins; i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        FES->SetNumOfSamples(i,samples);
        FES->SetEnergy(i,0.0);
        FES->SetError(i,0.0);
    }

// update FES
    if( GlobalMinSet ){
        // GPos.CreateVector(NCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << GPos[0];
        for(int i=1; i < Accumulator->GetNumberOfCoords(); i++){
            vout << "x" << GPos[i];
        }
        double glb_min = GetValue(GPos);
        vout << " (" << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            FES->SetEnergy(j,values[indj]-glb_min);
        }
    } else {
        // search for global minimum
        GPos.CreateVector(NCVs);
        bool   first = true;
        double glb_min = 0.0;
        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            int samples = Accumulator->GetNumberOfABFSamples(j);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = values[indj];
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                Accumulator->GetPoint(j,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "      Global minimum found at: ";
        vout << GPos[0];
        for(size_t i=1; i < NCVs; i++){
            vout << "x" << GPos[i];
        }
        vout << " (" << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            FES->SetEnergy(j,values[indj]-glb_min);
        }
    }

    GPosSet = true;

    vout << "      SigmaF2   = " << setprecision(5) << FES->GetSigmaF2() << endl;
    if( IncludeGluedBins ){
        vout << "      SigmaF2 (including glued bins) = " << setprecision(5) << FES->GetSigmaF2(true) << endl;
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetValue(const CSimpleVector<double>& position)
{
    CSimpleVector<double>   kff;
    kff.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(position,kff);
    double energy = CSciBlas::dot(kff,GPRModel);
    return(energy);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetMeanForce(const CSimpleVector<double>& position,size_t icoord)
{
    CSimpleVector<double>   kff2;
    kff2.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff2(position,icoord,kff2);
    double mf = CSciBlas::dot(kff2,GPRModel);
    return(mf);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetRMSR(size_t cv)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> ipos;
    ipos.CreateVector(NCVs);

    double rmsr = 0.0;

    #pragma omp parallel for firstprivate(ipos) reduction(+:rmsr)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        double mfi = Accumulator->GetValue(cv,i,EABF_DG_VALUE);
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

bool CABFIntegratorGPR::WriteMFInfo(const CSmallString& name)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> ipos;
    CSimpleVector<double> mfi;
    CSimpleVector<double> mfp;

    ipos.CreateVector(NCVs);
    mfi.CreateVector(GPRSize);
    mfp.CreateVector(GPRSize);

    // calculate
    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        Accumulator->GetPoint(i,ipos);
        for(size_t k=0; k < NCVs; k++){
            mfi[indi*NCVs+k] = Accumulator->GetValue(k,i,EABF_DG_VALUE);
            mfp[indi*NCVs+k] = GetMeanForce(ipos,k);
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
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(size_t c=0; c < NCVs; c++){
            ofs << format("%20.16f ")%ipos[c];
        }

        for(size_t k=0; k < NCVs; k++){
            ofs << format(" %20.16f %20.16f")%mfi[indi*NCVs+k]%mfp[indi*NCVs+k];
        }

        ofs << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::FilterByMFZScore(double zscore,CVerboseStr& vout)
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
    sig2.CreateVector(NCVs);
    maxi.CreateVector(NCVs);
    maxzscore.CreateVector(NCVs);
    mferror2.CreateVector(GPRSize);

    flags.Set(1);

    vout << high;
    vout << "      Precalculating MF errors ..." << endl;

    CSimpleVector<double> ipos;
    ipos.CreateVector(NCVs);

    // precalculate values
    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(size_t k=0; k < NCVs; k++){
            double diff2 = Accumulator->GetValue(k,i,EABF_DG_VALUE) - GetMeanForce(ipos,k);
            diff2 *= diff2;
            mferror2[indi*NCVs+k] = diff2;
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
                for(size_t k=0; k < NCVs; k++){
                    sig2[k] += mferror2[indi*NCVs+k];
                }
                count++;
            }
        }

        if( count == 0 ) break;
        for(size_t k=0; k < NCVs; k++){
            sig2[k] /= count;
        }

        // filter
        bool first = true;
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            if( flags[i] != 0 ){

                for(size_t k=0; k < NCVs; k++){
                    double diff2 = mferror2[indi*NCVs+k];
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

        for(size_t k=0; k < NCVs; k++){
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
            Accumulator->SetNumberOfABFSamples(i,0);
            outliers++;
        }
    }

    vout << high;
    vout << "      Number of outliers = " << outliers << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CABFIntegratorGPR::GetLogML(void)
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

double CABFIntegratorGPR::GetLogPL(void)
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

void CABFIntegratorGPR::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("GetLogMLDerivatives requires K+Sigma inverted matrix");
    }

    if( Accumulator == NULL ) {
        RUNTIME_ERROR("ABF accumulator is not set");
    }
    if( FES == NULL ) {
        RUNTIME_ERROR("FES is not set");
    }
    if( NCVs <= 0 ){
        RUNTIME_ERROR("NCVs <= NULL");
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
    size_t soffset = 0;
    size_t noffset;
    size_t woffset;

    if( SplitNCorr ){
        noffset = soffset+1;
        woffset = noffset+NCVs;
    } else {
        noffset = soffset+1;
        woffset = noffset+1;
    }

    for(size_t prm=0; prm < flags.size(); prm++){
        // shall we calc der?
        if( flags[prm] == false ) {
            continue;
        }

        // calc Kder
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( (prm >= noffset) && (prm < woffset) ){
            size_t cv = prm - noffset;
            CalcKderWRTNCorr(cv);
        } else {
            size_t cv = prm - woffset;
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

void CABFIntegratorGPR::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("GetLogPLDerivatives requires K+Sigma inverted matrix");
    }

    if( Accumulator == NULL ) {
        RUNTIME_ERROR("ABF accumulator is not set");
    }
    if( FES == NULL ) {
        RUNTIME_ERROR("FES is not set");
    }
    if( NCVs <= 0 ){
        RUNTIME_ERROR("NCVs <= NULL");
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
    size_t soffset = 0;
    size_t noffset;
    size_t woffset;

    if( SplitNCorr ){
        noffset = soffset+1;
        woffset = noffset+NCVs;
    } else {
        noffset = soffset+1;
        woffset = noffset+1;
    }

    for(size_t prm=0; prm < flags.size(); prm++){
        // shall we calc der?
        if( flags[prm] == false ) {
            continue;
        }

        // calc Kder
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( (prm >= noffset) && (prm < woffset) ){
            size_t cv = prm - noffset;
            CalcKderWRTNCorr(cv);
        } else {
            size_t cv = prm - woffset;
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

void CABFIntegratorGPR::CalcKderWRTSigmaF2(void)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NCVs,NCVs);

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            // calc Kblock
            if( UseNumDiff ){
                GetKernelDer2Num(ipos,jpos,kblock);
            } else {
                GetKernelDer2Ana(ipos,jpos,kblock);
            }

            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    // SigmaF2 cannot be zero
                    Kder[indi*NCVs+ii][indj*NCVs+jj] = kblock[ii][jj]/SigmaF2;
                }
            }
        }
    }

    if( UseZeroPoint ){
        CSimpleVector<double> kder;
        kder.CreateVector(NCVs);

        #pragma omp parallel for firstprivate(ipos,kder)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            Accumulator->GetPoint(i,ipos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumI(ipos,GPos,kder);
            } else {
                GetKernelDerAnaI(ipos,GPos,kder);
            }

            // distribute to matrix
            for(size_t ii=0; ii < NCVs; ii++){
                Kder[indi*NCVs+ii][GPRSize-1] = kder[ii]/SigmaF2;
            }
        }

        #pragma omp parallel for firstprivate(jpos,kder)
        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            // calc Kder
            if( UseNumDiff ){
                GetKernelDerNumJ(GPos,jpos,kder);
            } else {
                GetKernelDerAnaJ(GPos,jpos,kder);
            }

            // distribute to matrix
            for(size_t jj=0; jj < NCVs; jj++){
                Kder[GPRSize-1][indj*NCVs+jj] = kder[jj]/SigmaF2;
            }
        }

        Kder[GPRSize-1][GPRSize-1] = GetKernelValue(GPos,GPos)/SigmaF2;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CalcKderWRTNCorr(size_t cv)
{
    Kder.SetZero();

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        for(size_t ii=0; ii < NCVs; ii++){
            double er = Accumulator->GetValue(ii,i,EABF_DG_ERROR);
            if( SplitNCorr ){
                if( cv == ii ) Kder[indi*NCVs+ii][indi*NCVs+ii] += er*er;
            } else {
                Kder[indi*NCVs+ii][indi*NCVs+ii] += er*er;
            }
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CalcKderWRTWFac(size_t cv)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NCVs,NCVs);

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            GetKernelDer2AnaWFac(ipos,jpos,cv,kblock);

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NCVs; ii++){
                for(size_t jj=0; jj < NCVs; jj++){
                    // SigmaF2 cannot be zero
                    Kder[indi*NCVs+ii][indj*NCVs+jj] = kblock[ii][jj];
                }
            }
        }
    }

    if( UseZeroPoint ){
        RUNTIME_ERROR("not implemented");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CABFIntegratorGPR::GetVar(CSimpleVector<double>& lpos)
{
    CSimpleVector<double>   kff;
    CSimpleVector<double>   ik;

    kff.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(lpos,kff);
    CSciBlas::gemv(1.0,KS,kff,0.0,ik);
    double cov = GetKernelValue(lpos,lpos)  - CSciBlas::dot(kff,ik);
    return(cov);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& lrcov,double& rrvar)
{
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

void CABFIntegratorGPR::CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }

    vout << "   Calculating FES error ..." << endl;
    CSmallTime st;
    st.GetActualTime();

    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    double  vargp = GetVar(gpos);
    int     nbatches = 0;

    #pragma omp parallel for firstprivate(jpos)
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        Accumulator->GetPoint(j,jpos);

        double varfc,covfg;
        GetCovVar(gpos,jpos,covfg,varfc);

        // cout << varfc << " " << covfg << " " << vargp << endl;

        double error = varfc + vargp - 2.0*covfg;
        if( error > 0 ){
            error = sqrt(error);
        } else {
            error = 0.0;
        }
        FES->SetError(j,error);

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
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CalculateErrorsFromCov(CVerboseStr& vout)
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
            double value = FES->GetEnergy(j);
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
        FES->SetError(j,error);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorGPR::CalculateCovs(CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }
    if( GPRSize == 0 ){
        RUNTIME_ERROR("GPRSize == 0");
    }
    if( NCVs == 0 ){
        RUNTIME_ERROR("NCVs == 0");
    }

// ------------------------------------
        vout << "   Calculating covariances ..." << endl;
    if( UseNumDiff ) {
        vout << "      Creating K+Sigma (numeric differentation) ..." << endl;
    } else {
        vout << "      Creating K+Sigma ..." << endl;
    }
        vout << "         Dim    = " << GPRSize << " x " << GPRSize << endl;
    CreateKS();

    CSimpleVector<double> ipos;
    ipos.CreateVector(NCVs);

    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

// ------------------------------------
    CSimpleVector<double> kff;
    kff.CreateVector(GPRSize);

    size_t nvals = NumOfValues;
    if( GlobalMinSet ){
        nvals++; // include explicitly set global minumum
    }

    CFortranMatrix  Kr;
    Kr.CreateMatrix(GPRSize,nvals);

        vout << "      Constructing kff ..." << endl;
        vout << "         Dim    = " << GPRSize << " x " << nvals << endl;

    #pragma omp parallel for firstprivate(ipos,kff)
    for(size_t indi=0; indi < NumOfValues; indi++){
        size_t i = ValueMap[indi];
        Accumulator->GetPoint(i,ipos);
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
        Accumulator->GetPoint(i,ipos);
        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            Accumulator->GetPoint(j,jpos);
            Cov[indi][indj] = GetKernelValue(ipos,jpos);

        }
    }
    if( GlobalMinSet ){
        #pragma omp parallel for firstprivate(ipos)
        for(size_t indi=0; indi < NumOfValues; indi++){
            size_t i = ValueMap[indi];
            Accumulator->GetPoint(i,ipos);
            Cov[indi][nvals-1] = GetKernelValue(ipos,GPos);
            Cov[nvals-1][indi] = Cov[indi][nvals-1];
        }
        Cov[nvals-1][nvals-1] = GetKernelValue(GPos,GPos);
    }

    RunBlasLapackPar();
    CSciBlas::gemm(-1.0,'T',Kl,'N',Kr,1.0,Cov);

//    for(size_t indi=0; indi < NumOfValues; indi++){
//        for(size_t indj=0; indj < NumOfValues; indj++){
//            cout << Cov[indi][indj] << " ";
//        }
//        cout << endl;
//    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::ReduceFES(const std::vector<bool>& keepcvs,double temp,CEnergySurface* p_rsurf)
{
    if( p_rsurf == NULL ){
        ES_ERROR("p_rsurf == NULL");
        return(false);
    }
    if( keepcvs.size() != NCVs ){
        RUNTIME_ERROR("keepcvs.size() != NCVs");
    }

    bool only_variances = false;

    p_rsurf->Allocate(FES,keepcvs);
    p_rsurf->Clear();

    CSimpleVector<int>    midx;
    midx.CreateVector(NCVs);

    CSimpleVector<int>    ridx;
    ridx.CreateVector(p_rsurf->GetNumberOfCoords());

    const double R = 1.98720425864083e-3;

    CSimpleVector<size_t>  IdxMap;
    IdxMap.CreateVector(NumOfValues);

// calculate weights
    for(size_t indi=0; indi < NumOfValues; indi++){
        size_t mbin = ValueMap[indi];
        double ene = FES->GetEnergy(mbin);
        FES->GetIPoint(mbin,midx);
        FES->ReduceIPoint(keepcvs,midx,ridx);
        size_t rbin = p_rsurf->IPoint2Bin(ridx);
        IdxMap[indi] = rbin;
        double w = exp(-ene/(R*temp));
        p_rsurf->SetEnergy(rbin,p_rsurf->GetEnergy(rbin) + w);
        p_rsurf->SetNumOfSamples(rbin,1);
    }

// calculate errors
    for(size_t rbin = 0; rbin < (size_t)p_rsurf->GetNumberOfPoints(); rbin++){
        double err = 0.0;
        // err is now variance
        for(size_t indi=0; indi < NumOfValues; indi++){
            if( IdxMap[indi] != rbin ) continue;
            size_t mbini = ValueMap[indi];
            double enei = FES->GetEnergy(mbini);
            double wi = exp(-enei/(R*temp));
            for(size_t indj=0; indj < NumOfValues; indj++){
                if( IdxMap[indj] != rbin ) continue;
                if( only_variances ){
                    if( indi != indj ) continue;
                }
                size_t mbinj = ValueMap[indj];
                double enej = FES->GetEnergy(mbinj);
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
    for(size_t rbin = 0; rbin < (size_t)p_rsurf->GetNumberOfPoints(); rbin++){
        double w = p_rsurf->GetEnergy(rbin);
        p_rsurf->SetEnergy(rbin,-R*temp*log(w));
    }

// move global minimum
   double gmin = p_rsurf->GetGlobalMinimumValue();
   p_rsurf->ApplyOffset(-gmin);

    return(true);
}

//------------------------------------------------------------------------------
