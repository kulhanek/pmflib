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

#include <ABFEnthalpyGPR.hpp>
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

CABFEnthalpyGPR::CABFEnthalpyGPR(void)
{
    Accumulator         = NULL;
    FES                 = NULL;

    GPRSize             = 0;
    NumOfUsedBins       = 0;
    NumOfValues         = 0;
    NCVs                = 0;
    NumOfBins           = 0;

    SigmaF2             = 15.0;
    NCorr               = 1.0;

    IncludeError        = false;
    GlobalMinSet        = false;
    GPosSet             = false;

    Method              = EGPRLA_LU;
    Kernel              = EGPRK_ARDSE;

    NumOfThreads        = 1;

}

//------------------------------------------------------------------------------

CABFEnthalpyGPR::~CABFEnthalpyGPR(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFEnthalpyGPR::SetInputABFAccumulator(CABFAccumulator* p_accu)
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

void CABFEnthalpyGPR::SetOutputHESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetSigmaF2(double sigf2)
{
    SigmaF2 = sigf2;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetNCorr(const CSmallString& spec)
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
    }
    NCorr = last_ncorr;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetNCorr(double value)
{
    NCorr = value;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetWFac(const CSmallString& spec)
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

void CABFEnthalpyGPR::SetWFac(CSimpleVector<double>& wfac)
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

void CABFEnthalpyGPR::SetWFac(size_t cvind, double value)
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

void CABFEnthalpyGPR::SetIncludeError(bool set)
{
    IncludeError = set;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetLAMethod(EGPRLAMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::SetLAMethod(const CSmallString& method)
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

void CABFEnthalpyGPR::SetKernel(const CSmallString& kernel)
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

void CABFEnthalpyGPR::SetGlobalMin(const CSmallString& spec)
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

void CABFEnthalpyGPR::SetGlobalMin(const CSimpleVector<double>& pos)
{
    GlobalMinSet = true;
    GPos = pos;
    GPosSet = true;
}

//------------------------------------------------------------------------------

CSimpleVector<double> CABFEnthalpyGPR::GetGlobalMin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("")
    }
    return(GPos);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFEnthalpyGPR::Interpolate(CVerboseStr& vout,bool nostat)
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
    GPRSize = NumOfUsedBins;

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
        vout << format("      NCorr     = %10.4f")%NCorr << endl;

    for(size_t k=0; k < NCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    // finalize HES
    CalculateEnergy(vout);
    if( IncludeError ){
        CalculateCovs(vout);
        CalculateErrorsFromCov(vout);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::PrintExecInfo(CVerboseStr& vout)
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

void CABFEnthalpyGPR::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFEnthalpyGPR::TrainGP(CVerboseStr& vout)
{
    vout << "   Creating K+Sigma and Y ..." << endl;
    vout << "      Kernel    = " << GetKernelName() << endl;
    vout << "      Dim       = " << GPRSize << " x " << GPRSize << endl;

    Mean = 0.0;
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        double mf = Accumulator->GetValue(i,EABF_H_VALUE);
        Mean += mf;
    }
    Mean /= (double)GPRSize;

// construct Y
    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        double mf = Accumulator->GetValue(i,EABF_H_VALUE);
        Y[indi] = mf - Mean;
    }

// construct KS
    CreateKS();

    RunBlasLapackPar();

// inverting the K+Sigma
    int result = 0;
    switch(Method){
        case(EGPRLA_LU):
                GPRModel = Y;
                vout << "   Training GPR by LU ..." << endl;
                result = CSciLapack::solvleLU(KS,GPRModel,logdetK);
                if( result != 0 ) return(false);
            break;
        case(EGPRLA_LL):
                GPRModel = Y;
                vout << "   Training GPR by LL ..." << endl;
                result = CSciLapack::solvleLL(KS,GPRModel,logdetK);
                if( result != 0 ) return(false);
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

const CSmallString CABFEnthalpyGPR::GetKernelName(void)
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

void CABFEnthalpyGPR::CreateKS(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    // main kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            KS[indi][indj] = GetKernelValue(ipos,jpos);
        }
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        double er = Accumulator->GetValue(i,EABF_H_ERROR);
        KS[indi][indi] += er*er*NCorr;
    }
}

//------------------------------------------------------------------------------

void CABFEnthalpyGPR::CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& kff)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    // main kernel matrix
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t j = SampledMap[indj];

        Accumulator->GetPoint(j,jpos);
        kff[indj] = GetKernelValue(ip,jpos);
    }

}

//------------------------------------------------------------------------------

double CABFEnthalpyGPR::GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
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

void CABFEnthalpyGPR::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating enthalpy ..." << endl;

// create map for bins with calculated energy and error
    NumOfValues = 0;
    for(size_t i=0; i < NumOfBins; i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        if( samples <= 0 ) continue;
        NumOfValues++;
    }
    ValueMap.CreateVector(NumOfValues);

    size_t indi = 0;
    for(size_t i=0; i < NumOfBins; i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        if( samples <= 0 ) continue;
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

// basic HES update
    for(size_t i=0; i < NumOfBins; i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        FES->SetNumOfSamples(i,samples);
        FES->SetEnergy(i,0.0);
        FES->SetError(i,0.0);
    }

// update HES
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
        vout << "      Offset    = " << glb_min << endl;

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
        vout << "      Offset    = " << glb_min << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            FES->SetEnergy(j,values[indj]-glb_min);
        }
    }

    GPosSet = true;

    vout << "      SigmaF2   = " << setprecision(5) << FES->GetSigmaF2() << endl;
}

//------------------------------------------------------------------------------

double CABFEnthalpyGPR::GetValue(const CSimpleVector<double>& position)
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

void CABFEnthalpyGPR::CalculateErrorsFromCov(CVerboseStr& vout)
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

void CABFEnthalpyGPR::CalculateCovs(CVerboseStr& vout)
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
    vout << "      Creating K+Sigma ..." << endl;
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
}

//------------------------------------------------------------------------------
