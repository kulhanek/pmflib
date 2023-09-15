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
    NumOfValues         = 0;

    IncludeError        = false;
}

//------------------------------------------------------------------------------

CSmootherGPR::~CSmootherGPR(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSmootherGPR::SetAccumulator(CPMFAccumulatorPtr accu)
{
    if( accu == NULL ) return;                 // no-accu

    CGPRKernel::SetAccumulator(accu);

    NumOfSigmaF2 = 1;
    NumOfCoVar   = 0;
    NumOfNCorr   = 1;
    NumOfSigmaN2 = 1;
}

//------------------------------------------------------------------------------

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

    if( Accu == NULL ){
        SetAccumulator(p_prx->GetAccu());
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

void CSmootherGPR::SetIncludeError(bool set)
{
    IncludeError = set;
}

//------------------------------------------------------------------------------

void CSmootherGPR::PrepForMFInfo(void)
{
    IncludeError = true;
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
    if( SigmaN2.GetLength() == 0 ){
        RUNTIME_ERROR("sigman2 is not set");
    }

    // kernel setup
    SetupKernel();

    // number of data points
    GPRSize = 0;
    for(size_t i=0; i < EneProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( EneProxyItems[i]->GetNSamples(ibin) > 0 ) GPRSize++;
        }
    }

    // create sampled map
    SampledMap.resize(GPRSize);
    EneProxyMap.resize(GPRSize);
    size_t ind = 0;
    for(size_t i=0; i < EneProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( EneProxyItems[i]->GetNSamples(ibin) <= 0 ) continue;
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
    for(size_t k=0; k < NumOfSigmaF2; k++ ){
        vout << format("      SigmaF2#%-2d= %10.4f")%(k+1)%SigmaF2[k] << endl;
    }
    for(size_t k=0; k < NumOfCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(size_t k=0; k < NumOfNCorr; k++ ){
        vout << format("      NCorr#%-2d  = %10.4f")%(k+1)%NCorr[k] << endl;
    }
    for(size_t k=0; k < NumOfSigmaN2; k++ ){
        vout << format("      SigmaN2#%-2d= %10.4e")%(k+1)%SigmaN2[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    if( ! nostat ){
        // and log of marginal likelihood
        vout << "      logML     = " << setprecision(5) << GetLogML() << endl;
        if( NeedInv || UseInv ){
            // and log of pseudo-likelihood
            vout << "      logPL     = " << setprecision(5) << GetLogPL() << endl;
        }
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

bool CSmootherGPR::WriteMFInfo(const CSmallString& name)
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
    for(size_t indi=0; indi < GPRSize; indi++){
        CEnergyProxyPtr  item = EneProxyItems[EneProxyMap[indi]];
        size_t           ibin = SampledMap[indi];

        EneSurface->GetPoint(ibin,ipos);
        mfi[indi]  = item->GetValue(ibin,E_PROXY_VALUE);
        mfie[indi] = item->GetValue(ibin,E_PROXY_ERROR);    // sigma
        mfp[indi]  = GetValue(ipos);
        mfpe[indi] = 0.0;
        for(size_t indv=0; indv < NumOfValues; indv++){
            if( ValueMap[indv] == ibin ){
                mfpe[indi] = sqrt(Cov[indv][indv]);  // convert to sigma
            }
        }
    }

    ofstream ofs(name);
    if( ! ofs ){
        CSmallString error;
        error << "unable to open file '" << name << "' for energy";
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

        ofs << format(" %20.16f %20.16f %20.16f %20.16f")%mfi[indi]%mfie[indi]%mfp[indi]%mfpe[indi];

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
            KS[indi][indj] = SigmaF2[0]*GetKernelValue(ipos,jpos);
        }
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        size_t          ibin = SampledMap[indi];
        CEnergyProxyPtr item = EneProxyItems[EneProxyMap[indi]];
        double er = item->GetValue(ibin,E_PROXY_ERROR);
        // use only sigmaN2[0]
        KS[indi][indi] += er*er*NCorr[0] + SigmaN2[0];
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
        kff[indj] = SigmaF2[0]*GetKernelValue(ip,jpos);
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
            if( EneProxyItems[i]->GetNSamples(ibin) <= 0 ) continue;
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
            int nsamples = EneProxyItems[i]->GetNSamples(ibin);
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
    if( EneSurface->IsGlobalMinSet() ){

        CSimpleVector<double> gpos;

        gpos = EneSurface->GetGlobalMinPos();
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << endl;

        EneSurface->FindGlobalMinBin();


        gpos = EneSurface->GetGlobalMinPos();
        vout << "      Closest bin found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }

        double glb_min = EneSurface->GetGlobalMinEnergy();
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->SetEnergy(j,EneSurface->GetEnergy(j)-glb_min);
        }
    } else {
        // search for global minimum
        EneSurface->FindGlobalMin();

        double                glb_min = EneSurface->GetGlobalMinEnergy();
        CSimpleVector<double> gpos    = EneSurface->GetGlobalMinPos();

        vout << "      Global minimum found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurface->SetEnergy(j,EneSurface->GetEnergy(j)-glb_min);
        }
    }

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
    size_t iglb_bin = EneSurface->GetGlobalMinBin();
    size_t iglb = 0;

    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        if( j == iglb_bin ){
            iglb = indj;
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
            Cov[indi][indj] = SigmaF2[0]*GetKernelValue(ipos,jpos);

        }
    }

    RunBlasLapackPar();
    CSciBlas::gemm(-1.0,'T',Kl,'N',Kr,1.0,Cov);

    KS = KSInv;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

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
            CalcKderWRTSigmaF2();
        } else if( prm == 1 ){
            CalcKderWRTNCorr();
        } else if ( (prm >=2) && (prm < NumOfCVs+2) ) {
            size_t cv = prm - 2;
            CalcKderWRTWFac(cv);
        } else {
            size_t cv = prm - (NumOfCVs+2);
            CalcKderWRTSigmaN2(cv);
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

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

            Kder[indi][indj] = GetKernelValue(ipos,jpos);
        }
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
                Kder[indi][indj] = SigmaF2[0]*GetKernelValueWFacDer(ipos,jpos,cv);
            } else {
                // FIXME - check validity - this avoids division by zero in GetKernelValueWFacDer
                Kder[indi][indj] = 0.0;
            }
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
        Kder[indi][indi] = er*er;
    }
}

//------------------------------------------------------------------------------

void CSmootherGPR::CalcKderWRTSigmaN2(size_t cv)
{
    Kder.SetZero();

    // ignore cv index

    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        Kder[indi][indi] = 1.0;
    }
}

//------------------------------------------------------------------------------
