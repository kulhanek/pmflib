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

    IncludeError        = false;
    IncludeGluedBins    = false;
    NoEnergy            = false;
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

void CIntegratorGPR::SetAccumulator(CPMFAccumulatorPtr accu)
{
    if( accu == NULL ) return;                 // no-accu

    CGPRKernel::SetAccumulator(accu);

    NumOfSigmaF2 = 1;
    NumOfCoVar   = 0;
    NumOfNCorr   = 1;
    NumOfSigmaN2 = NumOfCVs;
}

//------------------------------------------------------------------------------

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

    if( Accu == NULL ){
        SetAccumulator(p_proxy->GetAccu());
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CIntegratorGPR::SetIncludeError(bool iset)
{
    IncludeError = iset;
    if( FastErrors == false ){
        NeedInv |= iset;
    }
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetNoEnergy(bool iset)
{
    NoEnergy = iset;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::IncludeGluedAreas(bool iset)
{
    IncludeGluedBins = iset;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::PrepForHyprmsGrd(bool iset)
{
   NeedInv |= iset;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetCalcLogPL(bool iset)
{
   NeedInv |= iset;
}

//------------------------------------------------------------------------------

void CIntegratorGPR::SetFastError(bool iset)
{
    FastErrors = iset;
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

    SetupKernel();

    // number of data points
    NumOfUsedBins = 0;
    for(size_t i=0; i < DerProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( DerProxyItems[i]->GetNSamples(ibin) > 0 ) NumOfUsedBins++;
        }
    }
    GPRSize = NumOfUsedBins * NumOfCVs;

    // create sampled map
    SampledMap.resize(NumOfUsedBins);
    DerProxyMap.resize(NumOfUsedBins);
    size_t ind = 0;
    for(size_t i=0; i < DerProxyItems.size(); i++){
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( DerProxyItems[i]->GetNSamples(ibin) <= 0 ) continue;
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
                CalculateErrors(vout);
            }
        }
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CIntegratorGPR::TrainGP(CVerboseStr& vout)
{
    if( IsNumDiffEnabled() ) {
        vout << "   Creating K+Sigma and Y (numeric differentiation) ..." << endl;
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

            GetKernelDerIJ(ipos,jpos,kblock);

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    KS[indi*NumOfCVs+ii][indj*NumOfCVs+jj] = SigmaF2[0]*kblock[ii][jj];
                }
            }
        }
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        CEnergyDerProxyPtr  item = DerProxyItems[DerProxyMap[indi]];
        size_t              ibin = SampledMap[indi];
        for(size_t ii=0; ii < NumOfCVs; ii++){
            double er = item->GetValue(ibin,ii,E_PROXY_ERROR);
            KS[indi*NumOfCVs+ii][indi*NumOfCVs+ii] += er*er*NCorr[0] + SigmaN2[ii];
        }
    }
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
        GetKernelDerJ(ip,jpos,kder);

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff[indj*NumOfCVs+jj] = SigmaF2[0]*kder[jj];
        }
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
        GetKernelDerIJ(ip,jpos,kblock);

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff2[indj*NumOfCVs+jj] = SigmaF2[0]*kblock[icoord][jj];
        }
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
            int samples = DerProxyItems[i]->GetNSamples(ibin);
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
            int nsamples = DerProxyItems[i]->GetNSamples(ibin);
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

    GetKernelDerIJ(position,position,kblock);

    // FIXME - validate
    double var = SigmaF2[0]*kblock[icoord][icoord] - CSciBlas::dot(kff2,ik);
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
    double cov = SigmaF2[0]*GetKernelValue(lpos,lpos) - CSciBlas::dot(kff,ik);
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

    lrcov = SigmaF2[0]*GetKernelValue(lpos,rpos) - CSciBlas::dot(kffl,ik);
    rrvar = SigmaF2[0]*GetKernelValue(rpos,rpos) - CSciBlas::dot(kffr,ik);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CIntegratorGPR::CalculateErrors(CVerboseStr& vout)
{
    if( NumOfValues == 0 ){
        RUNTIME_ERROR("NumOfValues == 0");
    }

    vout << "   Calculating FES error ..." << endl;
    CSmallTime st;
    st.GetActualTime();

    CSimpleVector<double> gpos = EneSurface->GetGlobalMinPos();

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
    if( IsNumDiffEnabled() ) {
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

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

            GetKernelDerIJ(ipos,jpos,kblock);

            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    // SigmaF2 cannot be zero
                    Kder[indi*NumOfCVs+ii][indj*NumOfCVs+jj] = kblock[ii][jj];
                }
            }
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

            GetKernelDerIJWFacDer(ipos,jpos,cv,kblock);

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    Kder[indi*NumOfCVs+ii][indj*NumOfCVs+jj] = SigmaF2[0]*kblock[ii][jj];
                }
            }
        }
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
