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

    ConstrainedTK       = false;
    IncludeError        = false;
    NoEnergy            = false;
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

void CGHSIntegratorGPR0A::SetAccumulator(CPMFAccumulatorPtr accu)
{
    if( accu == NULL ) return;                 // no-accu

    CGPRKernel::SetAccumulator(accu);

    NumOfSigmaF2 = 3;
    NumOfCoVar   = 3;
    NumOfNCorr   = 0;
    NumOfSigmaN2 = 3*NumOfCVs;
}

//------------------------------------------------------------------------------

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

void CGHSIntegratorGPR0A::SetNoEnergy(bool set)
{
    NoEnergy = set;
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

    if( Accu == NULL ){
        RUNTIME_ERROR("no input data - Accu");
    }
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

    SetupKernel();

    // number of data points
    NumOfUsedBins = 0;
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( Accu->GetNumOfSamples(ibin) > 0 ) NumOfUsedBins++;
    }
    GPRSize = 3 * NumOfUsedBins * NumOfCVs;

    // create sampled map
    SampledMap.resize(NumOfUsedBins);
    size_t ind = 0;
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( Accu->GetNumOfSamples(ibin) <= 0 ) continue;
        SampledMap[ind] = ibin;
        ind++;
    }

    // init GPR arrays
    GPRModel.CreateVector(GPRSize);
    Y.CreateVector(GPRSize);
    KS.CreateMatrix(GPRSize,GPRSize);
    TK.CreateMatrix(3,3);
    TKder.CreateMatrix(3,3);

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
        // and log of marginal likelihood
        vout << "      logML     = " << setprecision(5) << GetLogML() << endl;
        if( NeedInv || UseInv ){
            // and log of pseudo-likelihood
            vout << "      logPL     = " << setprecision(5) << GetLogPL() << endl;
        }
        vout << "      >>>>>>>>" << endl;
        // and finally some statistics
        for(size_t k=0; k < NumOfCVs; k++ ){
            vout << "      dG/dx RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k,0) << endl;
        }
        vout << "      >>>>>>>>" << endl;
        for(size_t k=0; k < NumOfCVs; k++ ){
            vout << "      dH/dx RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k,1) << endl;
        }
        vout << "      >>>>>>>>" << endl;
        for(size_t k=0; k < NumOfCVs; k++ ){
            vout << "    -TdS/dx RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k,2) << endl;
        }
    }

    // finalize EneSurface if requested
    if( ! NoEnergy ){
        CalculateEnergy(vout);
        if( IncludeError ){
            if( FastErrors ){
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

bool CGHSIntegratorGPR0A::TrainGP(CVerboseStr& vout)
{
    if( IsNumDiffEnabled() ) {
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

// construct TK and KS
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
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            GetKernelDerIJ(ipos,jpos,kblock);

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
        Accu->GetPoint(jbin,jpos);

        GetKernelDerJ(ip,jpos,kder);

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff[offset0 + indj*NumOfCVs+jj] = TK[0][task]*kder[jj];
            kff[offset1 + indj*NumOfCVs+jj] = TK[1][task]*kder[jj];
            kff[offset2 + indj*NumOfCVs+jj] = TK[2][task]*kder[jj];
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& kff2,int task)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];
        Accu->GetPoint(jbin,jpos);

        GetKernelDerIJ(ip,jpos,kblock);

        // distribute to vector
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kff2[offset0 + indj*NumOfCVs+jj] = TK[0][task]*kblock[icoord][jj];
            kff2[offset1 + indj*NumOfCVs+jj] = TK[1][task]*kblock[icoord][jj];
            kff2[offset2 + indj*NumOfCVs+jj] = TK[2][task]*kblock[icoord][jj];
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPR0A::CalculateEnergy(CVerboseStr& vout)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("no input data - Accu");
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

    vout << "   Calculating dG, dH, -TdS surfaces ..." << endl;
    vout << "      >>>>>>>>" << endl;

// calculate energies
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

// basic EneSurface update
    #pragma omp parallel for
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        int nsamples = Accu->GetNumOfSamples(ibin);
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
        Accu->GetPoint(jbin,jpos);
        GSurface->SetEnergy(jbin,GetValue(jpos,0));
        HSurface->SetEnergy(jbin,GetValue(jpos,1));
        SSurface->SetEnergy(jbin,GetValue(jpos,2));
    }

// update FES
    if( GSurface->IsGlobalMinSet() ){

        CSimpleVector<double> gpos;

        gpos = GSurface->GetGlobalMinPos();
        vout << "      dG(x) Global minimum provided at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << endl;

        GSurface->FindGlobalMinBin();

        gpos = GSurface->GetGlobalMinPos();
        vout << "      Closest bin found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }

        double glb_min = GSurface->GetGlobalMinEnergy();
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];
            GSurface->SetEnergy(j,GSurface->GetEnergy(j)-glb_min);
        }
    } else {
        // search for global minimum
        GSurface->FindGlobalMin();

        double                glb_min = GSurface->GetGlobalMinEnergy();
        CSimpleVector<double> gpos    = GSurface->GetGlobalMinPos();

        vout << "      dG(x) Global minimum found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];
            GSurface->SetEnergy(j,GSurface->GetEnergy(j)-glb_min);
        }
    }

    int gmin_bin = GSurface->GetGlobalMinBin();

    vout << "      dG(x) SigmaF2   = " << setprecision(5) << GSurface->GetSigmaF2() << endl;
    vout << "      dG(x) SigmaF    = " << setprecision(5) << GSurface->GetSigmaF() << endl;

    vout << "      >>>>>>>>" << endl;
    double glb_min = HSurface->GetEnergy(gmin_bin);
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        int jbin = SampledMap[indj];
        HSurface->SetEnergy(jbin,HSurface->GetEnergy(jbin)-glb_min);
    }
    vout << "      dH(x) SigmaF2   = " << setprecision(5) << HSurface->GetSigmaF2() << endl;
    vout << "      dH(x) SigmaF    = " << setprecision(5) << HSurface->GetSigmaF() << endl;

    vout << "      >>>>>>>>" << endl;
    glb_min = SSurface->GetEnergy(gmin_bin);
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        int jbin = SampledMap[indj];
        SSurface->SetEnergy(jbin,SSurface->GetEnergy(jbin)-glb_min);
    }
    vout << "    -TdS(x) SigmaF2   = " << setprecision(5) << SSurface->GetSigmaF2() << endl;
    vout << "    -TdS(x) SigmaF    = " << setprecision(5) << SSurface->GetSigmaF() << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGHSIntegratorGPR0A::GetValue(const CSimpleVector<double>& position,int task)
{
    CSimpleVector<double>   kff;
    kff.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(position,kff,task);
    return( CSciBlas::dot(kff,GPRModel) );
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPR0A::GetVar(CSimpleVector<double>& lpos,int task)
{
    if( KSInverted != true ) {
        RUNTIME_ERROR("KS must be inverted!");
    }

    CSimpleVector<double>   kff;
    CSimpleVector<double>   ik;

    kff.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(lpos,kff,task);
    CSciBlas::gemv(1.0,KS,kff,0.0,ik);
    double cov = TK[task][task]*GetKernelValue(lpos,lpos) - CSciBlas::dot(kff,ik);
    return(cov);
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& lrcov,double& rrvar,int task)
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

    CreateKff(rpos,kffr,task);
    CreateKff(lpos,kffl,task);

    CSciBlas::gemv(1.0,KS,kffr,0.0,ik);

    lrcov = TK[task][task]*GetKernelValue(lpos,rpos) - CSciBlas::dot(kffl,ik);
    rrvar = TK[task][task]*GetKernelValue(rpos,rpos) - CSciBlas::dot(kffr,ik);
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

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("GetLogMLDerivatives requires K+Sigma inverted matrix");
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
        // 0<3; 3<6; 6<6+NumOfCVs; 6+NumOfCVs < 3*NumOfCVs + 6+NumOfCVs
        if( (prm >= 0) && (prm < 3) ){
            // sigmaf2
            size_t idx = prm - 0;
            CalcKderWRTSigmaF2(idx);
        } else if( (prm >= 3) && (prm < 6) ){
            // covar
            size_t idx = prm - 3;
            CalcKderWRTCoVar(idx);
        } else if( (prm >= 6) && (prm < 6+NumOfCVs) ){
            // wfac
            size_t cv = prm - 6;
            CalcKderWRTWFac(cv);
        } else if( (prm >= 6+NumOfCVs) && (prm < 3*NumOfCVs + 6+NumOfCVs) ){
            // sigman2
            size_t idx = prm - (6+NumOfCVs);
            CalcKderWRTSigmaN2(idx);
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

void CGHSIntegratorGPR0A::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    if( ! (NeedInv || UseInv) ){
        RUNTIME_ERROR("GetLogPLDerivatives requires K+Sigma inverted matrix");
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
        // 0<3; 3<6; 6<6+NumOfCVs; 6+NumOfCVs < 3*NumOfCVs + 6+NumOfCVs
        if( (prm >= 0) && (prm < 3) ){
            // sigmaf2
            size_t idx = prm - 0;
            CalcKderWRTSigmaF2(idx);
        } else if( (prm >= 3) && (prm < 6) ){
            // covar
            size_t idx = prm - 3;
            CalcKderWRTCoVar(idx);
        } else if( (prm >= 6) && (prm < 6+NumOfCVs) ){
            // wfac
            size_t cv = prm - 6;
            CalcKderWRTWFac(cv);
        } else if( (prm >= 6+NumOfCVs) && (prm < 3*NumOfCVs + 6+NumOfCVs) ){
            // sigman2
            size_t idx = prm - (6+NumOfCVs);
            CalcKderWRTSigmaN2(idx);
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

void CGHSIntegratorGPR0A::CalcKderWRTSigmaF2(size_t idx)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    CreateTKDerSigmaF2(idx);

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            GetKernelDerIJ(ipos,jpos,kblock);

            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    Kder[offset0 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TKder[0][0]*kblock[ii][jj];
                    Kder[offset0 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TKder[0][1]*kblock[ii][jj];
                    Kder[offset0 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TKder[0][2]*kblock[ii][jj];

                    Kder[offset1 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TKder[1][0]*kblock[ii][jj];
                    Kder[offset1 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TKder[1][1]*kblock[ii][jj];
                    Kder[offset1 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TKder[1][2]*kblock[ii][jj];

                    Kder[offset2 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TKder[2][0]*kblock[ii][jj];
                    Kder[offset2 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TKder[2][1]*kblock[ii][jj];
                    Kder[offset2 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TKder[2][2]*kblock[ii][jj];
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CreateTKDerSigmaF2(size_t idx)
{
    TKder.SetZero();
    TKder[idx][idx] = 1.0;

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

        CSciBlas::gemm(1.0,A,TKder,0.0,B);
        CSciBlas::gemm(1.0,'N',B,'T',A,0.0,TKder);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalcKderWRTCoVar(size_t idx)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    CreateTKDerCoVar(idx);

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            GetKernelDerIJ(ipos,jpos,kblock);

            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    Kder[offset0 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TKder[0][0]*kblock[ii][jj];
                    Kder[offset0 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TKder[0][1]*kblock[ii][jj];
                    Kder[offset0 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TKder[0][2]*kblock[ii][jj];

                    Kder[offset1 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TKder[1][0]*kblock[ii][jj];
                    Kder[offset1 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TKder[1][1]*kblock[ii][jj];
                    Kder[offset1 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TKder[1][2]*kblock[ii][jj];

                    Kder[offset2 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TKder[2][0]*kblock[ii][jj];
                    Kder[offset2 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TKder[2][1]*kblock[ii][jj];
                    Kder[offset2 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TKder[2][2]*kblock[ii][jj];
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CreateTKDerCoVar(size_t idx)
{
    TKder.SetZero();

    if( idx == 0 ) {
        TKder[0][1] = 1.0;
        TKder[1][0] = 1.0;
    }
    if( idx == 1 ) {
        TKder[0][2] = 1.0;
        TKder[2][0] = 1.0;
    }
    if( idx == 2 ) {
        TKder[1][2] = 1.0;
        TKder[2][1] = 1.0;
    }

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

        CSciBlas::gemm(1.0,A,TKder,0.0,B);
        CSciBlas::gemm(1.0,'N',B,'T',A,0.0,TKder);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalcKderWRTWFac(size_t cv)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    int offset0 = 0*NumOfUsedBins*NumOfCVs;
    int offset1 = 1*NumOfUsedBins*NumOfCVs;
    int offset2 = 2*NumOfUsedBins*NumOfCVs;

    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            GetKernelDerIJWFacDer(ipos,jpos,cv,kblock);

            // distribute to main kernel matrix
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    Kder[offset0 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TK[0][0]*kblock[ii][jj];
                    Kder[offset0 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TK[0][1]*kblock[ii][jj];
                    Kder[offset0 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TK[0][2]*kblock[ii][jj];

                    Kder[offset1 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TK[1][0]*kblock[ii][jj];
                    Kder[offset1 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TK[1][1]*kblock[ii][jj];
                    Kder[offset1 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TK[1][2]*kblock[ii][jj];

                    Kder[offset2 + indi*NumOfCVs+ii][offset0 + indj*NumOfCVs+jj] = TK[2][0]*kblock[ii][jj];
                    Kder[offset2 + indi*NumOfCVs+ii][offset1 + indj*NumOfCVs+jj] = TK[2][1]*kblock[ii][jj];
                    Kder[offset2 + indi*NumOfCVs+ii][offset2 + indj*NumOfCVs+jj] = TK[2][2]*kblock[ii][jj];
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalcKderWRTSigmaN2(size_t idx)
{
    Kder.SetZero();

    #pragma omp parallel for
    for(size_t task=0; task < 3; task++){
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            for(size_t ii=0; ii < NumOfCVs; ii++){
                if( (task*NumOfCVs + ii) == idx ){
                    Kder[task*NumOfUsedBins+indi*NumOfCVs+ii][task*NumOfUsedBins+indi*NumOfCVs+ii] = 1.0;
                }
            }
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPR0A::PrepForMFInfo(void)
{
    NeedInv = true;
}

//------------------------------------------------------------------------------

bool CGHSIntegratorGPR0A::WriteMFInfo(const CSmallString& name,int task)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    if( (task < 0) || (task >= 3) ){
        ES_ERROR("out-of-range task");
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

    CEnergyDerProxyPtr proxy = NULL;

    switch(task){
        case(0):
            proxy = GDerProxy;
        break;
        case(1):
            proxy = HDerProxy;
        break;
        case(2):
            proxy = SDerProxy;
        break;
    }

    if( proxy == NULL ){
        RUNTIME_ERROR("no proxy");
    }

    // calculate
    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];

        Accu->GetPoint(ibin,ipos);
        for(size_t k=0; k < NumOfCVs; k++){
            mfi[indi*NumOfCVs+k] = proxy->GetValue(ibin,k,E_PROXY_VALUE);
            double mfe = proxy->GetValue(ibin,k,E_PROXY_ERROR);
            mfie[indi*NumOfCVs+k] = mfe;            // this is a sigma

            mfp[indi*NumOfCVs+k] = GetMeanForce(ipos,k,task);
            double mfv           = GetMeanForceVar(ipos,k,task);   // this is a variance, sigma^2
            mfpe[indi*NumOfCVs+k] = sqrt(mfv);                  // get a sigma
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
        size_t ibin = SampledMap[indi];

        Accu->GetPoint(ibin,ipos);

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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGHSIntegratorGPR0A::GetMeanForce(const CSimpleVector<double>& position,size_t icoord,int task)
{
    CSimpleVector<double>   kff2;
    kff2.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff2(position,icoord,kff2,task);
    double mf = CSciBlas::dot(kff2,GPRModel);
    return(mf);
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPR0A::GetMeanForceVar(const CSimpleVector<double>& position,size_t icoord,int task)
{
    if( KSInverted != true ) {
        RUNTIME_ERROR("KS must be inverted!");
    }

    CSimpleVector<double>   kff2;
    CSimpleVector<double>   ik;

    kff2.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff2(position,icoord,kff2,task);
    CSciBlas::gemv(1.0,KS,kff2,0.0,ik);     // KS must be inverted

    CFortranMatrix kblock;
    kblock.CreateMatrix(NumOfCVs,NumOfCVs);

    GetKernelDerIJ(position,position,kblock);

    double var = TK[task][task]*kblock[icoord][icoord] - CSciBlas::dot(kff2,ik);
    return(var);
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPR0A::GetRMSR(size_t cv,int task)
{
    if( NumOfBins == 0 ){
        RUNTIME_ERROR("number of bins is not > 0");
    }

    CEnergyDerProxyPtr proxy = NULL;

    switch(task){
        case(0):
            proxy = GDerProxy;
        break;
        case(1):
            proxy = HDerProxy;
        break;
        case(2):
            proxy = SDerProxy;
        break;
    }

    if( proxy == NULL ){
        RUNTIME_ERROR("no proxy");
    }

    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

    double rmsr = 0.0;

    #pragma omp parallel for firstprivate(ipos) reduction(+:rmsr)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t  ibin = SampledMap[indi];

        Accu->GetPoint(ibin,ipos);

        double mfi = proxy->GetValue(ibin,cv,E_PROXY_VALUE);
        double mfp = GetMeanForce(ipos,cv,task);
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPR0A::CalculateErrors(CVerboseStr& vout)
{
    vout << "   Calculating dG(x) error ..." << endl;
    CalculateErrors(vout,0);
    vout << "      gG(x) RMSError  = " << setprecision(5) << GSurface->GetRMSError() << endl;
    vout << "      dG(x) MaxError  = " << setprecision(5) << GSurface->GetMaxError() << endl;

    vout << "   Calculating dH(x) error ..." << endl;
    CalculateErrors(vout,1);
    vout << "      dH(x) RMSError  = " << setprecision(5) << HSurface->GetRMSError() << endl;
    vout << "      dH(x) MaxError  = " << setprecision(5) << HSurface->GetMaxError() << endl;

    vout << "   Calculating -TdS(x) error ..." << endl;
    CalculateErrors(vout,2);
    vout << "    -TdS(x) RMSError  = " << setprecision(5) << SSurface->GetRMSError() << endl;
    vout << "    -TdS(x) MaxError  = " << setprecision(5) << SSurface->GetMaxError() << endl;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalculateErrors(CVerboseStr& vout,int task)
{
    if( NumOfUsedBins == 0 ){
        RUNTIME_ERROR("NumOfUsedBins == 0");
    }

    CEnergySurfacePtr surf = NULL;

    switch(task){
        case(0):
            surf = GSurface;
        break;
        case(1):
            surf = HSurface;
        break;
        case(2):
            surf = SSurface;
        break;
    }

    if( surf == NULL ){
        RUNTIME_ERROR("no surf");
    }

    CSmallTime st;
    st.GetActualTime();

    CSimpleVector<double> gpos = GSurface->GetGlobalMinPos();

    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    double  vargp = GetVar(gpos,task);
    int     nbatches = 0;

    #pragma omp parallel for firstprivate(jpos)
    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t j = SampledMap[indj];

        Accu->GetPoint(j,jpos);

        double varfc,covfg;
        GetCovVar(gpos,jpos,covfg,varfc,task);

        double error = varfc + vargp - 2.0*covfg;
        if( error > 0 ){
            error = sqrt(error);
        } else {
            error = 0.0;
        }
        surf->SetError(j,error);

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
                int comp = nbatches*100 / NumOfUsedBins;
                vout << format("      completed %6d/%6d - %2d%%")%nbatches%NumOfUsedBins%comp << endl;
                st = ct;
            }
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalculateErrorsFromCov(CVerboseStr& vout)
{
    vout << "   Calculating dG(x) error ..." << endl;
    CalculateCovs(vout,0);
    CalculateErrorsFromCov(vout,0);
    vout << "         gG(x) RMSError  = " << setprecision(5) << GSurface->GetRMSError() << endl;
    vout << "         dG(x) MaxError  = " << setprecision(5) << GSurface->GetMaxError() << endl;

    vout << "   Calculating dH(x) error ..." << endl;
    CalculateCovs(vout,1);
    CalculateErrorsFromCov(vout,1);
    vout << "         dH(x) RMSError  = " << setprecision(5) << HSurface->GetRMSError() << endl;
    vout << "         dH(x) MaxError  = " << setprecision(5) << HSurface->GetMaxError() << endl;

    vout << "   Calculating -TdS(x) error ..." << endl;
    CalculateCovs(vout,2);
    CalculateErrorsFromCov(vout,2);
    vout << "       -TdS(x) RMSError  = " << setprecision(5) << SSurface->GetRMSError() << endl;
    vout << "       -TdS(x) MaxError  = " << setprecision(5) << SSurface->GetMaxError() << endl;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalculateErrorsFromCov(CVerboseStr& vout,int task)
{
    if( NumOfUsedBins == 0 ){
        RUNTIME_ERROR("NumOfUsedBins == 0");
    }

    CEnergySurfacePtr surf = NULL;

    switch(task){
        case(0):
            surf = GSurface;
        break;
        case(1):
            surf = HSurface;
        break;
        case(2):
            surf = SSurface;
        break;
    }

    if( surf == NULL ){
        RUNTIME_ERROR("no surf");
    }

    // global minimum
    size_t iglb = GSurface->GetGlobalMinBin();

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t j = SampledMap[indj];

        double varfc = Cov[indj][indj];
        double covfg = Cov[indj][iglb];
        double vargp = Cov[iglb][iglb];

        double error = varfc + vargp - 2.0*covfg;
        if( error > 0 ){
            error = sqrt(error);
        } else {
            error = 0.0;
        }
        surf->SetError(j,error);
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPR0A::CalculateCovs(CVerboseStr& vout,int task)
{
    if( NumOfUsedBins == 0 ){
        RUNTIME_ERROR("NumOfUsedBins == 0");
    }
    if( GPRSize == 0 ){
        RUNTIME_ERROR("GPRSize == 0");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("NumOfCVs == 0");
    }

// ------------------------------------
        vout << "      Calculating covariances ..." << endl;
    if( IsNumDiffEnabled() ) {
        vout << "         Creating K+Sigma (numeric differentation) ..." << endl;
    } else {
        vout << "         Creating K+Sigma ..." << endl;
    }
        vout << "            Dim    = " << GPRSize << " x " << GPRSize << endl;

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

    size_t nvals = NumOfUsedBins;

    CFortranMatrix  Kr;
    Kr.CreateMatrix(GPRSize,nvals);

        vout << "         Constructing kff ..." << endl;
        vout << "            Dim    = " << GPRSize << " x " << nvals << endl;

    #pragma omp parallel for firstprivate(ipos,kff)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        Accu->GetPoint(i,ipos);
        CreateKff(ipos,kff,task);
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
            vout << format("         Solving (K+Sigma)^(-1)*kff by LU ...") << endl;
            result = CSciLapack::solvleLU(KS,Kr);
            if( result != 0 ) return;
            break;
        case(EGPRLA_LL):
            vout << format("         Solving (K+Sigma)^(-1)*kff by LL ...") << endl;
            result = CSciLapack::solvleLL(KS,Kr);
            if( result != 0 ) return;
            break;
            break;
    default:
        INVALID_ARGUMENT("unsupported method");
    }

        vout << "         Calculating Cov ..." << endl;
        vout << "             Dim    = " << nvals << " x " << nvals << endl;

// ------------------------------------
    Cov.CreateMatrix(nvals,nvals);

    #pragma omp parallel for firstprivate(ipos,jpos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        Accu->GetPoint(i,ipos);
        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t j = SampledMap[indj];
            Accu->GetPoint(j,jpos);
            Cov[indi][indj] = TK[task][task]*GetKernelValue(ipos,jpos);

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


