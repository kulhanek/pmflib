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

#include <GHSIntegratorGPRcC.hpp>
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

CGHSIntegratorGPRcC::CGHSIntegratorGPRcC(void)
{
    GDerProxy           = NULL;
    HEneProxy           = NULL;
    SDerProxy           = NULL;
    GSurface            = NULL;
    HSurface            = NULL;
    SSurface            = NULL;

    GPRSize             = 0;
    NumOfUsedBins       = 0;

    NoEnergy            = false;

    KSInverted          = false;

    HMean               = 0.0;
}

//------------------------------------------------------------------------------

CGHSIntegratorGPRcC::~CGHSIntegratorGPRcC(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPRcC::SetAccumulator(CPMFAccumulatorPtr accu)
{
    if( accu == NULL ) return;                 // no-accu

    CGPRKernel::SetAccumulator(accu);

    NumOfSigmaF2 = 3;
    NumOfCoVar   = 0;
    NumOfNCorr   = 0;
    NumOfSigmaN2 = 3;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::SetGDerProxy(CEnergyDerProxyPtr p_proxy)
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

void CGHSIntegratorGPRcC::SetHEneProxy(CEnergyProxyPtr p_proxy)
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

    HEneProxy = p_proxy;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::SetSDerProxy(CEnergyDerProxyPtr p_proxy)
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

void CGHSIntegratorGPRcC::SetOutputFES(CEnergySurfacePtr p_surf)
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

void CGHSIntegratorGPRcC::SetOutputHES(CEnergySurfacePtr p_surf)
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

void CGHSIntegratorGPRcC::SetOutputSES(CEnergySurfacePtr p_surf)
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

void CGHSIntegratorGPRcC::SetIncludeError(bool iset)
{
    if( iset == true ){
        RUNTIME_ERROR("error calculation is not possible with this GPR integrator");
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::SetNoEnergy(bool iset)
{
    NoEnergy = iset;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::PrepForHyprmsGrd(bool iset)
{
   NeedInv |= iset;
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::SetCalcLogPL(bool iset)
{
   NeedInv |= iset;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CGHSIntegratorGPRcC::Integrate(CVerboseStr& vout,bool nostat)
{
    PrintExecInfo(vout);

    if( Accu == NULL ){
        RUNTIME_ERROR("no input data - Accu");
    }
    if( GDerProxy == NULL ){
        RUNTIME_ERROR("no input data - FES");
    }
    if( HEneProxy == NULL ){
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
        RUNTIME_ERROR("SigmaF2 is not set");
    }
    if( WFac.GetLength() != NumOfCVs ){
        RUNTIME_ERROR("WFac is not set");
    }
    if( SigmaN2.GetLength() != 3*NumOfCVs ){
        RUNTIME_ERROR("SigmaN2 is not set");
    }
    if( NumOfCVs != 1 ){
        RUNTIME_ERROR("NumOfCVs must be one");
    }

    // setup kernel
    SetupKernel();

    // number of data points
    NumOfUsedBins = 0;
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( Accu->GetNumOfSamples(ibin) > 0 ) NumOfUsedBins++;
    }
    GPRSize = 3 * NumOfUsedBins;

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

    // print hyperparameters
        vout        << "   Hyperparameters ..." << endl;
    for(size_t k=0; k < 3; k++ ){
        vout << format("      SigmaF2#%-2d= %10.4f")%(k+1)%SigmaF2[k] << endl;
    }
    for(size_t k=0; k < 1; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(size_t k=0; k < 3; k++ ){
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
        vout << "      dH RMSR     " << "    = " << setprecision(5) << GetRMSR(0,1) << endl;
        vout << "      >>>>>>>>" << endl;
        for(size_t k=0; k < NumOfCVs; k++ ){
            vout << "    -TdS/dx RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k,2) << endl;
        }
    }

    // finalize EneSurface if requested
    if( ! NoEnergy ){
        CalculateEnergy(vout);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CGHSIntegratorGPRcC::TrainGP(CVerboseStr& vout)
{
    if( IsNumDiffEnabled() ) {
        vout << "   Creating K+Sigma and Y (numeric differentiation) ..." << endl;
    } else {
        vout << "   Creating K+Sigma and Y ..." << endl;
    }
        vout << "      Kernel    = " << GetKernelName() << endl;
        vout << "      Dim       = " << GPRSize << " x " << GPRSize << endl;

// HMean
    HMean = 0.0;
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        HMean += HEneProxy->GetValue(ibin,E_PROXY_VALUE);
    }
    HMean /= (double)NumOfUsedBins;

// construct Y
    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Y[offset0 + indi] = GDerProxy->GetValue(ibin,0,E_PROXY_VALUE);
        Y[offset1 + indi] = HEneProxy->GetValue(ibin,E_PROXY_VALUE)-HMean;
        Y[offset2 + indi] = SDerProxy->GetValue(ibin,0,E_PROXY_VALUE);
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

void CGHSIntegratorGPRcC::CalculateEnergy(CVerboseStr& vout)
{
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
        HSurface->SetEnergy(jbin,GetValue(jpos,1)+HMean);
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

void CGHSIntegratorGPRcC::CreateKS(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(3,3);

    KS.SetZero();

    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    // generate main kernel block
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            CreateTK(ipos,jpos,kblock);

            // distribute to main kernel matrix

            KS[offset0 + indi][offset0 + indj] = kblock[0][0];
            KS[offset0 + indi][offset1 + indj] = kblock[0][1];
            KS[offset0 + indi][offset2 + indj] = kblock[0][2];

            KS[offset1 + indi][offset0 + indj] = kblock[1][0];
            KS[offset1 + indi][offset1 + indj] = kblock[1][1];
            KS[offset1 + indi][offset2 + indj] = kblock[1][2];

            KS[offset2 + indi][offset0 + indj] = kblock[2][0];
            KS[offset2 + indi][offset1 + indj] = kblock[2][1];
            KS[offset2 + indi][offset2 + indj] = kblock[2][2];
        }
    }

// error of data points;
    int offsetn0 = 0;
    int offsetn1 = offsetn0 + 1;
    int offsetn2 = offsetn1 + 1;

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        KS[offset0 + indi][offset0 + indi] += SigmaN2[offsetn0];
        KS[offset1 + indi][offset1 + indi] += SigmaN2[offsetn1];
        KS[offset2 + indi][offset2 + indi] += SigmaN2[offsetn2];
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CreateTK(const CSimpleVector<double>& ip, const CSimpleVector<double>& jp, CFortranMatrix& kblock)
{
    kblock.SetZero();

    CSimpleVector<double> kderi;
    CSimpleVector<double> kderj;
    kderi.CreateVector(NumOfCVs);
    kderj.CreateVector(NumOfCVs);

    CFortranMatrix kderij;
    kderij.CreateMatrix(NumOfCVs,NumOfCVs);

    double kval = GetKernelValue(ip,jp);
    GetKernelDerI(ip,jp,kderi);
    GetKernelDerJ(ip,jp,kderj);
    GetKernelDerIJ(ip,jp,kderij);

    double A = SigmaF2[0];
    double B = SigmaF2[1];
    double C = SigmaF2[2];

    kblock[0][0] = B*kval+C*kderij[0][0];
    kblock[1][0] = C*kderj[0];
    kblock[2][0] = B*kval;

    kblock[0][1] = C*kderi[0];
    kblock[1][1] = A*kval+C*kval;
    kblock[2][1] = -A*kderi[0];

    kblock[0][2] = B*kval;
    kblock[1][2] = -A*kderj[0];
    kblock[2][2] = B*kval+A*kderij[0][0];
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CreateTKInt(const CSimpleVector<double>& ip, const CSimpleVector<double>& jp, CFortranMatrix& kblock)
{
    kblock.SetZero();

    CSimpleVector<double> kderi;
    CSimpleVector<double> kderj;
    kderi.CreateVector(NumOfCVs);
    kderj.CreateVector(NumOfCVs);

    CFortranMatrix kderij;
    kderij.CreateMatrix(NumOfCVs,NumOfCVs);

    double kval = GetKernelValue(ip,jp);
    double kint = GetKernelIntI(ip,jp);

    GetKernelDerI(ip,jp,kderi);
    GetKernelDerJ(ip,jp,kderj);
    GetKernelDerIJ(ip,jp,kderij);

    double A = SigmaF2[0];
    double B = SigmaF2[1];
    double C = SigmaF2[2];

    kblock[0][0] = B*kint+C*kderj[0];
    kblock[0][1] = C*kval;
    kblock[0][2] =  B*kint;

    kblock[1][0] = C*kderj[0];
    kblock[1][1] = A*kval+C*kval;
    kblock[1][2] = -A*kderj[0];

    kblock[2][0] = B*kint;
    kblock[2][1] = -A*kval;
    kblock[2][2] = B*kint+A*kderj[0];
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGHSIntegratorGPRcC::GetValue(const CSimpleVector<double>& position,int task)
{
    CSimpleVector<double>   kff;
    kff.CreateVector(GPRSize);

    RunBlasLapackSeq();

    CreateKff(position,kff,task);
    return( CSciBlas::dot(kff,GPRModel) );
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPRcC::GetTrainingValue(const CSimpleVector<double>& position,size_t icoord,int task)
{
    if( (task < 0) || (task >= 3) ){
        ES_ERROR("out-of-range task");
        return(false);
    }

    CSimpleVector<double>   kff2;
    kff2.CreateVector(GPRSize);

    RunBlasLapackSeq();

    double mf = 0.0;

    CreateKff2(position,icoord,kff2,task);
    mf = CSciBlas::dot(kff2,GPRModel);

    if( task == 1 ){
        mf += HMean;
    }

    return(mf);
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPRcC::GetTrainingValueVar(const CSimpleVector<double>& position,size_t icoord,int task)
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

    double var = 0.0;
    var = TK[task][task]*GetKernelValue(position,position) - CSciBlas::dot(kff2,ik);
    return(var);
}

//------------------------------------------------------------------------------

double CGHSIntegratorGPRcC::GetRMSR(size_t cv,int task)
{
    if( NumOfBins == 0 ){
        RUNTIME_ERROR("number of bins is not > 0");
    }

    double rmsr = 0.0;

    if( (task == 0) || (task == 2) ){
        CEnergyDerProxyPtr proxy = NULL;

        switch(task){
            case(0):
                proxy = GDerProxy;
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

        #pragma omp parallel for firstprivate(ipos) reduction(+:rmsr)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t  ibin = SampledMap[indi];

            Accu->GetPoint(ibin,ipos);

            double mfi = proxy->GetValue(ibin,cv,E_PROXY_VALUE);
            double mfp = GetTrainingValue(ipos,cv,task);
            double diff = mfi - mfp;
            rmsr += diff*diff;
        }
    } else if (task == 1 ) {
        CEnergyProxyPtr proxy = NULL;

        switch(task){
            case(1):
                proxy = HEneProxy;
            break;
        }

        if( proxy == NULL ){
            RUNTIME_ERROR("no proxy");
        }

        CSimpleVector<double> ipos;
        ipos.CreateVector(NumOfCVs);

        #pragma omp parallel for firstprivate(ipos) reduction(+:rmsr)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t  ibin = SampledMap[indi];

            Accu->GetPoint(ibin,ipos);

            double mfi = proxy->GetValue(ibin,E_PROXY_VALUE);
            double mfp = GetTrainingValue(ipos,cv,task);
            double diff = mfi - mfp;
            rmsr += diff*diff;
        }
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

void CGHSIntegratorGPRcC::CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& kff,int task)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(3,3);

    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];
        Accu->GetPoint(jbin,jpos);

        CreateTKInt(ip,jpos,kblock);

        kff[offset0 + indj] = kblock[task][0];
        kff[offset1 + indj] = kblock[task][1];
        kff[offset2 + indj] = kblock[task][2];
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& kff2,int task)
{
    CSimpleVector<double> jpos;
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(3,3);

    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    for(size_t indj=0; indj < NumOfUsedBins; indj++){
        size_t jbin = SampledMap[indj];
        Accu->GetPoint(jbin,jpos);

        CreateTK(ip,jpos,kblock);

        kff2[offset0 + indj] = kblock[task][0];
        kff2[offset1 + indj] = kblock[task][1];
        kff2[offset2 + indj] = kblock[task][2];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPRcC::PrepForMFInfo(void)
{
    NeedInv = true;
}

//------------------------------------------------------------------------------

bool CGHSIntegratorGPRcC::WriteMFInfo(const CSmallString& name,int task)
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

    if( (task == 0) || (task == 2) ){
        CEnergyDerProxyPtr proxy = NULL;

        switch(task){
            case(0):
                proxy = GDerProxy;
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

                mfp[indi*NumOfCVs+k] = GetTrainingValue(ipos,k,task);
                double mfv           = GetTrainingValueVar(ipos,k,task);   // this is a variance, sigma^2
                mfpe[indi*NumOfCVs+k] = sqrt(mfv);                  // get a sigma
            }
        }
    } else if( task == 1 ) {
        CEnergyProxyPtr proxy = NULL;

        switch(task){
            case(1):
                proxy = HEneProxy;
            break;
        }

        if( proxy == NULL ){
            RUNTIME_ERROR("no proxy");
        }

        CSimpleVector<double> ipos;
        ipos.CreateVector(NumOfCVs);

        // calculate
        #pragma omp parallel for firstprivate(ipos)
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t ibin = SampledMap[indi];

            Accu->GetPoint(ibin,ipos);
            size_t k = 0;
            mfi[indi*NumOfCVs+k] = proxy->GetValue(ibin,E_PROXY_VALUE);
            double mfe = proxy->GetValue(ibin,E_PROXY_ERROR);
            mfie[indi*NumOfCVs+k] = mfe;            // this is a sigma

            mfp[indi*NumOfCVs+k] = GetTrainingValue(ipos,k,task);
            double mfv           = GetTrainingValueVar(ipos,k,task);   // this is a variance, sigma^2
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

        if( task == 1 ){
            size_t k = 0;
            ofs << format(" %20.16f %20.16f %20.16f %20.16f")%mfi[indi*NumOfCVs+k]%mfie[indi*NumOfCVs+k]%mfp[indi*NumOfCVs+k]%mfpe[indi*NumOfCVs+k];
        } else {
            for(size_t k=0; k < NumOfCVs; k++){
                ofs << format(" %20.16f %20.16f %20.16f %20.16f")%mfi[indi*NumOfCVs+k]%mfie[indi*NumOfCVs+k]%mfp[indi*NumOfCVs+k]%mfpe[indi*NumOfCVs+k];
            }
        }
        ofs << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGHSIntegratorGPRcC::GetLogML(void)
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

double CGHSIntegratorGPRcC::GetLogPL(void)
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

void CGHSIntegratorGPRcC::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
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
        // 0<3; 3<3+NumOfCVs; 3+NumOfCVs < 3*NumOfCVs + 3+NumOfCVs
        if( (prm >= 0) && (prm < 3) ){
            // sigmaf2
            size_t idx = prm - 0;
            CalcKderWRTSigmaF2(idx);
        } else if( (prm >= 3) && (prm < 3+NumOfCVs) ){
            // wfac
            size_t cv = prm - 3;
            CalcKderWRTWFac(cv);
        } else if( (prm >= 3+NumOfCVs) && (prm < 3*NumOfCVs + 3+NumOfCVs) ){
            // sigman2
            size_t idx = prm - (3+NumOfCVs);
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

void CGHSIntegratorGPRcC::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
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
        // 0<3; 3<3+NumOfCVs; 3+NumOfCVs < 3*NumOfCVs + 3+NumOfCVs
        if( (prm >= 0) && (prm < 3) ){
            // sigmaf2
            size_t idx = prm - 0;
            CalcKderWRTSigmaF2(idx);
        } else if( (prm >= 3) && (prm < 3+NumOfCVs) ){
            // wfac
            size_t cv = prm - 3;
            CalcKderWRTWFac(cv);
        } else if( (prm >= 3+NumOfCVs) && (prm < 3*NumOfCVs + 3+NumOfCVs) ){
            // sigman2
            size_t idx = prm - (3+NumOfCVs);
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGHSIntegratorGPRcC::CalcKderWRTSigmaF2(size_t idx)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(3,3);

    Kder.SetZero();

    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    // generate main kernel block
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            CreateTKDerSigmaF2(ipos,jpos,kblock,idx);

            Kder[offset0 + indi][offset0 + indj] = kblock[0][0];
            Kder[offset0 + indi][offset1 + indj] = kblock[0][1];
            Kder[offset0 + indi][offset2 + indj] = kblock[0][2];

            Kder[offset1 + indi][offset0 + indj] = kblock[1][0];
            Kder[offset1 + indi][offset1 + indj] = kblock[1][1];
            Kder[offset1 + indi][offset2 + indj] = kblock[1][2];

            Kder[offset2 + indi][offset0 + indj] = kblock[2][0];
            Kder[offset2 + indi][offset1 + indj] = kblock[2][1];
            Kder[offset2 + indi][offset2 + indj] = kblock[2][2];
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CreateTKDerSigmaF2(const CSimpleVector<double>& ip, const CSimpleVector<double>& jp, CFortranMatrix& kblock, size_t idx)
{
    kblock.SetZero();

    CSimpleVector<double> kderi;
    CSimpleVector<double> kderj;
    kderi.CreateVector(NumOfCVs);
    kderj.CreateVector(NumOfCVs);

    CFortranMatrix kderij;
    kderij.CreateMatrix(NumOfCVs,NumOfCVs);

    double kval = GetKernelValue(ip,jp);
    GetKernelDerI(ip,jp,kderi);
    GetKernelDerJ(ip,jp,kderj);
    GetKernelDerIJ(ip,jp,kderij);

    double A = 0;
    if( idx == 0 ) A = 1.0;

    double B = 0;
    if( idx == 1 ) B = 1.0;

    double C = 0;
    if( idx == 2 ) C = 1.0;

    kblock[0][0] = B*kval+C*kderij[0][0];
    kblock[1][0] = C*kderj[0];
    kblock[2][0] = B*kval;

    kblock[0][1] = C*kderi[0];
    kblock[1][1] = A*kval+C*kval;
    kblock[2][1] = -A*kderi[0];

    kblock[0][2] = B*kval;
    kblock[1][2] = -A*kderj[0];
    kblock[2][2] = B*kval+A*kderij[0][0];
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CalcKderWRTWFac(size_t cv)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NumOfCVs);
    jpos.CreateVector(NumOfCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(3,3);

    Kder.SetZero();

    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    // generate main kernel block
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t ibin = SampledMap[indi];
        Accu->GetPoint(ibin,ipos);

        for(size_t indj=0; indj < NumOfUsedBins; indj++){
            size_t jbin = SampledMap[indj];
            Accu->GetPoint(jbin,jpos);

            CreateTKDerWFac(ipos,jpos,kblock,cv);

            Kder[offset0 + indi][offset0 + indj] = kblock[0][0];
            Kder[offset0 + indi][offset1 + indj] = kblock[0][1];
            Kder[offset0 + indi][offset2 + indj] = kblock[0][2];

            Kder[offset1 + indi][offset0 + indj] = kblock[1][0];
            Kder[offset1 + indi][offset1 + indj] = kblock[1][1];
            Kder[offset1 + indi][offset2 + indj] = kblock[1][2];

            Kder[offset2 + indi][offset0 + indj] = kblock[2][0];
            Kder[offset2 + indi][offset1 + indj] = kblock[2][1];
            Kder[offset2 + indi][offset2 + indj] = kblock[2][2];
        }
    }
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CreateTKDerWFac(const CSimpleVector<double>& ip, const CSimpleVector<double>& jp, CFortranMatrix& kblock, int cv)
{
    kblock.SetZero();

    CSimpleVector<double> kderi;
    CSimpleVector<double> kderj;
    kderi.CreateVector(NumOfCVs);
    kderj.CreateVector(NumOfCVs);

    CFortranMatrix kderij;
    kderij.CreateMatrix(NumOfCVs,NumOfCVs);

    double kval = GetKernelValueWFacDer(ip,jp,cv);
    GetKernelDerIWFacDer(ip,jp,cv,kderi);
    GetKernelDerJWFacDer(ip,jp,cv,kderj);
    GetKernelDerIJWFacDer(ip,jp,cv,kderij);

    double A = SigmaF2[0];
    double B = SigmaF2[1];
    double C = SigmaF2[2];

    kblock[0][0] = B*kval+C*kderij[0][0];
    kblock[1][0] = C*kderj[0];
    kblock[2][0] = B*kval;

    kblock[0][1] = C*kderi[0];
    kblock[1][1] = A*kval+C*kval;
    kblock[2][1] = -A*kderi[0];

    kblock[0][2] = B*kval;
    kblock[1][2] = -A*kderj[0];
    kblock[2][2] = B*kval+A*kderij[0][0];
}

//------------------------------------------------------------------------------

void CGHSIntegratorGPRcC::CalcKderWRTSigmaN2(size_t idx)
{
    if( (idx < 0) || (idx > 2) ) RUNTIME_ERROR("idx out-of-range");

    Kder.SetZero();

    int offset0 = 0;
    int offset1 = offset0 + NumOfUsedBins;
    int offset2 = offset1 + NumOfUsedBins;

    #pragma omp parallel for
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        if( idx == 0 ) Kder[offset0 + indi][offset0 + indi] = 1.0;
        if( idx == 1 ) Kder[offset1 + indi][offset1 + indi] = 1.0;
        if( idx == 2 ) Kder[offset2 + indi][offset2 + indi] = 1.0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


