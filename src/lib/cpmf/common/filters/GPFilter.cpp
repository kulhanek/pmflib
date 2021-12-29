// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <GPFilter.hpp>
#include <iostream>
#include <SciLapack.hpp>
#include <SciBlas.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGPFilter::CGPFilter(void)
{
    NumOfThreads    = 1;
    GPRSize         = 0;
    SigmaF2         = 1.0;
    WFac            = 1.0;
    Noise           = 0.0;
    RCond           = 1e-7;
    Kernel          = EGPFK_ARDMC32;

}

//------------------------------------------------------------------------------

CGPFilter::~CGPFilter(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPFilter::RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out)
{
    // FIXME
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPFilter::SetFilter(double timestep,int framelen, double wfac, double noise,CVerboseStr& vout)
{
    SetTimeStep(timestep);

    PrintExecInfo(vout);

    vout << "   Creating K+Sigma ..." << endl;

    if( framelen <= 0 ) {
        RUNTIME_ERROR("framelen has to be greater than zero");
    }

    GPRSize = framelen;
    WFac    = wfac;
    Noise   = noise;

    vout << "      Kernel    = " << GetKernelName() << endl;
    vout << "      Dim       = " << GPRSize << " x " << GPRSize << endl;
    vout << "      WFac      = " << WFac << endl;
    vout << "      Noise     = " << Noise << endl;

    GPRModel.CreateVector(GPRSize);
    Y.CreateVector(GPRSize);
    KS.CreateMatrix(GPRSize,GPRSize);

// main kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos,kblock)
    for(size_t indi=0; indi < GPRSize; indi++){
        for(size_t indj=0; indj < GPRSize; indj++){
            KS[indi][indj] = GetKernelValue(indi,indj);
        }
    }

// error of data points
    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        KS[indi][indi] += Noise;
    }

    ofstream ofs("kernel.orig");
    for(size_t i=0; i < GPRSize; i++){
            ofs << i << " " << setprecision(16) << KS[70][i] << endl;
    }
    ofs.close();

    RunBlasLapackPar();

// inverting the K+Sigma
    int result = 0;
    vout << "   Inverting K+Sigma by SVD (divide and conquer driver) ..." << endl;
    int rank = 0;
    double realRCond = 0;
    result = CSciLapack::invSVD2(KS,logdetK,RCond,rank,realRCond);
    vout << "      Rank = " << rank << "; Info = " << result << "; Real rcond = " << scientific << realRCond << fixed << endl;

    if( result != 0 ){
        RUNTIME_ERROR("unable to invert K+Sigma");
    }

    result = CSciLapack::invSVD2(KS,logdetK,RCond,rank,realRCond);
    vout << "      Rank = " << rank << "; Info = " << result << "; Real rcond = " << scientific << realRCond << fixed << endl;

    ofstream ofs1("kernel.proc");
    for(size_t i=0; i < GPRSize; i++){
            ofs1 << i << " " << setprecision(16) << KS[70][i] << endl;
    }
    ofs1.close();

}

//------------------------------------------------------------------------------

void CGPFilter::TrainProcess(const CVectorDataPtr& in, size_t offset)
{
    if( in->GetLength()-offset < GPRSize ){
        RUNTIME_ERROR("not enough data to process")
    }

    Offset = 0.0;

    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        Y[indi] = in->GetRawDataField()[indi+offset];
    }

    // get average
    for(size_t indi=0; indi < GPRSize; indi++){
        Offset += Y[indi];
    }
    Offset /= GPRSize;

    #pragma omp parallel for
    for(size_t indi=0; indi < GPRSize; indi++){
        Y[indi] -= Offset;
    }

    // calculate weights
    RunBlasLapackPar();
    CSciBlas::gemv(1.0,KS,Y,0.0,GPRModel);
}

//------------------------------------------------------------------------------

void CGPFilter::PredictData(CVectorDataPtr& out, size_t offset)
{
    if( out->GetLength()-offset < GPRSize ){
        RUNTIME_ERROR("not enough data to process")
    }

    RunBlasLapackSeq();

    CSimpleVector<double> kff;
    kff.CreateVector(GPRSize);

    #pragma omp parallel for firstprivate(ipos)
    for(size_t indi=0; indi < GPRSize; indi++){
        GetKFF(kff,indi);
        double value = CSciBlas::dot(kff,GPRModel);
        out->GetRawDataField()[indi+offset] = value + Offset;
    }
}

//------------------------------------------------------------------------------

void CGPFilter::GetKFF(CSimpleVector<double>& kff, double indi)
{
    for(size_t indj=0; indj < GPRSize; indj++){
        kff[indj] = GetKernelValue(indi,indj);
    }
}

//------------------------------------------------------------------------------

void CGPFilter::SetKernel(const CSmallString& kernel)
{
    if( kernel == "ardse" ){
        Kernel = EGPFK_ARDSE;
    } else if( kernel == "ardmc52" ) {
        Kernel = EGPFK_ARDMC52;
    } else if( kernel == "ardmc32" ) {
        Kernel = EGPFK_ARDMC32;
    } else if( kernel == "default" ) {
        Kernel = EGPFK_ARDSE;
    } else {
        CSmallString error;
        error << "Specified kernel '" << kernel << "' is not supported. "
                 "Supported kernels are: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), ardmc32 (ARD Matern class 3/2), default(=ardse)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

const CSmallString CGPFilter::GetKernelName(void)
{
    switch(Kernel){
    case(EGPFK_ARDSE):
        return("ARD squared exponential");
    case(EGPFK_ARDMC52):
        return("ARD Matern class 5/2");
    case(EGPFK_ARDMC32):
        return("ARD Matern class 3/2");
    default:
        RUNTIME_ERROR("not implemented");
    }
}

//------------------------------------------------------------------------------

double CGPFilter::GetKernelValue(double indi,double indj)
{
    double scdist  = fabs((indi-indj)/WFac);
    double scdist2 = scdist * scdist;

    // get kernel value
    switch(Kernel){
        case(EGPFK_ARDSE): {
            return(SigmaF2*exp(-0.5*scdist2));
        }
        break;
        case(EGPFK_ARDMC52):{
            return(SigmaF2*(1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2)*exp(-sqrt(5.0)*scdist));
        }
        break;
        case(EGPFK_ARDMC32):{
            return(SigmaF2*(1.0+sqrt(3.0)*scdist)*exp(-sqrt(3.0)*scdist));
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

double CGPFilter::GetLogML(void)
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

double CGPFilter::GetLogPL(void)
{
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

void CGPFilter::PrintExecInfo(CVerboseStr& vout)
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

void CGPFilter::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CGPFilter::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
