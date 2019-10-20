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

    SplitNCorr          = false;
    SigmaF2             = 15.0;

    IncludeError        = false;
    IncludeGluedBins    = false;
    NoEnergy            = false;
    GlobalMinSet        = false;

    UseAnalyticalK      = true;
    Method              = EGPRINV_LU;
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

    if( (int)sncorrs.size() > Accumulator->GetNumberOfCoords() ){
        CSmallString error;
        error << "too many ncorrs (" << sncorrs.size() << ") than required (" << Accumulator->GetNumberOfCoords() << ")";
        RUNTIME_ERROR(error);
    }

    if( (SplitNCorr == false) && (sncorrs.size() > 1) ){
        CSmallString error;
        error << "too many ncorrs (" << sncorrs.size() << ") but only one allowed without --splitncorr";
        RUNTIME_ERROR(error);
    }

    NCorr.CreateVector(Accumulator->GetNumberOfCoords());

    // parse values of ncorr
    double last_ncorr = 1.0;
    for(int i=0; i < (int)sncorrs.size(); i++){
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
    for(int i=sncorrs.size(); i < Accumulator->GetNumberOfCoords(); i++){
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
    if( ncorr.GetLength() != Accumulator->GetNumberOfCoords() ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    NCorr = ncorr;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetNCorr(int cvind, double value)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetNCorr");
    }
    if( (cvind < 0) || (cvind >= Accumulator->GetNumberOfCoords()) ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is NCorr initialized?
    if( NCorr.GetLength() <= 0 ){
        NCorr.CreateVector(Accumulator->GetNumberOfCoords());
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

    if( (int)swfacs.size() > Accumulator->GetNumberOfCoords() ){
        CSmallString error;
        error << "too many wfacs (" << swfacs.size() << ") than required (" << Accumulator->GetNumberOfCoords() << ")";
        RUNTIME_ERROR(error);
    }

    WFac.CreateVector(Accumulator->GetNumberOfCoords());

    // parse values of wfac
    double last_wfac = 3.0;
    for(int i=0; i < (int)swfacs.size(); i++){
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
    for(int i=swfacs.size(); i < Accumulator->GetNumberOfCoords(); i++){
        WFac[i] = last_wfac;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetWFac(CSimpleVector<double>& wfac)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }
    if( wfac.GetLength() != Accumulator->GetNumberOfCoords() ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    WFac = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetWFac(int cvind, double value)
{
    if( Accumulator == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }
    if( (cvind < 0) || (cvind >= Accumulator->GetNumberOfCoords()) ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is wfac initialized?
    if( WFac.GetLength() <= 0 ){
        WFac.CreateVector(Accumulator->GetNumberOfCoords());
    }
    WFac[cvind] = value;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetIncludeError(bool set)
{
    IncludeError = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetNumericK(bool set)
{
    UseAnalyticalK = !set;
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

void CABFIntegratorGPR::SetINVMehod(EGPRINVMethod set)
{
    Method = set;
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
    GPos.CreateVector(Accumulator->GetNumberOfCoords());
    stringstream str(sspec);
    for(int i=0; i < Accumulator->GetNumberOfCoords(); i++){
        str >> GPos[i];
        if( ! str ){
            CSmallString error;
            error << "unable to decode CV value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorGPR::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    InitHyprmsOpt();

    int ind = 0;
    int soffset = 0;
    int noffset;
    int woffset;

    if( SplitNCorr ){
        noffset = soffset+1;
        woffset = noffset+NCVs;
    } else {
        noffset = soffset+1;
        woffset = noffset+1;
    }

    for(int prm=0; prm < (int)flags.size(); prm++){
        // shall we calc der?
        if( flags[prm] == false ) {
            continue;
        }

        // calc Kder
        if( prm == 0 ){
            // sigmaf2
            CalcKderWRTSigmaF2();
        } else if( (prm >= noffset) && (prm < woffset) ){
            int cv = prm - noffset;
            CalcKderWRTNCorr(cv);
        } else {
            int cv = prm - woffset;
            CalcKderWRTWFac(cv);
        }

        // calc trace
        double tr = 0.0;
        #pragma omp parallel for reduction(+:tr)
        for(int i=0; i < GPRSize; i++){
            for(int j=0; j < GPRSize; j++){
                tr += (ATA[i][j]-K[i][j])*Kder[j][i];
            }
        }

        // finalize derivative
        der[ind] = 0.5*tr;
        ind++;
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::InitHyprmsOpt(void)
{
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

    ATA.CreateMatrix(GPRSize,GPRSize);
    Kder.CreateMatrix(GPRSize,GPRSize);

    // calc ATA matrix
    #pragma omp parallel for
    for(int i=0; i < GPRSize; i++){
        for(int j=0; j < GPRSize; j++){
            ATA[i][j] = GPRModel[i]*GPRModel[j];
        }
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

    #pragma omp parallel for firstprivate(ipos,jpos)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(int indj=0; indj < NumOfUsedBins; indj++){
            int j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            double arg = 0.0;
            for(int ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ipos[ii],jpos[ii]);
                double dd = CVLengths2[ii];
                arg += du*du/(2.0*dd);
            }
            arg = exp(-arg);

            // calculate block of second derivatives
            for(int ii=0; ii < NCVs; ii++){
                for(int jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ipos[ii],jpos[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ipos[jj],jpos[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    Kder[indi*NCVs+ii][indj*NCVs+jj] -= arg*du/dd;
                    if( (ii == jj) ){
                        Kder[indi*NCVs+ii][indj*NCVs+jj] += arg/CVLengths2[ii];
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CalcKderWRTNCorr(int cv)
{
    Kder.SetZero();

    #pragma omp parallel for
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        for(int ii=0; ii < NCVs; ii++){
            // get mean force and its error
            double er = Accumulator->GetValue(ii,i,EABF_MEAN_FORCE_ERROR);
            if( SplitNCorr ){
                if( cv == ii ) Kder[indi*NCVs+ii][indi*NCVs+ii] += er*er;
            } else {
                Kder[indi*NCVs+ii][indi*NCVs+ii] += er*er;
            }
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CalcKderWRTWFac(int cv)
{
    Kder.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    double wf = WFac[cv];
    double wd3 = 1.0/(CVLengths2[cv]*wf);
    double wd5 = wd3/CVLengths2[cv];

    // create kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(int indj=0; indj < NumOfUsedBins; indj++){
            int j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            double dc = Accumulator->GetCoordinate(cv)->GetDifference(ipos[cv],jpos[cv]);

            double arg = 0.0;
            for(int ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ipos[ii],jpos[ii]);
                double dd = CVLengths2[ii];
                arg += du*du/(2.0*dd);
            }
            arg = SigmaF2*exp(-arg);

            double argd = arg*dc*dc*wd3;

            // calculate block of second derivatives

            for(int ii=0; ii < NCVs; ii++){
                for(int jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ipos[ii],jpos[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ipos[jj],jpos[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    Kder[indi*NCVs+ii][indj*NCVs+jj] -= argd*du/dd;
                    if( (cv == ii) && (cv != jj) ) {
                        Kder[indi*NCVs+ii][indj*NCVs+jj] -= -2.0*arg*du*wd3/(CVLengths2[jj]);
                    } else if( (cv == jj) && (cv != ii) ) {
                        Kder[indi*NCVs+ii][indj*NCVs+jj] -= -2.0*arg*du*wd3/(CVLengths2[ii]);
                    } else if( (cv == ii) && (cv == jj) ) {
                        Kder[indi*NCVs+ii][indj*NCVs+jj] -= -4.0*arg*du*wd5;
                    }
                    if( (ii == jj) ){
                        Kder[indi*NCVs+ii][indj*NCVs+ii] += argd/CVLengths2[ii];
                        if( ii == cv ){
                            Kder[indi*NCVs+ii][indj*NCVs+ii] += -2.0*arg*wd3;
                        }
                    }
                }
            }
        }
    }
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
    CVLengths2.CreateVector(Accumulator->GetNumberOfCoords());
    for(int i=0; i < Accumulator->GetNumberOfCoords(); i++){
        double l = WFac[i]*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        CVLengths2[i] = l*l;
    }

    // number of data points
    NumOfUsedBins = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) > 0 ) NumOfUsedBins++;
    }
    NCVs = Accumulator->GetNumberOfCoords();
    GPRSize = NumOfUsedBins*NCVs;

    // create sampled map
    SampledMap.CreateVector(NumOfUsedBins);
    int ind = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;
        SampledMap[ind] = i;
        ind++;
    }

    // init GPR weights
    GPRModel.CreateVector(GPRSize);

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
        vout << format("      SigmaF2  = %10.4f")%SigmaF2 << endl;
    if( SplitNCorr ){
        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++ ){
            vout << format("      NCorr#%-2d = %10.4f")%(k+1)%NCorr[k] << endl;
        }
    } else {
        vout << format("      NCorr    = %10.4f")%NCorr[0] << endl;
    }

    for(int k=0; k < Accumulator->GetNumberOfCoords(); k++ ){
        vout << format("      WFac#%-2d  = %10.4f")%(k+1)%WFac[k] << endl;
    }

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    if( ! nostat ){
        // and finaly some statistics
        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++ ){
        vout << "   RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
        }

        // and marginal likelihood
        vout << "   logML = " << setprecision(5) << GetLogMarginalLikelihood() << endl;
    }

    // finalize FES if requested
    if( ! NoEnergy ){
        CalculateEnergy(vout);
        if( IncludeError ){
            CalculateErrors(GPos,vout);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::PrintExecInfo(CVerboseStr& vout)
{
#if defined(_OPENMP)
    {
        int ncpus = omp_get_max_threads();
        vout << "   OpenMP - number of threads: " << ncpus << endl;
    }
#else
    vout << "   No OpenMP - sequential mode." << endl;
#endif
    CSciLapack::PrintExecInfo(vout);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::TrainGP(CVerboseStr& vout)
{
    if( UseAnalyticalK ) {
        vout << "   Creating K+Sigma and Y ..." << endl;
    } else {
        vout << "   Creating K+Sigma and Y (numeric differentation) ..." << endl;
    }

    vout << "   Dim: " << GPRSize << " x " << GPRSize << endl;

    K.CreateMatrix(GPRSize,GPRSize);
    Y.CreateVector(GPRSize);

// construct K
    if( UseAnalyticalK ){
        AnalyticalK();
    } else {
        NumericalK();
    }

//// debug - print K+sigma matrix
//    for(int i=0; i < GPRSize; i++){
//        for(int j=0; j < GPRSize; j++){
//            vout << format("%12.5e ")%K[i][j];
//        }
//        vout << endl;
//    }

// construct Y
    #pragma omp parallel for
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];
        for(int ii=0; ii < NCVs; ii++){
            double mf = Accumulator->GetValue(ii,i,EABF_MEAN_FORCE_VALUE);
            Y[indi*NCVs+ii] = mf;
        }
    }



// inverting the K+Sigma
    int result = 0;
    switch(Method){
        case(EGPRINV_LU):
            vout << "   Inverting K+Sigma by LU ..." << endl;
            result = CSciLapack::invLU(K,logdetK);
            if( result != 0 ) return(false);
            break;
        case(EGPRINV_LL):
            vout << "   Inverting K+Sigma by LL ..." << endl;
            result = CSciLapack::invLL(K,logdetK);
            if( result != 0 ) return(false);
            break;
        case(EGPRINV_SVD):{
            vout << "   Inverting K+Sigma by SVD (divide and conquer driver) ..." << endl;
            int rank = 0;
            result = CSciLapack::invSVD2(K,logdetK,RCond,rank);
            vout << "   Rank = " << rank << "; Info = " << result << endl;
            if( result != 0 ) return(false);
            }
            break;
        case(EGPRINV_SVD2):{
            vout << "   Inverting K+Sigma by SVD2 (simple driver) ..." << endl;
            int rank = 0;
            result = CSciLapack::invSVD1(K,logdetK,RCond,rank);
            vout << "   Rank = " << rank << "; Info = " << result << endl;
            if( result != 0 ) return(false);
            }
            break;
    default:
        INVALID_ARGUMENT("unsupported method");
    }

// calculate weights
    vout << "   Calculating weights B ..." << endl;
    CSciBlas::gemv(1.0,K,Y,0.0,GPRModel);

    return( true );
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::AnalyticalK(void)
{
    K.SetZero();

    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    // main kernel matrix
    #pragma omp parallel for firstprivate(ipos,jpos)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(int indj=0; indj < NumOfUsedBins; indj++){
            int j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            double arg = 0.0;
            for(int ii=0; ii < NCVs; ii++){
                double du = Accumulator->GetCoordinate(ii)->GetDifference(ipos[ii],jpos[ii]);
                double dd = CVLengths2[ii];
                arg += du*du/(2.0*dd);
            }
            arg = SigmaF2*exp(-arg);

            // calculate block of second derivatives

            for(int ii=0; ii < NCVs; ii++){
                for(int jj=0; jj < NCVs; jj++){
                    double du = Accumulator->GetCoordinate(ii)->GetDifference(ipos[ii],jpos[ii]) *
                                Accumulator->GetCoordinate(jj)->GetDifference(ipos[jj],jpos[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    K[indi*NCVs+ii][indj*NCVs+jj] -= arg*du/dd;
                    if( (ii == jj) ){
                        K[indi*NCVs+ii][indj*NCVs+jj] += arg/CVLengths2[ii];
                    }
                }
            }
        }
    }

// error of data points
    #pragma omp parallel for
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];
        for(int ii=0; ii < NCVs; ii++){
            double er = Accumulator->GetValue(ii,i,EABF_MEAN_FORCE_ERROR);
            if( SplitNCorr ){
                K[indi*NCVs+ii][indi*NCVs+ii] += er*er*NCorr[ii];
            } else {
                K[indi*NCVs+ii][indi*NCVs+ii] += er*er*NCorr[0];
            }
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::NumericalK(void)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> jpos;
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);

    CFortranMatrix kblock;
    kblock.CreateMatrix(NCVs,NCVs);

    // main kernel matrix
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        for(int indj=0; indj < NumOfUsedBins; indj++){
            int j = SampledMap[indj];

            Accumulator->GetPoint(j,jpos);

            // calc Kblock
            NumericalKBlock(ipos,jpos,kblock);

            // distribute to main kernel matrix
            for(int ii=0; ii < NCVs; ii++){
                for(int jj=0; jj < NCVs; jj++){
                    K[indi*NCVs+ii][indj*NCVs+jj] = kblock[ii][jj];
                }
            }
        }
    }

// error of data points
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];
        for(int ii=0; ii < NCVs; ii++){
            double er = Accumulator->GetValue(ii,i,EABF_MEAN_FORCE_ERROR);
            if( SplitNCorr ){
                K[indi*NCVs+ii][indi*NCVs+ii] += er*er*NCorr[ii];
            } else {
                K[indi*NCVs+ii][indi*NCVs+ii] += er*er*NCorr[0];
            }
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::NumericalKBlock(CSimpleVector<double>& ip,CSimpleVector<double>& jp,CFortranMatrix& kblok)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NCVs);
    CSimpleVector<double> tjp;
    tjp.CreateVector(NCVs);

    double  dh = 1e-5;
    double  v1,v2,v3,v4;

    kblok.SetZero();

    // off diagonal elements
    for(int ii=0; ii < NCVs; ii++) {
        for(int jj=0; jj < NCVs; jj++) {

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

            kblok[ii][jj] = (v1 - v2 - v3 + v4)/(4.0*dh*dh);
        }
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetKernelValue(CSimpleVector<double>& ip,CSimpleVector<double>& jp)
{
    double arg = 0.0;
    for(int ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        arg += du*du/(2.0*dd);
    }
    return(SigmaF2*exp(-arg));
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorGPR::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating FES ..." << endl;

// create map for bins with calculated energy and error
    NumOfValues = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        if( IncludeGluedBins ){
            if( samples == 0 ) continue;
        } else {
            if( samples <= 0 ) continue;
        }
        NumOfValues++;
    }
    ValueMap.CreateVector(NumOfValues);

    int indi = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
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
    for(int indj=0; indj < NumOfValues; indj++){
        int j = ValueMap[indj];
        Accumulator->GetPoint(j,jpos);
        values[indj] = GetValue(jpos);
    }

// basic FES update
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        FES->SetNumOfSamples(i,samples);
        FES->SetEnergy(i,0.0);
        FES->SetError(i,0.0);
    }

// update FES
    if( GlobalMinSet ){
        // GPos.CreateVector(NCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating FES ..." << endl;
        vout << "       Global minimum provided at: ";
        vout << GPos[0];
        for(int i=1; i < Accumulator->GetNumberOfCoords(); i++){
            vout << "x" << GPos[i];
        }
        double glb_min = GetValue(GPos);
        vout << " (" << glb_min << ")" << endl;

        for(int indj=0; indj < NumOfValues; indj++){
            int j = ValueMap[indj];
            FES->SetEnergy(j,values[indj]-glb_min);
        }
    } else {
        // search for global minimum
        GPos.CreateVector(NCVs);
        bool   first = true;
        double glb_min = 0.0;
        for(int indj=0; indj < NumOfValues; indj++){
            int j = ValueMap[indj];
            double value = values[indj];
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                Accumulator->GetPoint(j,GPos);
            }
        }

   //   vout << "   Calculating FES ..." << endl;
        vout << "       Global minimum found at: ";
        vout << GPos[0];
        for(int i=1; i < Accumulator->GetNumberOfCoords(); i++){
            vout << "x" << GPos[i];
        }
        vout << " (" << glb_min << ")" << endl;

        for(int indj=0; indj < NumOfValues; indj++){
            int j = ValueMap[indj];
            FES->SetEnergy(j,values[indj]-glb_min);
        }
    }

    vout << "   SigmaF2 = " << setprecision(5) << FES->GetSigmaF2() << endl;
    if( IncludeGluedBins ){
        vout << "   SigmaF2 (including glued bins) = " << setprecision(5) << FES->GetSigmaF2(true) << endl;
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetValue(const CSimpleVector<double>& position)
{
    double energy = 0.0;

    CSimpleVector<double> ipos;
    ipos.CreateVector(NCVs);

    // parallel execution was moved one level up
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        double arg = 0.0;
        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(position[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            arg += du*du/(2.0*dd);
        }
        arg = SigmaF2*exp(-arg);

        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(position[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            energy += GPRModel[indi*NCVs+ii]*(du/dd)*arg;
        }
    }

    return(energy);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetMeanForce(const CSimpleVector<double>& position,int icoord)
{
    double mf = 0.0;

    CSimpleVector<double> ipos;
    ipos.CreateVector(NCVs);

    // parallel execution was moved one level up
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        double arg = 0.0;
        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(position[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            arg += du*du/(2.0*dd);
        }
        arg = SigmaF2*exp(-arg);

        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(position[ii],ipos[ii])*
                        Accumulator->GetCoordinate(icoord)->GetDifference(position[icoord],ipos[icoord]);
            double dd = CVLengths2[ii]*CVLengths2[icoord];
            double der2 = -arg*du/dd;
            if( ii == icoord ){
                der2 += arg/CVLengths2[ii];
            }
            mf += GPRModel[indi*NCVs+ii]*der2;
        }
    }

    return(mf);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetRMSR(int cv)
{
    if( Accumulator->GetNumberOfBins() <= 0 ){
        ES_ERROR("number of bins is not > 0");
        return(0.0);
    }

    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    double rmsr = 0.0;

    #pragma omp parallel for firstprivate(jpos) reduction(+:rmsr)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,jpos);

        double mfi = Accumulator->GetValue(cv,i,EABF_MEAN_FORCE_VALUE);
        double mfp = GetMeanForce(jpos,cv);
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
    if( Accumulator->GetNumberOfBins() <= 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> jpos;
    CSimpleVector<double> mfi;
    CSimpleVector<double> mfp;

    jpos.CreateVector(NCVs);
    mfi.CreateVector(GPRSize);
    mfp.CreateVector(GPRSize);

    // calculate
    #pragma omp parallel for firstprivate(jpos)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];
        Accumulator->GetPoint(i,jpos);
        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            mfi[indi*NCVs+k] = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE);
            mfp[indi*NCVs+k] = GetMeanForce(jpos,k);
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
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,jpos);

        for(int c=0; c < Accumulator->GetNumberOfCoords(); c++){
            ofs << format("%20.16f ")%jpos[c];
        }

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            ofs << format(" %20.16f %20.16f")%mfi[indi*NCVs+k]%mfp[indi*NCVs+k];
        }

        ofs << endl;
    }


    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::FilterByMFZScore(double zscore,CVerboseStr& vout)
{
    if( Accumulator->GetNumberOfBins() <= 0 ){
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

    flags.CreateVector(Accumulator->GetNumberOfBins());
    sig2.CreateVector(NCVs);
    maxi.CreateVector(NCVs);
    maxzscore.CreateVector(NCVs);
    mferror2.CreateVector(GPRSize);

    flags.Set(1);

    vout << high;
    vout << "   Precalculating MF errors ..." << endl;

    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    // precalculate values
    #pragma omp parallel for firstprivate(jpos)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,jpos);

        for(int k=0; k < NCVs; k++){
            double diff2 = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(jpos,k);
            diff2 *= diff2;
            mferror2[indi*NCVs+k] = diff2;
        }
    }

    vout << "   Searching for MF outliers ..." << endl;
    vout << debug;

    bool testme = true;
    while( testme ){

        testme = false;

        double  count = 0;
        sig2.SetZero();

        // calc variances - we assume zero mean on errors
        for(int indi=0; indi < NumOfUsedBins; indi++){
            int i = SampledMap[indi];

            if( flags[i] != 0 ) {
                for(int k=0; k < NCVs; k++){
                    sig2[k] += mferror2[indi*NCVs+k];
                }
                count++;
            }
        }

        if( count == 0 ) break;
        for(int k=0; k < NCVs; k++){
            sig2[k] /= count;
        }

        // filter
        bool first = true;
        for(int indi=0; indi < NumOfUsedBins; indi++){
            int i = SampledMap[indi];

            if( flags[i] != 0 ){
                Accumulator->GetPoint(i,jpos);

                for(int k=0; k < NCVs; k++){
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

        for(int k=0; k < NCVs; k++){
            if( maxzscore[k] > zscore ){
                flags[maxi[k]] = 0;
                testme = true;
                vout << "   outlier found at " << maxi[k] << " for cv " << k << " with z-score " << sqrt(maxzscore[k]) << endl;
            }
        }
    }

    // apply limits
    int outliers = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( flags[i] == 0 ){
            Accumulator->SetNumberOfABFSamples(i,0);
            outliers++;
        }
    }

    // from now the integrator is in invalid state !!!
    vout << high;
    vout << "   Number of outliers = " << outliers << endl;
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetLogMarginalLikelihood(void)
{
    double ml = 0.0;

    // http://www.gaussianprocess.org/gpml/chapters/RW5.pdf
    // page 113

    ml -= CSciBlas::dot(Y,GPRModel);
    ml -= logdetK;
    ml -= GPRSize * log(2*M_PI);
    ml *= 0.5;

    return(ml);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetVar(CSimpleVector<double>& lpos)
{
    CSimpleVector<double> lk;
    CSimpleVector<double> ik;
    CSimpleVector<double> ipos;

    lk.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);
    ipos.CreateVector(NCVs);

    #pragma omp parallel for firstprivate(ipos)
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        double argl = 0.0;
        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            argl += du*du/(2.0*dd);
        }
        argl = SigmaF2*exp(-argl);

        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            lk[indi*NCVs + ii] = (du/dd)*argl;
        }
    }

    double cov = SigmaF2;

    CSciBlas::gemv(1.0,K,lk,0.0,ik);
    cov = cov  - CSciBlas::dot(lk,ik);

    return(cov);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& lrcov,double& rrvar)
{
    CSimpleVector<double> ipos;
    CSimpleVector<double> rk;
    CSimpleVector<double> lk;
    CSimpleVector<double> ik;

    ipos.CreateVector(NCVs);
    rk.CreateVector(GPRSize);
    lk.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);
    ik.SetZero();

    // parallel execution was moved one level up
    for(int indi=0; indi < NumOfUsedBins; indi++){
        int i = SampledMap[indi];

        Accumulator->GetPoint(i,ipos);

        double argl = 0.0;
        double argr = 0.0;
        double du,dd;
        for(int ii=0; ii < NCVs; ii++){
            dd = CVLengths2[ii];

            du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],ipos[ii]);
            argl += du*du/(2.0*dd);

            du = Accumulator->GetCoordinate(ii)->GetDifference(rpos[ii],ipos[ii]);
            argr += du*du/(2.0*dd);
        }
        argl = SigmaF2*exp(-argl);
        argr = SigmaF2*exp(-argr);

        for(int ii=0; ii < NCVs; ii++){
            dd = CVLengths2[ii];

            du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],ipos[ii]);
            lk[indi*NCVs + ii] = (du/dd)*argl;

            du = Accumulator->GetCoordinate(ii)->GetDifference(rpos[ii],ipos[ii]);
            rk[indi*NCVs + ii] = (du/dd)*argr;
        }
    }

// do not run threaded version of Lapack/Blas
    CSciLapack::SetNumThreadsLocal(1);

// threading is dissabled above
    CSciBlas::gemv(1.0,K,rk,0.0,ik);

// covariance
    double arg = 0.0;
    for(int ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],rpos[ii]);
        double dd = CVLengths2[ii];
        arg += du*du/(2.0*dd);
    }

// threading is dissabled above
    lrcov = SigmaF2*exp(-arg) - CSciBlas::dot(lk,ik);
    rrvar = SigmaF2           - CSciBlas::dot(rk,ik);
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
    for(int indj=0; indj < NumOfValues; indj++){
        int j = ValueMap[indj];
        Accumulator->GetPoint(j,jpos);

        double varfc,covfg;
        GetCovVar(gpos,jpos,covfg,varfc);

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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
