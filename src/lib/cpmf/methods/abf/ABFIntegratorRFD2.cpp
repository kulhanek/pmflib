// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//                       Martin Petrek, petrek@chemi.muni.cz
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

#include <ABFIntegratorRFD2.hpp>
#include <ABFAccumulator.hpp>
#include <EnergySurface.hpp>
#include <ErrorSystem.hpp>
#include <SciLapack.hpp>
#include <SciBlas.hpp>
#include <iomanip>

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegratorRFD2::CABFIntegratorRFD2(void)
{
    Accumulator = NULL;
    FES = NULL;

    NumOfVariables = 0;
    NumOfEquations = 0;

    Periodicity = false;
    FDLevel = 4;

    Method = ERFDLLS_SVD;

    RCond   = -1; // machine precision
}

//------------------------------------------------------------------------------

CABFIntegratorRFD2::~CABFIntegratorRFD2(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorRFD2::SetInputABFAccumulator(const CABFAccumulator* p_accu)
{
    Accumulator = p_accu;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD2::SetOutputFESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD2::SetFDPoints(int npts)
{
    FDLevel = npts;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD2::SetPeriodicity(bool set)
{
    Periodicity = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD2::SetLLSMehod(ERFDLLSMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD2::SetRCond(double rcond)
{
    RCond = rcond;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRFD2::Integrate(CVerboseStr& vout)
{
    if( Accumulator == NULL ) {
        ES_ERROR("ABF accumulator is not set");
        return(false);
    }
    if( FES == NULL ) {
        ES_ERROR("FES is not set");
        return(false);
    }

    if( Accumulator->GetNumberOfCoords() == 0 ) {
        ES_ERROR("number of coordinates is zero");
        return(false);
    }

    if( Accumulator->GetNumberOfCoords() != FES->GetNumberOfCoords() ){
        ES_ERROR("inconsistent ABF and FES - CVs");
        return(false);
    }
    if( Accumulator->GetNumberOfBins() != FES->GetNumberOfPoints() ){
        ES_ERROR("inconsistent ABF and FES - points");
        return(false);
    }

    if( BuildSystemOfEquations(vout) == false ) {
        return(false);
    }
    // make copy of A
    BA = A;
    BRhs = Rhs;
    if( SolveSystemOfEquations(vout) == false ) {
        return(false);
    }

    // release A and Rhs
    A.FreeMatrix();
    Rhs.FreeVector();

    // and finaly some statistics
    for(int k=0; k < Accumulator->GetNumberOfCoords(); k++ ){
    vout << "   RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
    }

// find global minimum
    double glb_min = X[0];
    for(size_t i=0; i < X.GetLength(); i++) {
        if(glb_min > X[i]) glb_min = X[i];
    }

// load data to FES
    for(int ipoint=0; ipoint < FES->GetNumberOfPoints(); ipoint++) {
        int x_index = XMap[ipoint];
        if(x_index >= 0) {
            double value = X[x_index]-glb_min;
            FES->SetEnergy(ipoint,value);
            FES->SetNumOfSamples(ipoint,Accumulator->GetNumberOfABFSamples(ipoint));
        }
    }

    vout << "   SigmaF2 = " << setprecision(5) << FES->GetSigmaF2() << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRFD2::BuildSystemOfEquations(CVerboseStr& vout)
{
    vout << "   Creating A and rhs ..." << endl;

// initializations
    IPoint.CreateVector(Accumulator->GetNumberOfCoords());
    XMap.CreateVector(Accumulator->GetNumberOfBins());

    XMap.Set(-1);

    NumOfVariables = 0;
    NumOfEquations = 0;

    BuildEquations(true);

    vout << "   Dim: " << NumOfEquations << " x " << NumOfVariables << endl;

    if(NumOfVariables <= 0) {
        CSmallString error;
        error << "number of function nodes is zero - system cannot be integrated";
        ES_ERROR(error);
        return(false);
    }

    if(NumOfEquations < NumOfVariables) {
        CSmallString error;
        error << "incomplete system - number of equations (" << CSmallString(NumOfEquations) << ") is smaller than number of variables (" << CSmallString(NumOfVariables) <<  ")";
        ES_ERROR(error);
        return(false);
    }

// allocate A and rhs
    if( A.CreateMatrix(NumOfEquations,NumOfVariables) == false ){
        CSmallString error;
        error << "unable to allocate A matrix (M: " << CSmallString(NumOfEquations) << ", N: " << CSmallString(NumOfVariables) << ")";
        ES_ERROR(error);
        return(false);
    }

    Rhs.CreateVector(NumOfEquations);
    RhsCv.CreateVector(NumOfEquations);

// build system of equations
    BuildEquations(false);

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD2::BuildEquations(bool trial)
{
    int LocIter = 0;

    if( ! trial ){
        A.SetZero();
    }

    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        // number of samples is controlled via GetFBinIndex

        Accumulator->GetIPoint(i,IPoint);

        for(int ifcoord=0; ifcoord < Accumulator->GetNumberOfCoords(); ifcoord++) {
            int ifbin1,ifbin2,ifbin3,ifbin4;

            ifbin1 = GetFBinIndex(IPoint,ifcoord,0);
            ifbin2 = GetFBinIndex(IPoint,ifcoord,1);
            ifbin3 = GetFBinIndex(IPoint,ifcoord,2);
            ifbin4 = GetFBinIndex(IPoint,ifcoord,3);

            const CColVariable* p_coord = Accumulator->GetCoordinate(ifcoord);
            double              diff = p_coord->GetBinWidth();

            switch(FDLevel) {
                case 3:
                    if((ifbin1 == -1) || (ifbin2 == -1) || (ifbin3 == -1)) break;

                    // map nodes to XMap -------------
                    if(XMap[ifbin1] < 0) {
                        XMap[ifbin1] = NumOfVariables;
                        NumOfVariables++;
                    }

                    if(XMap[ifbin2] < 0) {
                        XMap[ifbin2] = NumOfVariables;
                        NumOfVariables++;
                    }

                    if(XMap[ifbin3] < 0) {
                        XMap[ifbin3] = NumOfVariables;
                        NumOfVariables++;
                    }

                    // set A elements ----------------
                    if(! trial) {
                        double dfac = 1.0/(diff * 2.0);
                        A[LocIter][XMap[ifbin1]] = -3.0 * dfac;
                        A[LocIter][XMap[ifbin2]] = +4.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = -1.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin1,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                        A[LocIter][XMap[ifbin1]] = -1.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = +1.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin2,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                        A[LocIter][XMap[ifbin1]] = +1.0 * dfac;
                        A[LocIter][XMap[ifbin2]] = -4.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = +3.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin3,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                    } else {
                        NumOfEquations += 3;
                    }
                    break;
                case 4:
                    if((ifbin1 == -1) || (ifbin2 == -1) || (ifbin3 == -1) || (ifbin4 == -1)) break;

                    if(XMap[ifbin1] < 0) {
                        XMap[ifbin1] = NumOfVariables;
                        NumOfVariables++;
                    }

                    if(XMap[ifbin2] < 0) {
                        XMap[ifbin2] = NumOfVariables;
                        NumOfVariables++;
                    }

                    if(XMap[ifbin3] < 0) {
                        XMap[ifbin3] = NumOfVariables;
                        NumOfVariables++;
                    }

                    if(XMap[ifbin4] < 0) {
                        XMap[ifbin4] = NumOfVariables;
                        NumOfVariables++;
                    }

                    if(! trial) {
                        double dfac = 1.0/(diff * 6.0);
                        A[LocIter][XMap[ifbin1]] = -11.0 * dfac;
                        A[LocIter][XMap[ifbin2]] = +18.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = -9.0 * dfac;
                        A[LocIter][XMap[ifbin4]] = +2.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin1,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                        A[LocIter][XMap[ifbin1]] = -2.0 * dfac;
                        A[LocIter][XMap[ifbin2]] = -3.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = +6.0 * dfac;
                        A[LocIter][XMap[ifbin4]] = -1.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin2,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                        A[LocIter][XMap[ifbin1]] = +1.0 * dfac;
                        A[LocIter][XMap[ifbin2]] = -6.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = +3.0 * dfac;
                        A[LocIter][XMap[ifbin4]] = +2.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin3,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                        A[LocIter][XMap[ifbin1]] = -2.0 * dfac;
                        A[LocIter][XMap[ifbin2]] = +9.0 * dfac;
                        A[LocIter][XMap[ifbin3]] = -18.0 * dfac;
                        A[LocIter][XMap[ifbin4]] = +11.0 * dfac;
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin4,EABF_MEAN_FORCE_VALUE);
                        LocIter++;
                    } else {
                        NumOfEquations += 4;
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

//------------------------------------------------------------------------------

int CABFIntegratorRFD2::GetFBinIndex(const CSimpleVector<int>& position,int ifcoord,int offset) const
{
    int glbindex = 0;
    for(int i=0; i < Accumulator->GetNumberOfCoords(); i++) {
        const CColVariable* p_coord = Accumulator->GetCoordinate(i);
        int nbins = p_coord->GetNumberOfBins();
        int pos = position[i];
        if(i == ifcoord) {
            pos += offset;
            if(Periodicity == true) {
                if((p_coord->IsPeriodic() == true) && (pos < 0)) pos = nbins + pos;
                if((p_coord->IsPeriodic() == true) && (pos >= nbins)) pos = pos - nbins;
            }
            if(pos < 0) return(-1);
            if(pos >= nbins) return(-1);
        }
        glbindex = glbindex*nbins + pos;
    }

// check if we have sufficient number of samples
    if(Accumulator->GetNumberOfABFSamples(glbindex) <= 0) return(-1);

    return(glbindex);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRFD2::SolveSystemOfEquations(CVerboseStr& vout)
{
    int result = 0;

    switch(Method){
        case ERFDLLS_SVD:{
            vout << "   Solving least square problem by SVD ..." << endl;

            // solve least square problem via GELSD
            int rank = 0;
            double realRCond = 0.0;
            result = CSciLapack::gelsd(A,Rhs,RCond,rank,realRCond);
            vout << "      Rank = " << rank << "; Info = " << result << "; Real rcond = " << scientific << realRCond << fixed << endl;
            if( result != 0 ) return(false);
            }
        break;
        case ERFDLLS_QR:{
            vout << "   Solving least square problem by QR ..." << endl;

            // solve least square problem via GELS
            result = CSciLapack::gels(A,Rhs);
            if( result != 0 ) return(false);
            }
        break;
        default:
            INVALID_ARGUMENT("unsupported method");
    }

    X.CreateVector(NumOfVariables);

    // copy results to X
    for(int l=0; l < NumOfVariables; l++){
        X[l] = Rhs[l];
    }

    return(true);
}

//------------------------------------------------------------------------------

double CABFIntegratorRFD2::GetRMSR(int k)
{
    double rmsr = 0.0;
    double count = 0.0;

    Lhs.CreateVector(NumOfEquations);
    Lhs.SetZero();

    // multiply BA (copy of original A) and X and compare with BRhs (copy of Rhs)
    CSciBlas::gemv(1.0,BA,X,0.0,Lhs);

    // calculate error
    for(int i=0; i < NumOfEquations; i++){
        if( RhsCv[i] == k ){
            double err = Lhs[i] - BRhs[i];
            rmsr = rmsr + err*err;
            count++;
        }
    }

    if( count == 0 ) return(0.0);
    rmsr = rmsr / count;
    if( rmsr > 0 ){
        rmsr = sqrt(rmsr);
    }

    return(rmsr);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
