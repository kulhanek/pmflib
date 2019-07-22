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

#include <ABFIntegratorRFD.hpp>
#include <ABFAccumulator.hpp>
#include <EnergySurface.hpp>
#include <ErrorSystem.hpp>
#include <iomanip>

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegratorRFD::CABFIntegratorRFD(void)
{
    Accumulator = NULL;
    FES = NULL;

    NumOfVariables = 0;
    NumOfEquations = 0;
    NumOfNonZeros = 0;

    A = NULL;
    cA = NULL;
    LocIter = 0;

    Periodicity = false;
    FDLevel = 4;

    IntegratedRealm = EABF_MEAN_FORCE_VALUE;
}

//------------------------------------------------------------------------------

CABFIntegratorRFD::~CABFIntegratorRFD(void)
{
    ReleaseAllResources();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorRFD::SetInputABFAccumulator(const CABFAccumulator* p_accu)
{
    Accumulator = p_accu;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD::SetOutputFESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD::SetFDPoints(int npts)
{
    FDLevel = npts;
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD::SetPeriodicity(bool set)
{
    Periodicity = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRFD::Integrate(CVerboseStr& vout, bool errors)
{
    if( ! errors ){
        IntegratedRealm = EABF_MEAN_FORCE_VALUE;
    } else {
        IntegratedRealm = EABF_MEAN_FORCE_ERROR;
    }

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
        ReleaseAllResources();
        return(false);
    }
    if( SolveSystemOfEquations() == false ) {
        ReleaseAllResources();
        return(false);
    }

    vout << "   RMSR = " << setprecision(5) << GetRMSR() << endl;

// find global minimum
    double glb_min = X[0];
    for(int i=0; i < X.GetLength(); i++) {
        if(glb_min > X[i]) glb_min = X[i];
    }

// load data to FES
    for(int ipoint=0; ipoint < FES->GetNumberOfPoints(); ipoint++) {
        int x_index = XMap[ipoint];
        if(x_index >= 0) {
            double value = X[x_index]-glb_min;
            if( IntegratedRealm == EABF_MEAN_FORCE_VALUE ){
                FES->SetEnergy(ipoint,value);
            } else {
                FES->SetError(ipoint,value);
            }
            FES->SetNumOfSamples(ipoint,Accumulator->GetNumberOfABFSamples(ipoint));
        }
    }

    if( IntegratedRealm == EABF_MEAN_FORCE_VALUE ) {
        vout << "   SigmaF2 = " << setprecision(5) << FES->GetSigmaF2() << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRFD::BuildSystemOfEquations(CVerboseStr& vout)
{
// initializations
    IPoint.CreateVector(Accumulator->GetNumberOfCoords());
    XMap.CreateVector(Accumulator->GetNumberOfBins());

    XMap.Set(-1);

    NumOfVariables = 0;
    NumOfEquations = 0;
    NumOfNonZeros = 0;

    BuildEquations(true);

    vout << "   Number of variables           : " << NumOfVariables << std::endl;
    vout << "   Number of equations           : " << NumOfEquations << std::endl;
    vout << "   Number of non-zero A elements : " << NumOfNonZeros << std::endl;

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
    A = cs_spalloc(NumOfEquations,NumOfVariables,NumOfNonZeros,1,1);
    if(A == NULL) {
        CSmallString error;
        error << "unable to allocate A matrix (M: " << CSmallString(NumOfEquations) << ", N: " << CSmallString(NumOfVariables) << ", NZ: " << CSmallString(NumOfNonZeros) << ")";
        ES_ERROR(error);
        return(false);
    }

    Rhs.CreateVector(NumOfEquations);

// build system of equations
    BuildEquations(false);

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD::BuildEquations(bool trial)
{
    LocIter = 0;

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
                    if(trial == false) {
                        double dfac = 1.0 / (diff * 2.0);
                        cs_entry(A,LocIter,XMap[ifbin1],-3.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin2],+4.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],-1.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin1,IntegratedRealm);
                        LocIter++;
                        cs_entry(A,LocIter,XMap[ifbin1],-1.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],+1.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin2,IntegratedRealm);
                        LocIter++;
                        cs_entry(A,LocIter,XMap[ifbin1],+1.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin2],-4.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],+3.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin3,IntegratedRealm);
                        LocIter++;
                    } else {
                        NumOfEquations += 3;
                        NumOfNonZeros += 8;
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

                    if(trial == false) {
                        double dfac = 1.0 / (diff * 6.0);
                        cs_entry(A,LocIter,XMap[ifbin1],-11.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin2],+18.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],-9.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin4],+2.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin1,IntegratedRealm);
                        LocIter++;
                        cs_entry(A,LocIter,XMap[ifbin1],-2.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin2],-3.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],+6.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin4],-1.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin2,IntegratedRealm);
                        LocIter++;
                        cs_entry(A,LocIter,XMap[ifbin1],+1.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin2],-6.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],+3.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin4],+2.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin3,IntegratedRealm);
                        LocIter++;
                        cs_entry(A,LocIter,XMap[ifbin1],-2.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin2],+9.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin3],-18.0*dfac);
                        cs_entry(A,LocIter,XMap[ifbin4],+11.0*dfac);
                        Rhs[LocIter] = Accumulator->GetValue(ifcoord,ifbin4,IntegratedRealm);
                        LocIter++;
                    } else {
                        NumOfEquations += 4;
                        NumOfNonZeros += 16;
                    }
                    break;
                default:
                    break;
            }
        }
    }
}

//------------------------------------------------------------------------------

int CABFIntegratorRFD::GetFBinIndex(const CSimpleVector<int>& position,int ifcoord,int offset) const
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

bool CABFIntegratorRFD::SolveSystemOfEquations(void)
{
    cs* At;    // transposed A
    cs* AtA;   // AtA

    cA = cs_compress(A);
    if(cA == NULL) {
        CSmallString error;
        error << "unable to compress A matrix";
        ES_ERROR(error);
        return(false);
    }

    // keep A for RMSR determination
    cs_spfree(A);
    A = NULL;

    At = cs_transpose(cA,1);
    if(At == NULL) {
        CSmallString error;
        error << "unable to transpose A matrix";
        ES_ERROR(error);
        cs_spfree(cA);
        return(false);
    }

    AtA = cs_multiply(At,cA);
    if(At == NULL) {
        CSmallString error;
        error << "unable to create AtA matrix";
        ES_ERROR(error);
        cs_spfree(cA);
        cs_spfree(At);
        return(false);
    }

// skip - we needed it for RMSR
// release cA matrix
//    cs_spfree(cA);

// calculate right hand side vector
    X.CreateVector(AtA->n);

    X.SetZero();
    cs_gaxpy(At,Rhs,X);  // AtRhs = At * pRhs->x + AtRhs

// release At matrix
    cs_spfree(At);
    // keep rhs for rmsr calculation
//  Rhs.FreeVector();

// solve system of linear equations
    int order = 0;
    double tol = 0.01;

// solve Ax=b with LU
    if(cs_lusol(order,AtA,X,tol) == 0) {
        CSmallString error;
        error << "unable to solve system of equations";
        ES_ERROR(error);
        cs_spfree(AtA);
        return(false);
    }

// release AtA matrix
    cs_spfree(AtA);

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorRFD::ReleaseAllResources(void)
{
    XMap.FreeVector();
    IPoint.FreeVector();
    Rhs.FreeVector();
    X.FreeVector();
    if(A != NULL) cs_spfree(A);
    A = NULL;
    if(cA != NULL) cs_spfree(cA);
    cA = NULL;
}

//------------------------------------------------------------------------------

double CABFIntegratorRFD::GetRMSR(void)
{
    CSimpleVector<double> lhs;
    lhs.CreateVector(NumOfEquations);
    lhs.SetZero();

    // lhs = A * X + lhs
    if( cs_gaxpy(cA,X,lhs) == 0 ){
        ES_ERROR("cs_gaxpy(cA,X,lhs) - failed");
        return(0.0);
    }

    double rmsr = 0.0;
    double lv,rv,err;

    for(int i=0; i < NumOfEquations; i++){
        rv = Rhs[i];
        lv = lhs[i];
        err = lv-rv;
        rmsr = rmsr + err*err;
    }

    if( NumOfEquations <= 0 ) return(0.0);
    rmsr = rmsr / (double) NumOfEquations;
    if( rmsr > 0 ){
        rmsr = sqrt(rmsr);
    }
    return(rmsr);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
