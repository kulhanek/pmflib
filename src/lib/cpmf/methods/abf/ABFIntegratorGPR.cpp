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
#include <Lapack.hpp>

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegratorGPR::CABFIntegratorGPR(void)
{
    Accumulator = NULL;
    FES = NULL;

    Verbose = true;
    Periodicity = false;
    WidthOrder = 4;
    Method = EAGPR_SVD;

    RCond = -1; // machine precision
}

//------------------------------------------------------------------------------

CABFIntegratorGPR::~CABFIntegratorGPR(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorGPR::SetInputABFAccumulator(const CABFAccumulator* p_accu)
{
    Accumulator = p_accu;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetOutputFESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetVerbosity(bool set)
{
    Verbose = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetGaussianWidth(int order)
{
    WidthOrder = order;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetPeriodicity(bool set)
{
    Periodicity = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::Integrate(CVerboseStr& vout)
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
        ES_ERROR("inconsistent ABF and FES");
        return(false);
    }

    // GPR setup
    NumOfCVs = Accumulator->GetNumberOfCoords();
    NumOfGPRBins.CreateVector(NumOfCVs);

    double sfac = 5;

    NumOfGPRs = 1;
    for(int i=0; i < NumOfCVs; i++ ){
        NumOfGPRBins[i] = Accumulator->GetCoordinate(i)->GetNumberOfBins()/sfac;
        NumOfGPRs *= Accumulator->GetCoordinate(i)->GetNumberOfBins()/sfac;
    }
    Weights.CreateVector(NumOfGPRs);
    Weights.SetZero();

    Sigmas.CreateVector(NumOfCVs);
    for(int i=0; i < NumOfCVs; i++){
        Sigmas[i] = Accumulator->GetCoordinate(i)->GetRange()*WidthOrder/NumOfGPRBins[i];
        // cout << Sigmas[i] << endl;
    }

    // integrate
    if( IntegrateByLS(vout) == false ){
        ES_ERROR("unable to solve least square problem by svd");
        return(false);
    }

    vout << "   Calculating FES ..." << endl;

    // load data to FES
    CSimpleVector<double>   position;
    position.CreateVector(NumOfCVs);

    // find global minimum
    double glb_min = 0.0;
    for(unsigned int ipoint=0; ipoint < FES->GetNumberOfPoints(); ipoint++) {
        FES->GetPoint(ipoint,position);
        double value = GetValue(position);
        if( (ipoint == 0) || (glb_min > value) ){
            glb_min = value;
        }
    }

    // write results
    for(unsigned int ipoint=0; ipoint < FES->GetNumberOfPoints(); ipoint++) {
        FES->GetPoint(ipoint,position);
        double value = GetValue(position) - glb_min;
        FES->SetEnergy(ipoint,value);
        FES->SetNumOfSamples(ipoint,Accumulator->GetNumberOfABFSamples(ipoint));
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::IntegrateByLS(CVerboseStr& vout)
{
    vout << "   Creating A and rhs ..." << endl;

    // number of equations
    int neq = Accumulator->GetNumberOfBins()*NumOfCVs;

    CSimpleVector<double>   spos;           // best sampled bin
    spos.CreateVector(NumOfCVs);

    CSimpleVector<double>   ipos;           // best sampled bin
    ipos.CreateVector(NumOfCVs);

    CSimpleVector<double>   lpos;           // best sampled bin
    lpos.CreateVector(NumOfCVs);

    vout << "   Dim: " << neq << " x " << NumOfGPRs << endl;

    CFortranMatrix A;
    A.CreateMatrix(neq,NumOfGPRs);

    CVector rhs;
    int nrhs = std::max(neq,NumOfGPRs);
    rhs.CreateVector(nrhs);

    // calculate A and rhs
    int j=0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){

        Accumulator->GetPoint(i,ipos);
        // A
        for(int l=0; l < NumOfGPRs; l++){
            GetGPRPosition(l,lpos);
                  //      cout << lpos[0] << " " << lpos[1] << endl;
            double av = 1.0;
            for(int k=0; k < NumOfCVs; k++){
                double dvc = Accumulator->GetCoordinate(k)->GetDifference(ipos[k],lpos[k]);
                double sig = Sigmas[k];
                av *= exp( - dvc*dvc/(2.0*sig*sig) );
            }
            for(int k=0; k < NumOfCVs; k++){
                double dvc = Accumulator->GetCoordinate(k)->GetDifference(ipos[k],lpos[k]);
                double sig = Sigmas[k];
                double fc = -dvc/(sig*sig);  // switch to derivatives
                A[j+k][l] = av * fc;
            }
        }
        // rhs
        for(int k=0; k < NumOfCVs; k++){
            rhs[j+k] = Accumulator->GetIntegratedValue(k,i,false);
        }
        j += NumOfCVs;
    }

    int result = 0;

    switch(Method){
        case EAGPR_SVD:{
            vout << "   Solving least square problem by SVD ..." << endl;

            // solve least square problem via GELSD
            int rank = 0;
            result = CLapack::gelsd(A,rhs,RCond,rank);
            vout << "   Rank = " << rank << "; Info = " << result << endl;
            if( result != 0 ) return(false);
            }
        break;
        case EAGPR_QRLQ:{
            vout << "   Solving least square problem by QRLQ ..." << endl;

            // solve least square problem via GELS
            result = CLapack::gels(A,rhs);
            if( result != 0 ) return(false);
            }
        break;
        default:
            INVALID_ARGUMENT("unsupported method");
    }

    // copy results to Weights
    for(int l=0; l < NumOfGPRs; l++){
        Weights[l] = rhs[l];
    }

    // debug
    // cout << Weights[0] << endl;

    return( true );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorGPR::GetGPRPosition(unsigned int index,CSimpleVector<double>& position)
{
    for(int k=NumOfCVs-1; k >= 0; k--) {

        const CColVariable* p_coord = Accumulator->GetCoordinate(k);
        int ibin = index % NumOfGPRBins[k];
        // bin  = 0,...,NBins-1
        position[k] = p_coord->GetMinValue() + ((double)ibin + 0.5)*(p_coord->GetMaxValue()-p_coord->GetMinValue())/((double)NumOfGPRBins[k]);
        index = index / NumOfGPRBins[k];
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetValue(const CSimpleVector<double>& position)
{
    double                  energy = 0.0;
    CSimpleVector<double>   rbfpos;

    rbfpos.CreateVector(NumOfCVs);

    for(int i=0; i < NumOfGPRs; i++){
        GetGPRPosition(i,rbfpos);
        double value = Weights[i];
        for(int j=0; j < NumOfCVs; j++){
            double dvc = Accumulator->GetCoordinate(j)->GetDifference(position[j],rbfpos[j]);
            double sig = Sigmas[j];
            value *= exp( - dvc*dvc/(2.0*sig*sig) );
        }
        energy = energy + value;
    }
    return(energy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
