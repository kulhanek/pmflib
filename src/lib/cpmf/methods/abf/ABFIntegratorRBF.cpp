// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <ABFIntegratorRBF.hpp>
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

CABFIntegratorRBF::CABFIntegratorRBF(void)
{
    Accumulator = NULL;
    FES = NULL;

    Verbose = true;
    Periodicity = false;
    WidthOrder = 4;
    Limit = 500;
    Method = EARBF_SVD;

    RCond = -1; // machine precision
}

//------------------------------------------------------------------------------

CABFIntegratorRBF::~CABFIntegratorRBF(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorRBF::SetInputABFAccumulator(const CABFAccumulator* p_accu)
{
    Accumulator = p_accu;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetOutputFESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetVerbosity(bool set)
{
    Verbose = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetGaussianWidth(int order)
{
    WidthOrder = order;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetPeriodicity(bool set)
{
    Periodicity = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRBF::Integrate(CVerboseStr& vout)
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

    // RBF setup
    NumOfCVs = Accumulator->GetNumberOfCoords();
    NumOfRBFBins.CreateVector(NumOfCVs);

    double sfac = 5;

    NumOfRBFs = 1;
    for(int i=0; i < NumOfCVs; i++ ){
        NumOfRBFBins[i] = Accumulator->GetCoordinate(i)->GetNumberOfBins()/sfac;
        NumOfRBFs *= Accumulator->GetCoordinate(i)->GetNumberOfBins()/sfac;
    }
    Weights.CreateVector(NumOfRBFs);
    Weights.SetZero();

    Sigmas.CreateVector(NumOfCVs);
    for(int i=0; i < NumOfCVs; i++){
        Sigmas[i] = Accumulator->GetCoordinate(i)->GetRange()*WidthOrder/NumOfRBFBins[i];
        // cout << Sigmas[i] << endl;
    }

    // integrate
    if( IntegrateByLS(vout) == false ){
        ES_ERROR("unable to solve least square problem by svd");
        return(false);
    }

    vout << "   Calculating FES ..." << endl;

    // set-up FES
    FES->Allocate(Accumulator);

    // load data to FES
    CSimpleVector<double>   position;
    position.CreateVector(NumOfCVs);

    for(unsigned int i=0; i < FES->GetNumberOfPoints(); i++) {
        FES->GetPoint(i,position);
        double value = GetEnergy(position);
        FES->SetEnergy(i,value);
        FES->SetNumOfSamples(i,Accumulator->GetNumberOfABFSamples(i));
    }

    // find global minimum
    double glb_min = FES->GetEnergy(0);
    for(unsigned int i=0; i < FES->GetNumberOfPoints(); i++) {
        if(glb_min > FES->GetEnergy(i)) glb_min = FES->GetEnergy(i);
    }

    // move global minimum to zero
    for(unsigned int i=0; i < FES->GetNumberOfPoints(); i++) {
        FES->SetEnergy(i,FES->GetEnergy(i) - glb_min);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRBF::IntegrateByLS(CVerboseStr& vout)
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

    vout << "   Dim: " << neq << " x " << NumOfRBFs << endl;

    CFortranMatrix A;
    A.CreateMatrix(neq,NumOfRBFs);

    CVector rhs;
    int nrhs = std::max(neq,NumOfRBFs);
    rhs.CreateVector(nrhs);

    // calculate A and rhs
    int j=0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){

        Accumulator->GetPoint(i,ipos);
        // A
        for(int l=0; l < NumOfRBFs; l++){
            GetRBFPosition(l,lpos);
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
            rhs[j+k] = Accumulator->GetABFForceSum(k,i);
        }
        j += NumOfCVs;
    }

    int result = 0;

    switch(Method){
        case EARBF_SVD:{
            vout << "   Solving least square problem by SVD ..." << endl;

            // solve least square problem via GELSD
            int rank = 0;
            result = CLapack::gelsd(A,rhs,RCond,rank);
            vout << "   Rank = " << rank << "; Info = " << result << endl;
            if( result != 0 ) return(false);
            }
        break;
        case EARBF_QRLQ:{
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
    for(int l=0; l < NumOfRBFs; l++){
        Weights[l] = rhs[l];
    }

    cout << Weights[0] << endl;

    return( true );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorRBF::GetRBFPosition(unsigned int index,CSimpleVector<double>& position)
{
    for(int k=NumOfCVs-1; k >= 0; k--) {

        const CColVariable* p_coord = Accumulator->GetCoordinate(k);
        int ibin = index % NumOfRBFBins[k];
        // bin  = 0,...,NBins-1
        position[k] = p_coord->GetMinValue() + ((double)ibin + 0.5)*(p_coord->GetMaxValue()-p_coord->GetMinValue())/((double)NumOfRBFBins[k]);
        index = index / NumOfRBFBins[k];
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorRBF::GetEnergy(const CSimpleVector<double>& position)
{
    double                  energy = 0.0;
    CSimpleVector<double>   rbfpos;

    rbfpos.CreateVector(NumOfCVs);

    for(int i=0; i < NumOfRBFs; i++){
        GetRBFPosition(i,rbfpos);
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
