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
#include <iomanip>

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegratorRBF::CABFIntegratorRBF(void)
{
    Accumulator = NULL;
    FES = NULL;

    Periodicity = false;
    WFac        = 2.0;
    RFac        = 3.0;
    Overhang    = 2;
    Method      = EARBF_SVD;

    IntegrateErrors = false;

    RCond   = -1; // machine precision

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

void CABFIntegratorRBF::SetWFac(double wfac)
{
    WFac = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetRFac(double rfac)
{
    RFac = rfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetOverhang(int nrbfs)
{
    Overhang = nrbfs;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetPeriodicity(bool set)
{
    Periodicity = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRBF::Integrate(CVerboseStr& vout,bool errors)
{
    IntegrateErrors = errors;

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

    if( (unsigned int)Accumulator->GetNumberOfCoords() != FES->GetNumberOfCoords() ){
        ES_ERROR("inconsistent ABF and FES");
        return(false);
    }

    // RBF setup
    NumOfCVs = Accumulator->GetNumberOfCoords();
    NumOfRBFBins.CreateVector(NumOfCVs);

    if( RFac <= 0 ) RFac = 1.0;

    NumOfRBFs = 1;
    for(int i=0; i < NumOfCVs; i++ ){
        int nrbfbins = Accumulator->GetCoordinate(i)->GetNumberOfBins()/RFac;
        if( ! Accumulator->GetCoordinate(i)->IsPeriodic() ){
            nrbfbins += 2*Overhang;
        }
        NumOfRBFBins[i] = nrbfbins;
        NumOfRBFs *= nrbfbins;
    }
    Weights.CreateVector(NumOfRBFs);
    Weights.SetZero();

    Sigmas.CreateVector(NumOfCVs);
    for(int i=0; i < NumOfCVs; i++){
        Sigmas[i] = Accumulator->GetCoordinate(i)->GetRange()*WFac/NumOfRBFBins[i];
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
        if( ! IntegrateErrors ){
            FES->SetEnergy(ipoint,value);
        } else {
            FES->SetError(ipoint,value);
        }
        FES->SetNumOfSamples(ipoint,Accumulator->GetNumberOfABFSamples(ipoint));
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

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    CSimpleVector<double>   lpos;
    lpos.CreateVector(NumOfCVs);

    vout << "   Dim: " << neq << " x " << NumOfRBFs << endl;

    CFortranMatrix A;
    A.CreateMatrix(neq,NumOfRBFs);
    A.SetZero();

    CVector rhs;
    int nrhs = std::max(neq,NumOfRBFs);
    rhs.CreateVector(nrhs);
    rhs.SetZero();

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
            rhs[j+k] = Accumulator->GetIntegratedValue(k,i,IntegrateErrors);
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

    vout << "   RMSR = " << setprecision(5) << GetRMSR() << endl;

    // debug
    // cout << Weights[0] << endl;

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
        double start = p_coord->GetMinValue();
        int nbins = NumOfRBFBins[k];
        if( ! p_coord->IsPeriodic() ){
            nbins -= 2*Overhang;
        }
        double step = (p_coord->GetMaxValue()-p_coord->GetMinValue())/((double)nbins);
        if( ! p_coord->IsPeriodic() ){
            start -= Overhang*step; // only half
        }
        position[k] = start + ((double)ibin + 0.5)*step;
        index = index / NumOfRBFBins[k];
    }
}

//------------------------------------------------------------------------------

double CABFIntegratorRBF::GetValue(const CSimpleVector<double>& position)
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

//------------------------------------------------------------------------------

double CABFIntegratorRBF::GetRMSR(void)
{
    if( Accumulator->GetNumberOfBins() <= 0 ){
        ES_ERROR("number of bins is not > 0");
        return(0.0);
    }

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    CSimpleVector<double>   lpos;
    lpos.CreateVector(NumOfCVs);

    CSimpleVector<double>   der;
    der.CreateVector(NumOfCVs);

    double rmsr = 0.0;

    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){

        Accumulator->GetPoint(i,ipos);
        der.SetZero();

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
                der[k] += Weights[l] * av * fc;
            }
        }

        for(int k=0; k < NumOfCVs; k++){
            double diff = Accumulator->GetIntegratedValue(k,i,IntegrateErrors) - der[k];
            rmsr += diff*diff;
        }
    }

    rmsr /= Accumulator->GetNumberOfBins();
    if( rmsr > 0.0 ){
        rmsr = sqrt(rmsr);
    }

    return(rmsr);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
