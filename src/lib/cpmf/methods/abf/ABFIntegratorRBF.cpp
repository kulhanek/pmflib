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
#include <SciLapack.hpp>
#include <iomanip>
#include <boost/format.hpp>

using namespace std;
using namespace boost;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegratorRBF::CABFIntegratorRBF(void)
{
    Accumulator = NULL;
    FES = NULL;

    WFac1       = 2.0;
    RFac1       = 3.0;
    WFac2       = 0.0;
    RFac2       = 0.0;
    Overhang    = 2;
    Method      = ERBFLLS_SVD;

    RCond   = -1; // machine precision

    IncludeGluedBins = false;
}

//------------------------------------------------------------------------------

CABFIntegratorRBF::~CABFIntegratorRBF(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegratorRBF::SetInputABFAccumulator(CABFAccumulator* p_accu)
{
    Accumulator = p_accu;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetOutputFESurface(CEnergySurface* p_surf)
{
    FES = p_surf;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetWFac1(double wfac)
{
    WFac1 = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetWFac2(double wfac)
{
    WFac2 = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetRFac1(double rfac)
{
    RFac1 = rfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetRFac2(double rfac)
{
    RFac2 = rfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetOverhang(int nrbfs)
{
    Overhang = nrbfs;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::SetLLSMehod(ERBFLLSMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::IncludeGluedAreas(bool set)
{
    IncludeGluedBins = set;
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

    if( Accumulator->GetNumberOfCoords() != FES->GetNumberOfCoords() ){
        ES_ERROR("inconsistent ABF and FES - CVs");
        return(false);
    }
    if( Accumulator->GetNumberOfBins() != FES->GetNumberOfPoints() ){
        ES_ERROR("inconsistent ABF and FES - points");
        return(false);
    }

    // RBF setup
    NumOfCVs = Accumulator->GetNumberOfCoords();
    NumOfRBFBins.CreateVector(NumOfCVs);

    if( (RFac2 > 0.0) && (Accumulator->GetNumberOfCoords() > 2) ){
        RUNTIME_ERROR("rfac2 > 0 and ncvs > 2");
    }

   if( RFac1 <= 0 ) RFac1 = 1.0;

    NumOfRBFs = 1;
    for(int i=0; i < NumOfCVs; i++ ){
        int nrbfbins;
        nrbfbins = Accumulator->GetCoordinate(i)->GetNumberOfBins()/RFac1;
        if( (i == 1) && (RFac2 > 0.0) ){
            nrbfbins = Accumulator->GetCoordinate(i)->GetNumberOfBins()/RFac2;
        }
        if( ! Accumulator->GetCoordinate(i)->IsPeriodic() ){
            nrbfbins += 2*Overhang;
        }
        NumOfRBFBins[i] = nrbfbins;
        NumOfRBFs *= nrbfbins;
    }

    if( (WFac2 > 0.0) && (Accumulator->GetNumberOfCoords() > 2) ){
        RUNTIME_ERROR("wfac2 > 0 and ncvs > 2");
    }

    Weights.CreateVector(NumOfRBFs);
    Weights.SetZero();

    Sigmas.CreateVector(NumOfCVs);
    for(int i=0; i < NumOfCVs; i++){
        Sigmas[i] = Accumulator->GetCoordinate(i)->GetRange()*WFac1/NumOfRBFBins[i];
        // cout << Sigmas[i] << endl;
    }
    if( (WFac2 > 0.0) && (Accumulator->GetNumberOfCoords() == 2) ){
        int cv = 1;
        Sigmas[cv] = Accumulator->GetCoordinate(cv)->GetRange()*WFac2/NumOfRBFBins[cv];
    }

    // integrate
    if( IntegrateByLS(vout) == false ){
        ES_ERROR("unable to solve least square problem by svd");
        return(false);
    }

    vout << "   Calculating FES ..." << endl;

    // load data to FES
    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    // calculate energies
    double glb_min = 0.0;
    bool   first = true;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        FES->SetNumOfSamples(i,samples);
        double value = 0.0;
        FES->SetEnergy(i,value);
        if( IncludeGluedBins ){
            if( samples == 0 ) continue;
        } else {
            if( samples <= 0 ) continue;
        }
        Accumulator->GetPoint(i,ipos);
        value = GetValue(ipos);
        FES->SetEnergy(i,value);
        if( first || (glb_min > value) ){
            glb_min = value;
            first = false;
        }
    }

    // move to global minima
    for(int i=0; i < FES->GetNumberOfPoints(); i++) {
        if( IncludeGluedBins ){
            if( FES->GetNumOfSamples(i) == 0 ) continue;
        } else {
            if( FES->GetNumOfSamples(i) <= 0 ) continue;
        }
        double value = 0.0;
        value = FES->GetEnergy(i);
        value = value - glb_min;
        FES->SetEnergy(i,value);
    }

    for(int k=0; k < Accumulator->GetNumberOfCoords(); k++ ){
    vout << "   RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
    }

    vout << "   SigmaF2 = " << setprecision(5) << FES->GetSigmaF2() << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorRBF::IntegrateByLS(CVerboseStr& vout)
{
    vout << "   Creating A and rhs ..." << endl;

    // number of equations
    int neq = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) > 0 ) neq++;
    }
    neq = neq*NumOfCVs;

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
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

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
            rhs[j+k] = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE);
        }
        j += NumOfCVs;
    }

    int result = 0;

    switch(Method){
        case ERBFLLS_SVD:{
            vout << "   Solving least square problem by SVD ..." << endl;

            // solve least square problem via GELSD
            int rank = 0;
            result = CSciLapack::gelsd(A,rhs,RCond,rank);
            vout << "   Rank = " << rank << "; Info = " << result << endl;
            if( result != 0 ) return(false);
            }
        break;
        case ERBFLLS_QR:{
            vout << "   Solving least square problem by QR ..." << endl;

            // solve least square problem via GELS
            result = CSciLapack::gels(A,rhs);
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

double CABFIntegratorRBF::GetRMSR(int cv)
{
    if( Accumulator->GetNumberOfBins() <= 0 ){
        ES_ERROR("number of bins is not > 0");
        return(0.0);
    }

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    double rmsr = 0.0;
    double nsamples = 0.0;

    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){

        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,ipos);

        double mfi = Accumulator->GetValue(cv,i,EABF_MEAN_FORCE_VALUE);
        double mfp = GetMeanForce(ipos,cv);
        double diff = mfi - mfp;

        rmsr += diff*diff;
        nsamples++;
    }

    if( nsamples > 0 ){
        rmsr /= nsamples;
    }
    if( rmsr > 0.0 ){
        rmsr = sqrt(rmsr);
    }

    return(rmsr);
}

//------------------------------------------------------------------------------

bool CABFIntegratorRBF::WriteMFInfo(const CSmallString& name)
{
    if( Accumulator->GetNumberOfBins() <= 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    ofstream ofs(name);
    if( ! ofs ){
        CSmallString error;
        error << "unable to open file '" << name << "' for derivatives";
        ES_ERROR(error);
        return(false);
    }

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){

        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,ipos);

        for(int c=0; c < Accumulator->GetNumberOfCoords(); c++){
            ofs << format("%20.16f ")%ipos[c];
        }

        for(int k=0; k < NumOfCVs; k++){
            double mfi = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE);
            double mfp = GetMeanForce(ipos,k);
            ofs << format(" %20.16f %20.16f")%mfi%mfp;
        }

        ofs << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::FilterByMFFac(double mffac)
{
    if( Accumulator->GetNumberOfBins() <= 0 ){
        ES_ERROR("number of bins is not > 0");
        return;
    }

    // we work with variances
    mffac *= mffac;

    // we assume zero mean on errors
    CSimpleVector<double>   sig2;
    sig2.CreateVector(Accumulator->GetNumberOfCoords());
    double                  count = 0;
    sig2.SetZero();

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    // calc variances
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,ipos);

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double diff2 = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(ipos,k);
            diff2 *= diff2;
            sig2[k] += diff2;
        }
        count++;
    }

    if( count == 0 ) return;
    for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
        sig2[k] /= count;
    }

    CSimpleVector<int>  flags;
    flags.CreateVector(Accumulator->GetNumberOfBins());
    flags.Set(1);

    // filter
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,ipos);

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double diff2 = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(ipos,k);
            diff2 *= diff2;
            if( diff2 > mffac*sig2[k] ){
                flags[i] = 0;
            }
        }
    }

    // apply limits
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( flags[i] == 0 ){
            Accumulator->SetNumberOfABFSamples(i,0);
        }
    }

    // from now the integrator is in invalid state !!!
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::FilterByMFMaxError(double mfmaxerr)
{
    CSimpleVector<int>  flags;
    flags.CreateVector(Accumulator->GetNumberOfBins());
    flags.Set(1);

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    // filter
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,ipos);

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double err = fabs(Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(ipos,k));
            if( err > mfmaxerr ){
                flags[i] = 0;
            }
        }
    }

    // apply limits
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( flags[i] == 0 ){
            Accumulator->SetNumberOfABFSamples(i,0);
        }
    }

    // from now the integrator is in invalid state !!!
}

//------------------------------------------------------------------------------

void CABFIntegratorRBF::FilterByMFMaxError(double mfmaxerr1,double mfmaxerr2)
{
    if( Accumulator->GetNumberOfCoords() != 2 ){
        LOGIC_ERROR("only two CV can be combined with FilterByMFMaxError(double mfmaxerr1,double mfmaxerr2)");
    }

    CSimpleVector<int>  flags;
    flags.CreateVector(Accumulator->GetNumberOfBins());
    flags.Set(1);

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NumOfCVs);

    // filter
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,ipos);

        double err = fabs(Accumulator->GetValue(0,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(ipos,0));
        if( err > mfmaxerr1 ){
            flags[i] = 0;
        }

        err = fabs(Accumulator->GetValue(1,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(ipos,1));
        if( err > mfmaxerr2 ){
            flags[i] = 0;
        }

    }

    // apply limits
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( flags[i] == 0 ){
            Accumulator->SetNumberOfABFSamples(i,0);
        }
    }

    // from now the integrator is in invalid state !!!
}

//------------------------------------------------------------------------------

double CABFIntegratorRBF::GetMeanForce(const CSimpleVector<double>& ipos,int icoord)
{
    CSimpleVector<double>   lpos;
    lpos.CreateVector(NumOfCVs);

    double mf = 0.0;
    for(int l=0; l < NumOfRBFs; l++){
        GetRBFPosition(l,lpos);
              //      cout << lpos[0] << " " << lpos[1] << endl;
        double av = 1.0;
        for(int k=0; k < NumOfCVs; k++){
            double dvc = Accumulator->GetCoordinate(k)->GetDifference(ipos[k],lpos[k]);
            double sig = Sigmas[k];
            av *= exp( - dvc*dvc/(2.0*sig*sig) );
        }
        double dvc = Accumulator->GetCoordinate(icoord)->GetDifference(ipos[icoord],lpos[icoord]);
        double sig = Sigmas[icoord];
        double fc = -dvc/(sig*sig);  // switch to derivatives
        mf += Weights[l] * av * fc;
    }

    return(mf);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
