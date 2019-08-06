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

using namespace std;
using namespace boost;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegratorGPR::CABFIntegratorGPR(void)
{
    Accumulator = NULL;
    FES = NULL;

    SigmaF2 = 15.0;
    WFac1 = 3.0;
    WFac2 = 0.0;

    IncludeError = false;
    IncludeGluedBins = false;
    NoEnergy = false;

    Method = EGPRINV_LU;
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

void CABFIntegratorGPR::SetWFac1(double wfac)
{
    WFac1 = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetWFac2(double wfac)
{
    WFac2 = wfac;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetSigmaF2(double sigf2)
{
    SigmaF2 = sigf2;
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::SetIncludeError(bool set)
{
    IncludeError = set;
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
        ES_ERROR("inconsistent ABF and FES - CVs");
        return(false);
    }
    if( Accumulator->GetNumberOfBins() != FES->GetNumberOfPoints() ){
        ES_ERROR("inconsistent ABF and FES - points");
        return(false);
    }

    // GPR setup
    CVLengths2.CreateVector(Accumulator->GetNumberOfCoords());
    for(int i=0; i < Accumulator->GetNumberOfCoords(); i++){
        double l = WFac1*Accumulator->GetCoordinate(i)->GetRange()/Accumulator->GetCoordinate(i)->GetNumberOfBins();
        CVLengths2[i] = l*l;
    }
    if( (WFac2 > 0.0) && (Accumulator->GetNumberOfCoords() == 2) ){
        int cv = 1;
        double l = WFac2*Accumulator->GetCoordinate(cv)->GetRange()/Accumulator->GetCoordinate(cv)->GetNumberOfBins();
        CVLengths2[cv] = l*l;
    }
    if( (WFac2 > 0.0) && (Accumulator->GetNumberOfCoords()> 2) ){
        RUNTIME_ERROR("wfac2 > 0 and ncvs > 2");
    }

    // number of data points
    GPRSize = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) > 0 ) GPRSize++;
    }
    NCVs = Accumulator->GetNumberOfCoords();
    GPRSize = GPRSize*NCVs;

    // load data to FES
    ipos.CreateVector(NCVs);
    jpos.CreateVector(NCVs);
    gpos.CreateVector(NCVs);
    rk.CreateVector(GPRSize);
    lk.CreateVector(GPRSize);
    ik.CreateVector(GPRSize);
    GPRModel.CreateVector(GPRSize);

    // train GPR
    if( TrainGP(vout) == false ){
        ES_ERROR("unable to train GPR model");
        return(false);
    }

    // and finaly some statistics
    for(int k=0; k < Accumulator->GetNumberOfCoords(); k++ ){
    vout << "   RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
    }

    // and marginal likelihood
    vout << "   logML = " << setprecision(5) << GetLogMarginalLikelihood() << endl;

    if( ! NoEnergy ){
        vout << "   Calculating FES ..." << endl;

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
            Accumulator->GetPoint(i,jpos);
            value = GetValue(jpos);
            FES->SetEnergy(i,value);
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                gpos = jpos;
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

        vout << "   SigmaF2 = " << setprecision(5) << FES->GetSigmaF2() << endl;
        if( IncludeGluedBins ){
            vout << "   SigmaF2 (including glued bins) = " << setprecision(5) << FES->GetSigmaF2(true) << endl;
        }

        if( IncludeError ){
            CalculateErrors(gpos,vout);
        }
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegratorGPR::TrainGP(CVerboseStr& vout)
{
    vout << "   Creating K+Sigma and Y ..." << endl;
    vout << "   Dim: " << GPRSize << " x " << GPRSize << endl;

    K.CreateMatrix(GPRSize,GPRSize);
    K.SetZero();

    // create kernel matrix
    int indi = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue; // not sampled

        Accumulator->GetPoint(i,ipos);

        int indj = 0;
        for(int j=0; j < Accumulator->GetNumberOfBins(); j++){
            if( Accumulator->GetNumberOfABFSamples(j) <= 0 ) continue; // not sampled

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
            indj++;
        }
        indi++;
    }

// get mean forces and their variances
    Y.CreateVector(GPRSize);

    indi = 0.0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue; // not sampled

        for(int ii=0; ii < NCVs; ii++){

            // get mean force and its error
            double mf = Accumulator->GetValue(ii,i,EABF_MEAN_FORCE_VALUE);
            double er = Accumulator->GetValue(ii,i,EABF_MEAN_FORCE_ERROR);

            Y[indi*NCVs+ii] = mf;
            K[indi*NCVs+ii][indi*NCVs+ii] += er*er;
        }
        indi++;
    }

// inverting the K+Sigma
    int result = 0;
    switch(Method){
        case(EGPRINV_LU):
            vout << "   Inverting K+Sigma by LU ..." << endl;
            result = CSciLapack::inv1(K,logdetK);
            if( result != 0 ) return(false);
            break;
        case(EGPRINV_SVD):{
            vout << "   Inverting K+Sigma by SVD ..." << endl;
            int rank = 0;
            result = CSciLapack::inv2(K,logdetK,RCond,rank);
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CABFIntegratorGPR::GetValue(const CSimpleVector<double>& position)
{
    double energy = 0.0;

    int indi = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue; // not sampled

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
        indi++;
    }

    return(energy);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetMeanForce(const CSimpleVector<double>& position,int icoord)
{
    double mf = 0.0;

    int indi = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue; // not sampled

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
        indi++;
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

    double rmsr = 0.0;
    double nsamples = 0.0;

    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,jpos);

        double mfi = Accumulator->GetValue(cv,i,EABF_MEAN_FORCE_VALUE);
        double mfp = GetMeanForce(jpos,cv);
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

bool CABFIntegratorGPR::WriteMFInfo(const CSmallString& name)
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

    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,jpos);

        for(int c=0; c < Accumulator->GetNumberOfCoords(); c++){
            ofs << format("%20.16f ")%jpos[c];
        }

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double mfi = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE);
            double mfp = GetMeanForce(jpos,k);

            ofs << format(" %20.16f %20.16f")%mfi%mfp;
        }

        ofs << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::FilterByMFFac(double mffac)
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

    // calc variances
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,jpos);

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double diff2 = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(jpos,k);
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

        Accumulator->GetPoint(i,jpos);

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double diff2 = Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(jpos,k);
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

void CABFIntegratorGPR::FilterByMFMaxError(double mfmaxerr)
{
    CSimpleVector<int>  flags;
    flags.CreateVector(Accumulator->GetNumberOfBins());
    flags.Set(1);

    // filter
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,jpos);

        for(int k=0; k < Accumulator->GetNumberOfCoords(); k++){
            double err = fabs(Accumulator->GetValue(k,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(jpos,k));
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

void CABFIntegratorGPR::FilterByMFMaxError(double mfmaxerr1,double mfmaxerr2)
{
    if( Accumulator->GetNumberOfCoords() != 2 ){
        LOGIC_ERROR("only two CV can be combined with FilterByMFMaxError(double mfmaxerr1,double mfmaxerr2)");
    }

    CSimpleVector<int>  flags;
    flags.CreateVector(Accumulator->GetNumberOfBins());
    flags.Set(1);

    // filter
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue;

        Accumulator->GetPoint(i,jpos);

        double err = fabs(Accumulator->GetValue(0,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(jpos,0));
        if( err > mfmaxerr1 ){
            flags[i] = 0;
        }

        err = fabs(Accumulator->GetValue(1,i,EABF_MEAN_FORCE_VALUE) - GetMeanForce(jpos,1));
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

double CABFIntegratorGPR::GetCov(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos)
{
    int indi = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue; // not sampled

        Accumulator->GetPoint(i,ipos);

        double argl = 0.0;
        double argr = 0.0;
        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            argl += du*du/(2.0*dd);

            du = Accumulator->GetCoordinate(ii)->GetDifference(rpos[ii],ipos[ii]);
            argr += du*du/(2.0*dd);
        }
        argl = SigmaF2*exp(-argl);
        argr = SigmaF2*exp(-argr);

        for(int ii=0; ii < NCVs; ii++){
            double du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],ipos[ii]);
            double dd = CVLengths2[ii];
            lk[indi*NCVs + ii] = (du/dd)*argl;

            du = Accumulator->GetCoordinate(ii)->GetDifference(rpos[ii],ipos[ii]);
            rk[indi*NCVs + ii] = (du/dd)*argr;
        }
        indi++;
    }

    double cov = 0.0;
    for(int ii=0; ii < NCVs; ii++){
        double du = Accumulator->GetCoordinate(ii)->GetDifference(lpos[ii],rpos[ii]);
        double dd = CVLengths2[ii];
        cov += du*du/(2.0*dd);
    }
    cov = SigmaF2*exp(-cov);

    CSciBlas::gemv(1.0,K,rk,0.0,ik);
    cov = cov  - CSciBlas::dot(lk,ik);

    return(cov);
}

//------------------------------------------------------------------------------

double CABFIntegratorGPR::GetVar(CSimpleVector<double>& lpos)
{
    int indi = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        if( Accumulator->GetNumberOfABFSamples(i) <= 0 ) continue; // not sampled

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
        indi++;
    }

    double cov = SigmaF2;

    CSciBlas::gemv(1.0,K,lk,0.0,ik);
    cov = cov  - CSciBlas::dot(lk,ik);

    return(cov);
}

//------------------------------------------------------------------------------

void CABFIntegratorGPR::CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout)
{
    vout << "   Calculating FES error ..." << endl;
    CSmallTime st;
    st.GetActualTime();

    double vargp = GetVar(gpos);

    int totbatches = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        if( samples <= 0 ) continue;
        totbatches++;
    }

    int nbatches = 0;
    for(int i=0; i < Accumulator->GetNumberOfBins(); i++){
        int samples = Accumulator->GetNumberOfABFSamples(i);
        FES->SetError(i,0.0);
        if( samples <= 0 ) continue;
        Accumulator->GetPoint(i,jpos);
        double varfc = GetVar(jpos);
        double covfg = GetCov(jpos,gpos);
        // vout << varfc << " " << vargp << " " << covfg << endl;
        double error = varfc + vargp - 2.0*covfg;
        if( error > 0 ){
            error = sqrt(error);
        } else {
            error = 0.0;
        }
        FES->SetError(i,error);
        nbatches++;
        CSmallTime ct;
        ct.GetActualTime();
        if( (ct - st).GetSecondsFromBeginning() > 5*60 ){
            int comp = nbatches*100 / totbatches;
            vout << format("      completed %6d/%6d - %2d%%")%nbatches%totbatches%comp << endl;
            st = ct;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
