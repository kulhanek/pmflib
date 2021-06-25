// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <IntegratorRBF.hpp>
#include <ErrorSystem.hpp>
#include <FortranMatrix.hpp>
#include <Vector.hpp>
#include <algorithm>
#include <SciLapack.hpp>
#include <iomanip>
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

CIntegratorRBF::CIntegratorRBF(void)
{
    Overhang    = 2;
    Method      = ERBFLLS_SVD;

    RCond   = -1; // machine precision

    IncludeGluedBins = false;

    NumOfBins   = 0;
    NCVs        = 0;
    NumOfRBFs   = 0;
    NumOfValues = 0;

    GlobalMinSet    = false;
    NoEnergy        = false;

    NumOfThreads    = 1;
}

//------------------------------------------------------------------------------

CIntegratorRBF::~CIntegratorRBF(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CIntegratorRBF::SetInputEnergyDerProxy(CEnergyDerProxyPtr p_proxy)
{
    DerProxy = p_proxy;
    if( DerProxy ){
        Accu = DerProxy->GetAccu();
    } else {
        Accu = CPMFAccumulatorPtr();
    }

    if( Accu ){
        NCVs = Accu->GetNumOfCVs();
        NumOfBins = Accu->GetNumOfBins();
    } else {
        NCVs = 0;
        NumOfBins = 0;
    }
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetOutputES(CEnergySurfacePtr p_surf)
{
    EneSurf = p_surf;
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetWFac(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetWFac");
    }

    string          sspec(spec);
    vector<string>  swfacs;

    split(swfacs,sspec,is_any_of("x"),token_compress_on);

    if( swfacs.size() > NCVs ){
        CSmallString error;
        error << "too many wfacs (" << swfacs.size() << ") than required (" << NCVs << ")";
        RUNTIME_ERROR(error);
    }

    WFac.CreateVector(NCVs);

    // parse values of wfac
    double last_wfac = 3.0;
    for(size_t i=0; i < swfacs.size(); i++){
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
    for(size_t i=swfacs.size(); i < NCVs; i++){
        WFac[i] = last_wfac;
    }
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetRFac(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetRFac");
    }

    string          sspec(spec);
    vector<string>  srfacs;

    split(srfacs,sspec,is_any_of("x"),token_compress_on);

    if( srfacs.size() > NCVs ){
        CSmallString error;
        error << "too many rfacs (" << srfacs.size() << ") than required (" << NCVs << ")";
        RUNTIME_ERROR(error);
    }

    RFac.CreateVector(NCVs);

    // parse values of rfac
    double last_rfac = 1.0;
    for(size_t i=0; i < srfacs.size(); i++){
        stringstream str(srfacs[i]);
        str >> last_rfac;
        if( ! str ){
            CSmallString error;
            error << "unable to decode rfac value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        RFac[i] = last_rfac;
    }

    // pad the rest with the last value
    for(size_t i=srfacs.size(); i < NCVs; i++){
        RFac[i] = last_rfac;
    }
}
//------------------------------------------------------------------------------

void CIntegratorRBF::SetOverhang(int nrbfs)
{
    Overhang = nrbfs;
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetLLSMethod(ERBFLLSMethod set)
{
    Method = set;
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetLLSMethod(const CSmallString& method)
{
    if( method == "svd" ){
        SetLLSMethod(ERBFLLS_SVD);
    } else if( method == "qr" ) {
        SetLLSMethod(ERBFLLS_QR);
    } else if( method == "default" ) {
        SetLLSMethod(ERBFLLS_SVD);
    } else {
        CSmallString error;
        error << "Specified method '" << method << "' for linear algebra is not supported. "
                 "Supported methods are: svd (divide and conquer driver), qr, default (=svd)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetNoEnergy(bool set)
{
    NoEnergy = set;
}

//------------------------------------------------------------------------------

void CIntegratorRBF::IncludeGluedAreas(bool set)
{
    IncludeGluedBins = set;
}

//------------------------------------------------------------------------------

void CIntegratorRBF::SetGlobalMin(const CSmallString& spec)
{
    GlobalMinSet = true;
    string sspec(spec);
    if( Accu == NULL ){
        RUNTIME_ERROR("accumulator is not set for SetGlobalMin");
    }

    // remove "x" from the string
    replace (sspec.begin(), sspec.end(), 'x' , ' ');

    // parse values of CVs
    GPos.CreateVector(NCVs);
    stringstream str(sspec);
    for(size_t i=0; i < NCVs; i++){
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

bool CIntegratorRBF::Integrate(CVerboseStr& vout)
{
    PrintExecInfo(vout);

    if( Accu == NULL ) {
        ES_ERROR("ABF accumulator is not set");
        return(false);
    }
    if( EneSurf == NULL ) {
        ES_ERROR("EneSurf is not set");
        return(false);
    }

    if( Accu->GetNumOfCVs() == 0 ) {
        RUNTIME_ERROR("number of coordinates is zero");
    }
    if( Accu->GetNumOfCVs() != EneSurf->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent ABF and EneSurf - CVs");
    }
    if( Accu->GetNumOfBins() != EneSurf->GetNumOfPoints() ){
        RUNTIME_ERROR("inconsistent ABF and EneSurf - points");
    }
    if( NCVs == 0 ) {
        ES_ERROR("number of coordinates is zero");
        return(false);
    }
    if( RFac.GetLength() != NCVs ){
        ES_ERROR("rfac not set");
        return(false);
    }
    if( WFac.GetLength() != NCVs ){
        ES_ERROR("wfac not set");
        return(false);
    }

    // RBF setup
    NumOfRBFBins.CreateVector(NCVs);

    for(size_t i=0; i < NCVs; i++ ){
        if( RFac[i] <= 0 ) RFac[i] = 1.0;
    }

    NumOfRBFs = 1;
    for(size_t i=0; i < NCVs; i++ ){
        size_t nrbfbins;
        nrbfbins = Accu->GetCV(i)->GetNumOfBins()/RFac[i];
        if( ! Accu->GetCV(i)->IsPeriodic() ){
            nrbfbins += 2*Overhang;
        }
        NumOfRBFBins[i] = nrbfbins;
        NumOfRBFs *= nrbfbins;
    }

    Weights.CreateVector(NumOfRBFs);
    Weights.SetZero();

    Sigmas.CreateVector(NCVs);
    for(size_t i=0; i < NCVs; i++){
        if( NumOfRBFBins[i] == 0 ){
            RUNTIME_ERROR("NumOfRBFBins[i] is zero");
        }
        Sigmas[i] = Accu->GetCV(i)->GetRange()*WFac[i]/NumOfRBFBins[i];
        // cout << Sigmas[i] << endl;
    }

    // print hyperparameters
    vout << "   RBF parameters ..." << endl;
    for(size_t k=0; k < NCVs; k++ ){
        vout << format("      WFac#%-2d   = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(size_t k=0; k < NCVs; k++ ){
        vout << format("      RFac#%-2d   = %10.4f")%(k+1)%RFac[k] << endl;
    }

    // integrate
    if( IntegrateByLS(vout) == false ){
        ES_ERROR("unable to solve least square problem by svd");
        return(false);
    }

    for(size_t k=0; k < NCVs; k++ ){
    vout << "      RMSR CV#" << k+1 << " = " << setprecision(5) << GetRMSR(k) << endl;
    }

    if( ! NoEnergy ){
        CalculateEnergy(vout);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CIntegratorRBF::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   Calculating EneSurf ..." << endl;

// create map for bins with calculated energy and error
    NumOfValues = 0;
    for(size_t i=0; i < NumOfBins; i++){
        int samples = DerProxy->GetNumOfSamples(i);
        if( IncludeGluedBins ){
            if( samples == 0 ) continue;
        } else {
            if( samples <= 0 ) continue;
        }
        NumOfValues++;
    }
    ValueMap.CreateVector(NumOfValues);

    size_t indi = 0;
    for(size_t i=0; i < NumOfBins; i++){
        int samples = DerProxy->GetNumOfSamples(i);
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
    for(size_t indj=0; indj < NumOfValues; indj++){
        size_t j = ValueMap[indj];
        Accu->GetPoint(j,jpos);
        values[indj] = GetValue(jpos);
    }

// basic EneSurf update
    for(size_t i=0; i < NumOfBins; i++){
        int samples = DerProxy->GetNumOfSamples(i);
        EneSurf->SetNumOfSamples(i,samples);
        EneSurf->SetEnergy(i,0.0);
        EneSurf->SetError(i,0.0);
    }

// update EneSurf
    if( GlobalMinSet ){
        // GPos.CreateVector(NCVs) - is created in  SetGlobalMin
   //   vout << "   Calculating EneSurf ..." << endl;
        vout << "      Global minimum provided at: ";
        vout << GPos[0];
        for(int i=1; i < Accu->GetNumOfCVs(); i++){
            vout << "x" << GPos[i];
        }
        double glb_min = GetValue(GPos);
        vout << " (" << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurf->SetEnergy(j,values[indj]-glb_min);
        }
    } else {
        // search for global minimum
        GPos.CreateVector(NCVs);
        bool   first = true;
        double glb_min = 0.0;
        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            int samples = DerProxy->GetNumOfSamples(j);
            if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
            double value = values[indj];
            if( first || (glb_min > value) ){
                glb_min = value;
                first = false;
                Accu->GetPoint(j,GPos);
            }
        }

   //   vout << "   Calculating EneSurf ..." << endl;
        vout << "      Global minimum found at: ";
        vout << GPos[0];
        for(size_t i=1; i < NCVs; i++){
            vout << "x" << GPos[i];
        }
        vout << " (" << glb_min << ")" << endl;

        for(size_t indj=0; indj < NumOfValues; indj++){
            size_t j = ValueMap[indj];
            EneSurf->SetEnergy(j,values[indj]-glb_min);
        }
    }

    vout << "      SigmaF2   = " << setprecision(5) << EneSurf->GetSigmaF2() << endl;
    if( IncludeGluedBins ){
        vout << "      SigmaF2 (including glued bins) = " << setprecision(5) << EneSurf->GetSigmaF2(true) << endl;
    }
}

//------------------------------------------------------------------------------

void CIntegratorRBF::PrintExecInfo(CVerboseStr& vout)
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

void CIntegratorRBF::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CIntegratorRBF::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CIntegratorRBF::IntegrateByLS(CVerboseStr& vout)
{
    vout << "   Creating A and rhs ..." << endl;

    // number of equations
    NumOfUsedBins = 0;
    for(size_t i=0; i < NumOfBins; i++){
        if( DerProxy->GetNumOfSamples(i) > 0 ) NumOfUsedBins++;
    }
    NumOfEq = NumOfUsedBins*NCVs;

    // create sampled map
    SampledMap.CreateVector(NumOfUsedBins);
    size_t ind = 0;
    for(size_t i=0; i < NumOfBins; i++){
        if( DerProxy->GetNumOfSamples(i) <= 0 ) continue;
        SampledMap[ind] = i;
        ind++;
    }

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NCVs);

    CSimpleVector<double>   lpos;
    lpos.CreateVector(NCVs);

    vout << "      Dim: " << NumOfEq << " x " << NumOfRBFs << endl;

    CFortranMatrix A;
    A.CreateMatrix(NumOfEq,NumOfRBFs);
    A.SetZero();

    CVector rhs;
    size_t nrhs = std::max(NumOfEq,NumOfRBFs);
    rhs.CreateVector(nrhs);
    rhs.SetZero();

    // calculate A and rhs
    #pragma omp parallel for firstprivate(ipos,lpos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accu->GetPoint(i,ipos);
        // A
        for(size_t l=0; l < NumOfRBFs; l++){
            GetRBFPosition(l,lpos);
                  //      cout << lpos[0] << " " << lpos[1] << endl;
            double av = 1.0;
            for(size_t k=0; k < NCVs; k++){
                double dvc = Accu->GetCV(k)->GetDifference(ipos[k],lpos[k]);
                double sig = Sigmas[k];
                av *= exp( - dvc*dvc/(2.0*sig*sig) );
            }
            for(size_t k=0; k < NCVs; k++){
                double dvc = Accu->GetCV(k)->GetDifference(ipos[k],lpos[k]);
                double sig = Sigmas[k];
                double fc = -dvc/(sig*sig);  // switch to derivatives
                A[indi*NCVs+k][l] = av * fc;
            }
        }
        // rhs
        for(size_t k=0; k < NCVs; k++){
            rhs[indi*NCVs+k] = DerProxy->GetValue(i,k,E_PROXY_VALUE);
        }
    }

    int result = 0;

    RunBlasLapackPar();

    switch(Method){
        case ERBFLLS_SVD:{
            vout << "   Solving least square problem by SVD ..." << endl;

            // solve least square problem via GELSD
            int rank = 0;
            double realRCond = 0.0;
            result = CSciLapack::gelsd(A,rhs,RCond,rank,realRCond);
            vout << "      Rank = " << rank << "; Info = " << result << "; Real rcond = " << scientific << realRCond << fixed << endl;
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
    #pragma omp parallel for
    for(size_t l=0; l < NumOfRBFs; l++){
        Weights[l] = rhs[l];
    }

    return( true );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CIntegratorRBF::GetRBFPosition(size_t index,CSimpleVector<double>& position)
{
    for(long int k=NCVs-1; k >= 0; k--) {

        const CColVariablePtr p_coord = Accu->GetCV(k);
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

double CIntegratorRBF::GetValue(const CSimpleVector<double>& position)
{
    double                  energy = 0.0;

    CSimpleVector<double>   rbfpos;
    rbfpos.CreateVector(NCVs);

    // parallelized one level up
    for(size_t i=0; i < NumOfRBFs; i++){
        GetRBFPosition(i,rbfpos);
        double value = Weights[i];
        for(size_t j=0; j < NCVs; j++){
            double dvc = Accu->GetCV(j)->GetDifference(position[j],rbfpos[j]);
            double sig = Sigmas[j];
            value *= exp( - dvc*dvc/(2.0*sig*sig) );
        }
        energy = energy + value;
    }

    return(energy);
}

//------------------------------------------------------------------------------

double CIntegratorRBF::GetRMSR(size_t cv)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(0.0);
    }

    CSimpleVector<double>   ipos;
    ipos.CreateVector(NCVs);

    double rmsr = 0.0;

    #pragma omp parallel for firstprivate(ipos) reduction(+:rmsr)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accu->GetPoint(i,ipos);

        double mfi = DerProxy->GetValue(i,cv,E_PROXY_VALUE);
        double mfp = GetMeanForce(ipos,cv);
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

bool CIntegratorRBF::WriteMFInfo(const CSmallString& name)
{
    if( NumOfBins == 0 ){
        ES_ERROR("number of bins is not > 0");
        return(false);
    }

    CSimpleVector<double> jpos;
    CSimpleVector<double> mfi;
    CSimpleVector<double> mfp;

    jpos.CreateVector(NCVs);
    mfi.CreateVector(NumOfEq);
    mfp.CreateVector(NumOfEq);

    // calculate
    #pragma omp parallel for firstprivate(jpos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];
        Accu->GetPoint(i,jpos);
        for(size_t k=0; k < NCVs; k++){
            mfi[indi*NCVs+k] = DerProxy->GetValue(i,k,E_PROXY_VALUE);
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
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accu->GetPoint(i,jpos);

        for(size_t c=0; c < NCVs; c++){
            ofs << format("%20.16f ")%jpos[c];
        }

        for(size_t k=0; k < NCVs; k++){
            ofs << format(" %20.16f %20.16f")%mfi[indi*NCVs+k]%mfp[indi*NCVs+k];
        }

        ofs << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CIntegratorRBF::FilterByMFZScore(double zscore,CVerboseStr& vout)
{
    vout << high;
    vout << "   Running Z-test filter on MF errors ..." << endl;

    if( NumOfBins == 0 ){
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

    flags.CreateVector(NumOfBins);
    sig2.CreateVector(NCVs);
    maxi.CreateVector(NCVs);
    maxzscore.CreateVector(NCVs);
    mferror2.CreateVector(NumOfEq);

    flags.Set(1);

    vout << high;
    vout << "      Precalculating MF errors ..." << endl;

    CSimpleVector<double> jpos;
    jpos.CreateVector(NCVs);

    // precalculate values
    #pragma omp parallel for firstprivate(jpos)
    for(size_t indi=0; indi < NumOfUsedBins; indi++){
        size_t i = SampledMap[indi];

        Accu->GetPoint(i,jpos);

        for(size_t k=0; k < NCVs; k++){
            double diff2 = DerProxy->GetValue(i,k,E_PROXY_VALUE) - GetMeanForce(jpos,k);
            diff2 *= diff2;
            mferror2[indi*NCVs+k] = diff2;
        }
    }

    vout << "      Searching for MF outliers ..." << endl;
    vout << debug;

    bool testme = true;
    while( testme ){

        testme = false;

        double  count = 0;
        sig2.SetZero();

        // calc variances - we assume zero mean on errors
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            if( flags[i] != 0 ) {
                for(size_t k=0; k < NCVs; k++){
                    sig2[k] += mferror2[indi*NCVs+k];
                }
                count++;
            }
        }

        if( count == 0 ) break;
        for(size_t k=0; k < NCVs; k++){
            sig2[k] /= count;
        }

        // filter
        bool first = true;
        for(size_t indi=0; indi < NumOfUsedBins; indi++){
            size_t i = SampledMap[indi];

            if( flags[i] != 0 ){

                for(size_t k=0; k < NCVs; k++){
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

        for(size_t k=0; k < NCVs; k++){
            if( maxzscore[k] > zscore ){
                flags[maxi[k]] = 0;
                testme = true;
                vout << "   outlier found at " << maxi[k] << " for cv " << k << " with z-score " << sqrt(maxzscore[k]) << endl;
            }
        }
    }

    // apply limits
    size_t outliers = 0;
    for(size_t i=0; i < NumOfBins; i++){
        if( flags[i] == 0 ){
            DerProxy->SetNumOfSamples(i,0);
            outliers++;
        }
    }

    // from now the integrator is in invalid state !!!
    vout << high;
    vout << "      Number of outliers = " << outliers << endl;
}

//------------------------------------------------------------------------------

double CIntegratorRBF::GetMeanForce(const CSimpleVector<double>& ipos,size_t icoord)
{
    CSimpleVector<double>   lpos;
    lpos.CreateVector(NCVs);

    double mf = 0.0;
    for(size_t l=0; l < NumOfRBFs; l++){
        GetRBFPosition(l,lpos);
              //      cout << lpos[0] << " " << lpos[1] << endl;
        double av = 1.0;
        for(size_t k=0; k < NCVs; k++){
            double dvc = Accu->GetCV(k)->GetDifference(ipos[k],lpos[k]);
            double sig = Sigmas[k];
            av *= exp( - dvc*dvc/(2.0*sig*sig) );
        }
        double dvc = Accu->GetCV(icoord)->GetDifference(ipos[icoord],lpos[icoord]);
        double sig = Sigmas[icoord];
        double fc = -dvc/(sig*sig);  // switch to derivatives
        mf += Weights[l] * av * fc;
    }

    return(mf);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
