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

#include <math.h>
#include <errno.h>
#include "EntropyDer.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iomanip>
#include <Vector.hpp>
#include <FortranMatrix.hpp>
#include <GPFilter.hpp>
#include <SGFilter.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CEntropyDer)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEntropyDer::CEntropyDer(void)
{
    State       = 1;
    NSTLimit    = 0;
    NumOfBins   = 0;
    NumOfCVs    = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CEntropyDer::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-entropy-der (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgInAccuName() != "-") {
        vout << "# Input PMF accumulator (in)      : " << Options.GetArgInAccuName() << endl;
    } else {
        vout << "# Input PMF accumulator (in)      : - (standard input)" << endl;
    }
    if(Options.GetArgOutAccuName() != "-") {
        vout << "# Resulting PMF accumulator (out) : " << Options.GetArgOutAccuName() << endl;
    } else {
        vout << "# Resulting PMF accumulator (out) : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CEntropyDer::Run(void)
{
// load accumulator
    vout << endl;
    vout << format("%02d:Loading the input PMF accumulator ...")%State << endl;
    State++;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgInAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    InAccu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        InAccu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input PMF accumulator file");
        return(false);
    }
    InputFile.Close();

    NSTLimit    = InAccu->GetNSTLimit();
    NumOfBins   = InAccu->GetNumOfBins();
    NumOfCVs    = InAccu->GetNumOfCVs();

    vout << "   Done." << endl;

    InAccu->PrintInfo(vout);

// duplicate data
    OutAccu = InAccu->Duplicate();

    vout << endl;
    vout << format("%02d:Preparing input data ...")%State << endl;
    vout << "   Calculating Etot ..." << endl;
    GetEtot();

    vout << "   Numerical differentiation ..." << endl;
    //GetNativeData();
    //for(double wfac=10.0; wfac < 200; wfac += 10.0){
        GetICFByGPF(75.0);
    //}

    vout << "   Done." << endl;


    exit(1);

    vout << endl;
    vout << format("%02d:Processing data ...")%State << endl;
    State++;
    if( Options.IsOptKSKernelSet() || Options.IsOptKSWFacSet() ){
        vout << "   Employing kernel smoothing (KS) for data binning ..." << endl;
        vout << "      Kernel: " << Options.GetOptKSKernel() << endl;
        SetKSKernel(Options.GetOptKSKernel());
        vout << "      Wfac:   " << Options.GetOptKSWFac() << endl;
        SetKSWFac(Options.GetOptKSWFac());
        CalculatePPandPN_KS();
    } else {
        vout << "   Employing native approach for data binning ..."  << endl;
        CalculatePPandPN();
    }
    vout << "   Done." << endl;

// save results
    vout << endl;
    vout << format("%02d:Saving the resulting PMF accumulator ...")%State << endl;
    State++;

    if( OutputFile.Open(Options.GetArgOutAccuName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    try {
        OutAccu->Save(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the resulting PMF accumulator file");
        return(false);
    }
    OutputFile.Close();
    vout << "   Done." << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::GetEtot(void)
{
// input data
    CVectorDataPtr inepot  = InAccu->GetSectionData("TEPOT")->GetDataBlob();
    CVectorDataPtr inerst  = InAccu->GetSectionData("TERST")->GetDataBlob();
    CVectorDataPtr inekin  = InAccu->GetSectionData("TEKIN")->GetDataBlob();

    CVectorDataPtr outetot = InAccu->CreateSectionData("TETOT", "IG","R","S")->GetDataBlob();;

    for(size_t t=10; t < NSTLimit; t++){
        double epot = inepot->GetRawDataField()[t];
        double erst = inerst->GetRawDataField()[t];
        double ekin = inekin->GetRawDataField()[t];

        double etot = epot + erst + ekin;
        outetot->GetRawDataField()[t] = etot;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::GetICFBySGF(void)
{
    CSGFilterPtr sgfilter(new CSGFilter);
    sgfilter->SetFilter(InAccu->GetTimeStep(),5,3);

    CSimpleVector<double> cvsva;
    CSimpleVector<double> cvs1d;
    CSimpleVector<double> cvs2d;

    CSimpleVector<double> gfx;

    cvsva.CreateVector(NumOfCVs);
    cvs1d.CreateVector(NumOfCVs);
    cvs2d.CreateVector(NumOfCVs);

    gfx.CreateVector(NumOfCVs);

    CFortranMatrix zinvva;
    CFortranMatrix zinv1d;

    zinvva.CreateMatrix(NumOfCVs,NumOfCVs);
    zinv1d.CreateMatrix(NumOfCVs,NumOfCVs);

// input data
    CPMFAccuDataPtr insdcvs     = InAccu->GetSectionData("TCVS");
    CPMFAccuDataPtr insdzinv    = InAccu->GetSectionData("TZINV");
    CVectorDataPtr  inetot      = InAccu->GetSectionData("TETOT")->GetDataBlob();
    CPMFAccuDataPtr insdmicf    = InAccu->GetSectionData("TBICF");

    CVectorDataPtr  outetot     = InAccu->CreateSectionData("USE_ETOT","IG","R","T")->GetDataBlob();
    CPMFAccuDataPtr outsdcvs    = InAccu->CreateSectionData("USE_CVS", "IG","R","S");
    CPMFAccuDataPtr outsdicf    = InAccu->CreateSectionData("USE_ICF", "IG","R","S");

    for(size_t t=2000; t < NSTLimit-5; t++){

        outetot->GetRawDataField()[t] = sgfilter->GetValue(inetot,t);

    // get force components
        for(size_t icv=0; icv < NumOfCVs; icv++){
            CVectorDataPtr incvs  = insdcvs->GetDataBlob(icv);

            cvsva[icv] = sgfilter->GetValue(incvs,t);
            cvs1d[icv] = sgfilter->Get1stDer(incvs,t);
            cvs2d[icv] = sgfilter->Get2ndDer(incvs,t);

            for(size_t jcv=0; jcv < NumOfCVs; jcv++){
                CVectorDataPtr inzinv  = insdzinv->GetDataBlob(icv,jcv);
                zinvva[icv][jcv] = sgfilter->GetValue(inzinv,t);
                zinv1d[icv][jcv] = sgfilter->Get1stDer(inzinv,t);
            }
        }

    // complete forces
        gfx.SetZero();
        for(size_t icv=0; icv < NumOfCVs; icv++){
            for(size_t jcv=0; jcv < NumOfCVs; jcv++){
                gfx[icv] = gfx[icv] + zinv1d[icv][jcv]*cvs1d[jcv] + zinvva[icv][jcv]*cvs2d[jcv];
            }
        }

    // apply biasing force
        for(size_t icv=0; icv < NumOfCVs; icv++){
            CVectorDataPtr inmicf  = insdmicf->GetDataBlob(icv);
            gfx[icv] = gfx[icv] - inmicf->GetRawDataField()[t+2];
        }

    // save
        for(size_t icv=0; icv < NumOfCVs; icv++){
            CVectorDataPtr outcvs  = outsdcvs->GetDataBlob(icv);
            outcvs->GetRawDataField()[t] = cvsva[icv];
            CVectorDataPtr outicf  = outsdicf->GetDataBlob(icv);
            outicf->GetRawDataField()[t] = gfx[icv];
        }
    }

    ofstream tout("icf");
    CVectorDataPtr inicf  = InAccu->GetSectionData("TICF")->GetDataBlob(0);
    CVectorDataPtr outicf = outsdicf->GetDataBlob(0);
    for(int t=200; t< 50000; t++){
        tout << inicf->GetRawDataField()[t] << " " <<  outicf->GetRawDataField()[t+2] << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::GetICFByGPF(double wfac)
{
    CGPFilterPtr gpfilter(new CGPFilter);

    int    gplen    = 301;
    double gpwfac   = 153;
    double gpnoise  = 0;


    gpfilter->SetKernel("ardmc32");
    gpfilter->SetFilter(InAccu->GetTimeStep(),gplen,gpwfac,gpnoise,vout);

    CSimpleVector<double> cvsva;
    CSimpleVector<double> cvs1d;
    CSimpleVector<double> cvs2d;

    CSimpleVector<double> gfx;

    cvsva.CreateVector(NumOfCVs);
    cvs1d.CreateVector(NumOfCVs);
    cvs2d.CreateVector(NumOfCVs);

    gfx.CreateVector(NumOfCVs);

    CFortranMatrix zinvva;
    CFortranMatrix zinv1d;

    zinvva.CreateMatrix(NumOfCVs,NumOfCVs);
    zinv1d.CreateMatrix(NumOfCVs,NumOfCVs);

// input data
    CPMFAccuDataPtr insdcvs     = InAccu->GetSectionData("TCVS");
    CPMFAccuDataPtr insdzinv    = InAccu->GetSectionData("TZINV");
    CVectorDataPtr  inetot      = InAccu->GetSectionData("TETOT")->GetDataBlob();
    CPMFAccuDataPtr insdmicf    = InAccu->GetSectionData("TBICF");

    CVectorDataPtr  outetot     = InAccu->CreateSectionData("USE_ETOT","IG","R","T")->GetDataBlob();
    CPMFAccuDataPtr outsdcvs    = InAccu->CreateSectionData("USE_CVS", "IG","R","S");
    CPMFAccuDataPtr outsdicf    = InAccu->CreateSectionData("USE_ICF", "IG","R","S");

    double totlogml = 0.0;
    double lsize    = 0.0;

   for(size_t t=10; t < NSTLimit-gplen; t += gplen){

    // get force components
        for(size_t icv=0; icv < NumOfCVs; icv++){
            CVectorDataPtr incvs  = insdcvs->GetDataBlob(icv);

            gpfilter->TrainProcess(incvs,t);
            double logml = gpfilter->GetLogML();
            totlogml += logml;
            lsize++;

            CVectorDataPtr outcvs  = outsdcvs->GetDataBlob(icv);
            gpfilter->PredictData(outcvs,t);
        }

//        gpfilter->TrainProcess(inetot,t);
//        double logml = gpfilter->GetLogPL();
//        totlogml += logml;
//        lsize++;

        gpfilter->PredictData(outetot,t);
   }

    totlogml /= lsize;
    vout << "   wfac = " << wfac << " avelogml= " << totlogml << endl;

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::SetKSWFac(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetKSWFac");
    }

    string          sspec(spec);
    vector<string>  swfacs;

    split(swfacs,sspec,is_any_of("x"),token_compress_on);

    if( swfacs.size() > NumOfCVs ){
        CSmallString error;
        error << "too many kswfacs (" << swfacs.size() << ") than required (" << NumOfCVs << ")";
        RUNTIME_ERROR(error);
    }

    KSWFac.CreateVector(NumOfCVs);

    // parse values of wfac
    double last_wfac = 1.0;
    for(size_t i=0; i < swfacs.size(); i++){
        stringstream str(swfacs[i]);
        str >> last_wfac;
        if( ! str ){
            CSmallString error;
            error << "unable to decode kswfac value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        KSWFac[i] = last_wfac;
    }

    // pad the rest with the last value
    for(size_t i=swfacs.size(); i < NumOfCVs; i++){
        KSWFac[i] = last_wfac;
    }
}

//------------------------------------------------------------------------------

void CEntropyDer::SetKSKernel(const CSmallString& kskernel)
{
    if( kskernel == "parabolic" ){
        KSKernel = EKSKT_PARABOLIC;
    } else if( kskernel == "triweight" ) {
        KSKernel = EKSKT_TRIWEIGHT;
    } else if( kskernel == "gaussian" ) {
        KSKernel = EKSKT_GAUSSIAN;
    } else if( kskernel == "default" ) {
        KSKernel = EKSKT_TRIWEIGHT;
    } else {
        CSmallString error;
        error << "Specified kernel '" << kskernel << "' is not supported. "
                 "Supported kernels are: parabolic, triweight, gaussian, and default(=triweight)";
        INVALID_ARGUMENT(error);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::CalculatePPandPN(void)
{
    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

// input data
    CPMFAccuDataPtr incvs  = InAccu->GetSectionData("USE_CVS");
    CPMFAccuDataPtr inetot = InAccu->GetSectionData("USE_ETOT");
    CPMFAccuDataPtr inicf  = InAccu->GetSectionData("USE_ICF");

// output data
    CPMFAccuDataPtr outnsamples    = OutAccu->CreateSectionData("TDS_NSAMPLES",  "AD","R","B");
    CPMFAccuDataPtr outmicf        = OutAccu->CreateSectionData("MTDS_ICF",      "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2icf       = OutAccu->CreateSectionData("M2TDS_ICF",     "M2","R","M","TDS_NSAMPLES","MTDS_ICF");
    CPMFAccuDataPtr outmetot       = OutAccu->CreateSectionData("MTDS_ETOT",     "WA","R","B","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2etot      = OutAccu->CreateSectionData("M2TDS_ETOT",    "M2","R","B","TDS_NSAMPLES","MTDS_ETOT");
    CPMFAccuDataPtr outmpp         = OutAccu->CreateSectionData("MTDS_PP",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pp        = OutAccu->CreateSectionData("M2TDS_PP",      "M2","R","M","TDS_NSAMPLES","MTDS_PP");
    CPMFAccuDataPtr outmpn         = OutAccu->CreateSectionData("MTDS_PN",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pn        = OutAccu->CreateSectionData("M2TDS_PN",      "M2","R","M","TDS_NSAMPLES","MTDS_PN");

    for(size_t t=10; t < NSTLimit; t++){
        incvs->GetDataBlob(t,ipos);
        int gi0 = InAccu->GetGlobalIndex(ipos);
        if( gi0 == -1 ) continue;

        // cout << t << " " << gi0 << " " << ipos[0] << endl;

        double etot = inetot->GetData(t);

        // increase number of samples
        double n    = outnsamples->GetData(gi0);
        n++;
        outnsamples->SetData(gi0,n);

        double inv = 1.0 / n;

        double metot  = outmetot->GetData(gi0);
        double m2etot = outm2etot->GetData(gi0);

        double detot1 = etot - metot;
        metot = metot + detot1 * inv;
        double detot2 = etot - metot;
        m2etot = m2etot + detot1 * detot2;

        outmetot->SetData(gi0,metot);
        outm2etot->SetData(gi0,m2etot);

        for(size_t icv=0; icv < NumOfCVs; icv++){
            double icf  = inicf->GetData(t,icv);

            double pp   = etot+icf;
            double pn   = etot-icf;

        // ---------------------------------------
            double micf  = outmicf->GetData(gi0,icv);
            double m2icf = outm2icf->GetData(gi0,icv);

            double dicf1 = icf - micf;
            micf  = micf + dicf1 * inv;
            double dicf2 = icf - micf;
            m2icf = m2icf + dicf1 * dicf2;

            outmicf->SetData(gi0,icv,micf);
            outm2icf->SetData(gi0,icv,m2icf);

        // ---------------------------------------
            double mpp  = outmpp->GetData(gi0,icv);
            double m2pp = outm2pp->GetData(gi0,icv);

            double dpp1 = pp - mpp;
            mpp  = mpp + dpp1 * inv;
            double dpp2 = pp - mpp;
            m2pp = m2pp + dpp1 * dpp2;

            outmpp->SetData(gi0,icv,mpp);
            outm2pp->SetData(gi0,icv,m2pp);

        // ---------------------------------------
            double mpn  = outmpn->GetData(gi0,icv);
            double m2pn = outm2pn->GetData(gi0,icv);

            double dpn1 = pn - mpn;
            mpn  = mpn + dpn1 * inv;
            double dpn2 = pn - mpn;
            m2pn = m2pn + dpn1 * dpn2;

            outmpn->SetData(gi0,icv,mpn);
            outm2pn->SetData(gi0,icv,m2pn);
        }
    }
}

//------------------------------------------------------------------------------

void CEntropyDer::CalculatePPandPN_KS(void)
{
    CSimpleVector<double> weights;
    weights.CreateVector(NumOfBins);

    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

// input data
    CPMFAccuDataPtr incvs  = InAccu->GetSectionData("USE_CVS");
    CPMFAccuDataPtr inetot = InAccu->GetSectionData("USE_ETOT");
    CPMFAccuDataPtr inicf  = InAccu->GetSectionData("USE_ICF");

// output data
    CPMFAccuDataPtr outnsamples    = OutAccu->CreateSectionData("TDS_NSAMPLES",  "AD","R","B");
    CPMFAccuDataPtr outmicf        = OutAccu->CreateSectionData("MTDS_ICF",      "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2icf       = OutAccu->CreateSectionData("M2TDS_ICF",     "M2","R","M","TDS_NSAMPLES","MTDS_ICF");
    CPMFAccuDataPtr outmetot       = OutAccu->CreateSectionData("MTDS_ETOT",     "WA","R","B","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2etot      = OutAccu->CreateSectionData("M2TDS_ETOT",    "M2","R","B","TDS_NSAMPLES","MTDS_ETOT");
    CPMFAccuDataPtr outmpp         = OutAccu->CreateSectionData("MTDS_PP",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pp        = OutAccu->CreateSectionData("M2TDS_PP",      "M2","R","M","TDS_NSAMPLES","MTDS_PP");
    CPMFAccuDataPtr outmpn         = OutAccu->CreateSectionData("MTDS_PN",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pn        = OutAccu->CreateSectionData("M2TDS_PN",      "M2","R","M","TDS_NSAMPLES","MTDS_PN");

    for(size_t t=10; t < NSTLimit; t++){

        incvs->GetDataBlob(t,ipos);
        CalculateKSWeights(ipos,weights);

        double etot = inetot->GetData(t);

        for(size_t ibin=0; ibin < NumOfBins; ibin++){
        // increase number of samples
            double w    = weights[ibin];
            double n    = outnsamples->GetData(ibin);
            n = n + w;
            outnsamples->SetData(ibin,n);

            if( w <= 0.0 ) continue;

            double inv = w / n;

            double metot  = outmetot->GetData(ibin);
            double m2etot = outm2etot->GetData(ibin);

            double detot1 = etot - metot;
            metot = metot + detot1 * inv;
            double detot2 = etot - metot;
            m2etot = m2etot + detot1 * detot2;

            outmetot->SetData(ibin,metot);
            outm2etot->SetData(ibin,m2etot);

            for(size_t icv=0; icv < NumOfCVs; icv++){
                double icf  = inicf->GetData(t,icv);
                double pp   = etot+icf;
                double pn   = etot-icf;

            // ---------------------------------------
                double micf  = outmicf->GetData(ibin,icv);
                double m2icf = outm2icf->GetData(ibin,icv);

                double dicf1 = icf - micf;
                micf  = micf + dicf1 * inv;
                double dicf2 = icf - micf;
                m2icf = m2icf + dicf1 * dicf2;

                outmicf->SetData(ibin,icv,micf);
                outm2icf->SetData(ibin,icv,m2icf);

            // ---------------------------------------
                double mpp  = outmpp->GetData(ibin,icv);
                double m2pp = outm2pp->GetData(ibin,icv);

                double dpp1 = pp - mpp;
                mpp  = mpp + dpp1 * inv;
                double dpp2 = pp - mpp;
                m2pp = m2pp + w * dpp1 * dpp2;

                outmpp->SetData(ibin,icv,mpp);
                outm2pp->SetData(ibin,icv,m2pp);

            // ---------------------------------------
                double mpn  = outmpn->GetData(ibin,icv);
                double m2pn = outm2pn->GetData(ibin,icv);

                double dpn1 = pn - mpn;
                mpn  = mpn + dpn1 * inv;
                double dpn2 = pn - mpn;
                m2pn = m2pn + w * dpn1 * dpn2;

                outmpn->SetData(ibin,icv,mpn);
                outm2pn->SetData(ibin,icv,m2pn);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CEntropyDer::CalculateKSWeights(CSimpleVector<double>& jpos,CSimpleVector<double>& weights)
{
    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

    weights.SetZero();

    double tot_w = 0.0;

    for(size_t indi=0; indi < NumOfBins; indi++) {
        InAccu->GetPoint(indi,ipos);
        double u2 = 0.0;
        for(size_t icv=0; icv < NumOfCVs; icv++){
            double dx = (InAccu->GetCV(icv)->GetDifference(ipos[icv],jpos[icv])) / (KSWFac[icv] * InAccu->GetCV(icv)->GetBinWidth());
            u2 += dx*dx;
        }
        double w = GetKSKernelValue(u2);
        weights[indi] = w;
        tot_w += w;
    }

    if( tot_w <= 0.0 ) return;

    // normalize
    for(size_t indi=0; indi < NumOfBins; indi++) {
        weights[indi] /= tot_w;
    }
}

//------------------------------------------------------------------------------

double CEntropyDer::GetKSKernelValue(double u2)
{
    switch(KSKernel){
        case(EKSKT_PARABOLIC):
            if( u2 < 1.0 ){
                return( 3.0/4.0*(1-u2)*(1-u2) );
            }
            return(0.0);
        case(EKSKT_TRIWEIGHT):
            if( u2 < 1.0 ){
                return( 35.0/32.0*(1-u2)*(1-u2)*(1-u2) );
            }
            return(0.0);
        case(EKSKT_GAUSSIAN):
            return( 1.0/sqrt(2.0*M_PI)*exp(-0.5*u2) );
        default:
            RUNTIME_ERROR("not implemented");
            break;
    }
    return(0.0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-entropy-der terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

