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

    KSWfac.CreateVector(1);
    KSWfac[0] = 2.0;
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

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgInAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgOutAccuName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }
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

    InAccu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        InAccu->Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input PMF accumulator file");
        return(false);
    }
    NSTLimit    = InAccu->GetNSTLimit();
    NumOfBins   = InAccu->GetNumOfBins();
    NumOfCVs    = InAccu->GetNumOfCVs();

    vout << "   Done." << endl;

    InAccu->PrintInfo(vout);

// duplicate data
    OutAccu = InAccu->Duplicate();

    vout << endl;
    vout << format("%02d:Processing data ...")%State << endl;
    State++;
    // CalculatePPandPN();
    CalculatePPandPN_KS();
    vout << "   Done." << endl;

// save results
    vout << endl;
    vout << format("%02d:Saving the resulting PMF accumulator ...")%State << endl;
    State++;

    try {
        OutAccu->Save(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the resulting PMF accumulator file");
        return(false);
    }
    vout << "   Done." << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::CalculatePPandPN(void)
{
    CSimpleVector<double> ipos;
    ipos.CreateVector(NumOfCVs);

// input data
    CPMFAccuDataPtr incvs  = InAccu->GetSectionData("TCVS");
    CPMFAccuDataPtr inepot = InAccu->GetSectionData("TEPOT");
    CPMFAccuDataPtr inerst = InAccu->GetSectionData("TERST");
    CPMFAccuDataPtr inekin = InAccu->GetSectionData("TEKIN");
    CPMFAccuDataPtr inicf  = InAccu->GetSectionData("TICF");

// output data
    CPMFAccuDataPtr outnsamples    = OutAccu->CreateSectionData("TDS_NSAMPLES",  "AD","R","B");
    CPMFAccuDataPtr outmpp         = OutAccu->CreateSectionData("MTDS_PP",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pp        = OutAccu->CreateSectionData("M2TDS_PP",      "M2","R","M","TDS_NSAMPLES","MTDS_PP");
    CPMFAccuDataPtr outmpn         = OutAccu->CreateSectionData("MTDS_PN",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pn        = OutAccu->CreateSectionData("M2TDS_PN",      "M2","R","M","TDS_NSAMPLES","MTDS_PN");

    for(size_t t=10; t < NSTLimit; t++){
        incvs->GetData(t,ipos);
        int gi0 = InAccu->GetGlobalIndex(ipos);
        if( gi0 == -1 ) continue;

        // cout << t << " " << gi0 << " " << ipos[0] << endl;

        double epot = inepot->GetData(t);
        double erst = inerst->GetData(t);
        double ekin = inekin->GetData(t);

        // increase number of samples
        double n    = outnsamples->GetData(gi0);
        n++;
        outnsamples->SetData(gi0,n);

        double inv = 1.0 / n;

        for(size_t icv=0; icv < NumOfCVs; icv++){
            double icf  = inicf->GetData(t,icv);
            double pp   = epot+erst+ekin+icf;
            double pn   = epot+erst+ekin-icf;

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
    CPMFAccuDataPtr incvs  = InAccu->GetSectionData("TCVS");
    CPMFAccuDataPtr inepot = InAccu->GetSectionData("TEPOT");
    CPMFAccuDataPtr inerst = InAccu->GetSectionData("TERST");
    CPMFAccuDataPtr inekin = InAccu->GetSectionData("TEKIN");
    CPMFAccuDataPtr inicf  = InAccu->GetSectionData("TICF");

// output data
    CPMFAccuDataPtr outnsamples    = OutAccu->CreateSectionData("TDS_NSAMPLES",  "AD","R","B");
    CPMFAccuDataPtr outmpp         = OutAccu->CreateSectionData("MTDS_PP",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pp        = OutAccu->CreateSectionData("M2TDS_PP",      "M2","R","M","TDS_NSAMPLES","MTDS_PP");
    CPMFAccuDataPtr outmpn         = OutAccu->CreateSectionData("MTDS_PN",       "WA","R","M","TDS_NSAMPLES");
    CPMFAccuDataPtr outm2pn        = OutAccu->CreateSectionData("M2TDS_PN",      "M2","R","M","TDS_NSAMPLES","MTDS_PN");

    for(size_t t=10; t < NSTLimit; t++){

        incvs->GetData(t,ipos);
        CalculateKSWeights(ipos,weights);

        double epot = inepot->GetData(t);
        double erst = inerst->GetData(t);
        double ekin = inekin->GetData(t);

        for(size_t ibin=0; ibin < NumOfBins; ibin++){
        // increase number of samples
            double w    = weights[ibin];
            double n    = outnsamples->GetData(ibin);
            n = n + w;
            outnsamples->SetData(ibin,n);

            if( w <= 0.0 ) continue;

            double inv = w / n;

            for(size_t icv=0; icv < NumOfCVs; icv++){
                double icf  = inicf->GetData(t,icv);
                double pp   = epot+erst+ekin+icf;
                double pn   = epot+erst+ekin-icf;

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
            double dx = (InAccu->GetCV(icv)->GetDifference(ipos[icv],jpos[icv])) / (KSWfac[icv] * InAccu->GetCV(icv)->GetBinWidth());
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
    if( u2 < 1.0 ){
        return(35.0/32.0*(1-u2)*(1-u2)*(1-u2));
    }

    return(0.0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEntropyDer::Finalize(void)
{
    // close files if they are own by program
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# pmf-entropy terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

