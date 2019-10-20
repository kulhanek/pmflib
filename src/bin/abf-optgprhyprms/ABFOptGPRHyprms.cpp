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

#include <math.h>
#include <errno.h>
#include <ErrorSystem.hpp>
#include <EnergySurface.hpp>
#include "ABFOptGPRHyprms.hpp"
#include <iomanip>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFOptGPRHyprms)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFOptGPRHyprms::CABFOptGPRHyprms(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFOptGPRHyprms::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::high);
    } else {
        vout.Verbosity(CVerboseStr::low);
    }

    StartTime.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-optgprhyprms (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if( Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accu file (in)        : " << Options.GetArgABFAccuName()  << endl;
    } else {
        vout << "# ABF accu file (in)        : - (standard input)" << endl;
    }
    if( Options.GetArgGPRHyprmsName() != "-") {
        vout << "# Optimized GPR hyprms (out): " << Options.GetArgGPRHyprmsName() << endl;
    } else {
        vout << "# Optimized GPR hyprms (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if( Options.IsOptLoadHyprmsSet() ){
        vout << "# GPR hyperprms file    : " << Options.GetOptLoadHyprms() << endl;
    } else {
        vout << "# SigmaF2               : " << setprecision(3) << Options.GetOptSigmaF2() << endl;
        vout << "# NCorr                 : " << setprecision(3) << Options.GetOptNCorr() << endl;
        vout << "# Width factor wfac     : " << (const char*)Options.GetOptWFac() << endl;
    }
    vout << "# Split NCorr mode      : " << bool_to_str(Options.GetOptSplitNCorr()) << endl;
    vout << "# ------------------------------------------------" << endl;
    vout << "# Linear algebra        : " << Options.GetOptLAMethod() << endl;
    if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2") ){
    vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
    }
    vout << "# ------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }

    if( OutputFile.Open(Options.GetArgGPRHyprmsName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFOptGPRHyprms::Run(void)
{
    State = 1;

// load accumulator
    vout << endl;
    vout << format("%02d:Loading ABF accumulator: %s")%State%string(Options.GetArgABFAccuName()) << endl;
    State++;
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done" << endl;

// print summary
    // print CVS info
    Accumulator.PrintCVSInfo(vout);
    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accumulator.SetNCorr(Options.GetOptNCorr());
    FES.Allocate(&Accumulator);

    vout << endl;
    vout << format("%02d:ABF accumulator statistics ...")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done" << endl;

// run optimization
    InitOptimizer();
    if( ! Options.GetOptTest() ){
        Optimize();
    } else {
        Test();
    }

    vout << endl;
    vout << format("%02d:Saving optimized hyperparameters: %s")%State%string(Options.GetArgGPRHyprmsName()) << endl;
    State++;
    if( WriteHyperPrms(OutputFile) == false ){
        ES_ERROR("unable to save optimized hyperparameters");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::InitOptimizer(void)
{
    vout << endl;
    vout << format("%02d:Optimization of GPR hyperparameters ...")%State << endl;
    State++;
    CABFIntegratorGPR::PrintExecInfo(vout);

// what should be optimized?
    SplitNCorr      = Options.GetOptSplitNCorr();
    SigmaF2Enabled  = Options.GetOptSigmaF2Enabled();
    DecodeEList(Options.GetOptNCorrEnabled(),NCorrEnabled,"--enablencorr");
    DecodeEList(Options.GetOptWFacEnabled(),WFacEnabled,"--enablewfac");

// number of optimized parameters
    NumOfOptPrms = 0;
    NumOfPrms = 0;
// sigmaf2
    if( SigmaF2Enabled ) NumOfOptPrms++;
    NumOfPrms++;
// ncorr
    if( SplitNCorr ){
        for(size_t i=0; i < NCorrEnabled.size(); i++){
            if( NCorrEnabled[i] ) NumOfOptPrms++;
        }
        NumOfPrms += Accumulator.GetNumberOfCoords();
    } else {
        if( NCorrEnabled[0] ) NumOfOptPrms++;
        NumOfPrms++;
    }
// wfac
    for(size_t i=0; i < WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) NumOfOptPrms++;
    }
    NumOfPrms += Accumulator.GetNumberOfCoords();

    if( NumOfOptPrms == 0 ){
        RUNTIME_ERROR("no hyperparameter to optimize");
    }

    int     noptsteps   = Options.GetOptNOptSteps();
    double  termeps     = Options.GetOptTermEps();
    int     nlbfgscorr  = Options.GetOptNumOfLBFGSCorr();

// print header
    vout << "   Maximum number of L-BFGS steps  = " << noptsteps << endl;
    vout << "   Termination criteria            = " << std::scientific << termeps << std::fixed << endl;
    vout << "   Number of L-BFGS corrections    = " << nlbfgscorr << endl;
    vout << "   Total number of hyprms          = " << NumOfPrms << endl;
    vout << "   Number of optimized hyprms      = " << NumOfOptPrms << endl;

// initial values
    if( Options.IsOptLoadHyprmsSet() ){
        LoadGPRHyprms();
    } else {
        SigmaF2 = Options.GetOptSigmaF2();
        DecodeVList(Options.GetOptNCorr(),NCorr,"--ncorr",1.0);
        DecodeVList(Options.GetOptWFac(),WFac,"--wfac",3.0);
    }

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
        vout << format("      SigmaF2  = %10.4f")%SigmaF2 << endl;
    if( SplitNCorr ){
        for(int k=0; k < Accumulator.GetNumberOfCoords(); k++ ){
            vout << format("      NCorr#%-2d = %10.4f")%(k+1)%NCorr[k] << endl;
        }
    } else {
        vout << format("      NCorr    = %10.4f")%NCorr[0] << endl;
    }

    for(int k=0; k < Accumulator.GetNumberOfCoords(); k++ ){
        vout << format("      WFac#%-2d  = %10.4f")%(k+1)%WFac[k] << endl;
    }

    if( ! Options.GetOptTest() ){
        vout << endl;
        vout << "# step         logML";
        if( SigmaF2Enabled )    vout << "    SigmaF2";
        if( SplitNCorr ){
            for(size_t i=0; i < NCorrEnabled.size(); i++){
                if( NCorrEnabled[i] ) vout << format("    NCorr#%-1d")%(i+1);
            }
        } else {
        if( NCorrEnabled[0] )      vout << "      NCorr";
        }
        for(size_t i=0; i < WFacEnabled.size(); i++){
            if( WFacEnabled[i] ) vout << format("     Wfac#%-1d")%(i+1);
        }
        vout << endl;

        vout << "# ---- -------------";
        if( SigmaF2Enabled )    vout << " ----------";
        if( SplitNCorr ){
            for(size_t i=0; i < NCorrEnabled.size(); i++){
                if( NCorrEnabled[i] ) vout << " ----------";
            }
        } else {
            if( NCorrEnabled[0] )      vout << " ----------";
        }
        for(size_t i=0; i < WFacEnabled.size(); i++){
            if( WFacEnabled[i] ) vout << " ----------";
        }
        vout << endl;
    }

    Hyprms.CreateVector(NumOfOptPrms);
    HyprmsGrd.CreateVector(NumOfOptPrms);


// set initial hyprms
    int ind = 0;
    if( SigmaF2Enabled ){
        Hyprms[ind] = SigmaF2;
        ind++;
    }
    if( SplitNCorr ){
        for(int i=0; i < (int)NCorrEnabled.size(); i++){
            if( NCorrEnabled[i] ){
                Hyprms[ind] = NCorr[i];
                ind++;
            }
        }
    } else {
        if( NCorrEnabled[0] ){
            Hyprms[ind] = NCorr[0];
            ind++;
        }
    }

    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ){
            Hyprms[ind] = WFac[i];
            ind++;
        }
    }

// init L-BFGS
    // N*(2*M+1)+2*M
    Work.CreateVector(NumOfOptPrms*(2*nlbfgscorr+1)+2*nlbfgscorr);
    Work.SetZero();
    TmpXG.CreateVector(NumOfOptPrms);
    TmpXG.SetZero();
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname)
{
    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( (int)slist.size() > Accumulator.GetNumberOfCoords() ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << Accumulator.GetNumberOfCoords() << ")";
        RUNTIME_ERROR(error);
    }

    elist.resize(Accumulator.GetNumberOfCoords());

    // parse values
    bool last_st = false;
    for(int i=0; i < (int)slist.size(); i++){
        stringstream str(slist[i]);
        char letter;
        str >> letter;
        if( ! str ){
            CSmallString error;
            error << "unable to decode value for " << optionname << " at position: " << i+1;
            RUNTIME_ERROR(error);
        }
        if( (letter == 'T') || (letter == 't') ){
            last_st = true;
        } else {
            last_st = false;
        }
        elist[i] = last_st;
    }

    // pad the rest with the last value
    for(int i=slist.size(); i < Accumulator.GetNumberOfCoords(); i++){
        elist[i] = last_st;
    }
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::DecodeVList(const CSmallString& spec, CSimpleVector<double>& vlist,const CSmallString& optionname,double defv)
{
    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( (int)slist.size() > Accumulator.GetNumberOfCoords() ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << Accumulator.GetNumberOfCoords() << ")";
        RUNTIME_ERROR(error);
    }

    vlist.CreateVector(Accumulator.GetNumberOfCoords());

    // parse values
    double last_st = defv;
    for(int i=0; i < (int)slist.size(); i++){
        stringstream str(slist[i]);
        str >> last_st;
        if( ! str ){
            CSmallString error;
            error << "unable to decode value for " << optionname << " at position: " << i+1;
            RUNTIME_ERROR(error);
        }
        vlist[i] = last_st;
    }

    // pad the rest with the last value
    for(int i=slist.size(); i < Accumulator.GetNumberOfCoords(); i++){
        vlist[i] = last_st;
    }
}

//------------------------------------------------------------------------------

//        SUBROUTINE LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
//  C
//        INTEGER N,M,IPRINT(2),IFLAG
//        DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
//        DOUBLE PRECISION F,EPS,XTOL
//        LOGICAL DIAGCO

extern "C" void lbfgs_(int* N,int* M,double* X,double* F,double* G,int* DIAGCO,double* DIAG,
                       int* IPRINT,double* EPS,double* XTOL,double* W,int* IFLAG);

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::Optimize(void)
{
    int iflag;
    int diacgo = 0;
    int iprint[2];
    iprint[0] = -1;
    iprint[1] = 1;

    double xtol = 1e-15;
    iflag = 0;

    int     noptsteps   = Options.GetOptNOptSteps();
    double  termeps     = Options.GetOptTermEps();
    int     nlbfgscorr  = Options.GetOptNumOfLBFGSCorr();

    for(int istep=1; istep <= noptsteps; istep++ ){
        if( Options.GetOptNumeric() ){
            RunGPRNumerical();
        } else {
            RunGPRAnalytical();
        }
        WriteResults(istep);

        // we need to maximize logML, thus reverse curvature
        double rv = -logML;
        for(int i=0; i < NumOfOptPrms; i++){
            HyprmsGrd[i] = -HyprmsGrd[i];
        }

        // run L-BFGS
        lbfgs_(&NumOfOptPrms,&nlbfgscorr,Hyprms,&rv,HyprmsGrd,&diacgo,TmpXG,iprint,&termeps,&xtol,Work,&iflag);

        if( iflag == 0 ) break;
        if( iflag < 0 ) {
            CSmallString error;
            error << ">>> ERROR: Internal L-BFGS driver error! Code = " << iflag;
            RUNTIME_ERROR(error);
        }
    }

}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::Test(void)
{
    vout << endl;
    vout << format("%02d:Testing gradients ...")%State << endl;
    State++;

    CSimpleVector<double>   grd1;
    CSimpleVector<double>   grd2;

    vout << "   Analytical gradients ... " <<  flush;
    RunGPRAnalytical();
    grd1 = HyprmsGrd;
    vout << "done." << endl;

    vout << "   Numerical gradients ";
    if( Options.GetOptCD5() ){
        vout << "(CD: 5-points) ... " << flush;
    } else {
        vout << "(CD: 3-points) ... " << flush;
    }
    RunGPRNumerical();
    grd2 = HyprmsGrd;
    vout << "done." << endl;

    vout << endl;
    vout << "# idx       prms     analytical      numerical           diff" << endl;
    vout << "# --- ---------- -------------- -------------- --------------" << endl;

    double trh = 1e-5;

    int ind = 0;

    int soffset = 0;
    int noffset;
    int woffset;

    if( SplitNCorr ){
        noffset = soffset+1;
        woffset = noffset+NCorr.GetLength();
    } else {
        noffset = soffset+1;
        woffset = noffset+1;
    }

    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( HyprmsEnabled[prm] ){
            double diff = grd1[ind]-grd2[ind];
            string rel;
            if( fabs(grd1[ind]) > 1.0 ){
                // switch to relative scale
                diff /= grd2[ind];
                rel = " rel.err.";
            }
            if( fabs(diff) > trh ){
                vout << "<red><bold>";
            }
            vout << format("%5d ")%(ind+1);
            if( prm == 0 ){
                vout << "SigmaF2    ";
            } else if( (prm >= noffset) && (prm < woffset) ){
                if( SplitNCorr ){
               int cv = prm - noffset;
                vout << format("NCorr#%-2d   ")%(cv+1);
                } else {
                vout << "NCorr      ";
                }
            } else {
                int cv = prm - woffset;
                vout << format("WFac#%-2d    ")%(cv+1);
            }

            vout << format("%14.6e %14.6e %14.6e")%grd1[ind]%grd2[ind]%diff << rel << endl;
            if( fabs(diff) > trh ){
                vout << "</bold></red>";
            }
            ind++;
        }
    }
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::RunGPRAnalytical(void)
{
// setup parameters
    ScatterHyprms(Hyprms);

// calculate gradient and logML
    CABFIntegratorGPR gpr;
    logML = GetLogML(gpr);
    gpr.GetLogMLDerivatives(HyprmsEnabled,HyprmsGrd);
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::RunGPRNumerical(void)
{
    CSimpleVector<double>   tmp_prms;
    tmp_prms.CreateVector(NumOfOptPrms);

// value
    tmp_prms = Hyprms;
    ScatterHyprms(tmp_prms);
    {
        CABFIntegratorGPR gpr;
        logML = GetLogML(gpr);
    }

// derivatives
    double dh = 1e-3;

    for(int i=0; i < NumOfOptPrms; i++){
        if( Options.GetOptCD5() ){
            double v1,v2,v3,v4;
            tmp_prms = Hyprms;
            tmp_prms[i] -= 2.0*dh;
            ScatterHyprms(tmp_prms);
            {
                CABFIntegratorGPR gpr;
                v1 = GetLogML(gpr);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CABFIntegratorGPR gpr;
                v2 = GetLogML(gpr);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CABFIntegratorGPR gpr;
                v3 = GetLogML(gpr);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += 2.0*dh;
            ScatterHyprms(tmp_prms);
            {
                CABFIntegratorGPR gpr;
                v4 = GetLogML(gpr);
            }
            HyprmsGrd[i] = (v1-8.0*v2+8.0*v3-v4)/(12.0*dh);
        } else {
            double lv,rv;
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CABFIntegratorGPR gpr;
                lv = GetLogML(gpr);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CABFIntegratorGPR gpr;
                rv = GetLogML(gpr);
            }
            HyprmsGrd[i] = (rv-lv)/(2.0*dh);
        }

    }
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::ScatterHyprms(CSimpleVector<double>& hyprsm)
{
    // update parameters
    HyprmsEnabled.resize(NumOfPrms);

    int ind = 0;
    int i = 0;
// ---------------
    if( SigmaF2Enabled ){
        SigmaF2 = fabs(hyprsm[ind]);
        ind++;
        HyprmsEnabled[i] = true;
    } else {
        HyprmsEnabled[i] = false;
    }
    i++;
// ---------------
    if( SplitNCorr ){
        for(int k=0; k < (int)NCorrEnabled.size(); k++){
            if( NCorrEnabled[k] ){
                NCorr[k] = fabs(hyprsm[ind]);
                ind++;
                HyprmsEnabled[i] = true;
            } else {
                HyprmsEnabled[i] = false;
            }
            i++;
        }
    } else {
        if( NCorrEnabled[0] ){
            NCorr[0] = fabs(hyprsm[ind]);
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
// ---------------
    for(int k=0; k < (int)WFacEnabled.size(); k++){
        if( WFacEnabled[k] ){
            WFac[k] = fabs(hyprsm[ind]);
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
}

//------------------------------------------------------------------------------

double CABFOptGPRHyprms::GetLogML(CABFIntegratorGPR& gpr)
{
    gpr.SetInputABFAccumulator(&Accumulator);
    gpr.SetOutputFESurface(&FES);

    gpr.SetRCond(Options.GetOptRCond());

    gpr.SetIncludeError(false);
    gpr.SetNoEnergy(true);
    gpr.IncludeGluedAreas(false);

    gpr.SetSplitNCorr(SplitNCorr);

    if( Options.GetOptLAMethod() == "svd" ){
        gpr.SetINVMehod(EGPRINV_SVD);
    } else if( Options.GetOptLAMethod() == "svd2" ){
            gpr.SetINVMehod(EGPRINV_SVD2);
    } else if( Options.GetOptLAMethod() == "lu" ) {
        gpr.SetINVMehod(EGPRINV_LU);
    } else if( Options.GetOptLAMethod() == "ll" ) {
        gpr.SetINVMehod(EGPRINV_LL);
    } else if( Options.GetOptLAMethod() == "default" ) {
        // nothing to do - use default method set in constructor of integrator
    } else {
        INVALID_ARGUMENT("algorithm - not implemented");
    }

// setup integrator
    gpr.SetSigmaF2(SigmaF2);
    gpr.SetNCorr(NCorr);
    gpr.SetWFac(WFac);
    vout << high;
    gpr.Integrate(vout,true);
    vout << low;

// calculate logML
    return(gpr.GetLogMarginalLikelihood());
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::WriteResults(int istep)
{
    vout << format("%6d %13e")%istep%logML;
    if( SigmaF2Enabled ) vout << format(" %10.3f")%SigmaF2;
    if( SplitNCorr ){
        for(int i=0; i < (int)NCorrEnabled.size(); i++){
            if( NCorrEnabled[i] ) vout << format(" %10.3f")%NCorr[i];
        }
    } else {
        if( NCorrEnabled[0] )   vout << format(" %10.3f")%NCorr[0];
    }
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) vout << format(" %10.3f")%WFac[i];
    }
    vout << endl;
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = Accumulator.GetNumberOfBins();
    int    sampled = 0;
    int    glued = 0;
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( Accumulator.GetNumberOfABFSamples(ibin) > 0 ) {
            sampled++;
        }
        if( Accumulator.GetNumberOfABFSamples(ibin) < 0 ) {
            glued++;
        }
    }
    if( maxbins > 0 ){
        vout << "   Sampled area: "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
    if( glued > 0 ){
        vout << "   Glued area:   "
             << setw(6) << glued << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << glued/maxbins*100 <<"%" << endl;
    }
}

//------------------------------------------------------------------------------

bool CABFOptGPRHyprms::WriteHyperPrms(FILE* p_fout)
{
    if( p_fout == NULL ){
        return(false);
    }

    ScatterHyprms(Hyprms);

    if( fprintf(p_fout,"# GPR hyperparameters for abf-integrate\n") <= 0 ) return(false);
    if( fprintf(p_fout,"SigmaF2  = %10.4f\n",SigmaF2) <= 0 ) return(false);
    if( SplitNCorr ){
        for(int i=0; i < NCorr.GetLength(); i++ ){
            if( fprintf(p_fout,"NCorr#%-2d = %10.4f\n",i+1,NCorr[i]) <= 0 ) return(false);
        }
    } else {
        if( fprintf(p_fout,"NCorr    = %10.4f\n",NCorr[0]) <= 0 ) return(false);
    }
    for(int i=0; i < WFac.GetLength(); i++ ){
        if( fprintf(p_fout,"WFac#%-2d  = %10.4f\n",i+1,WFac[i]) <= 0 ) return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFOptGPRHyprms::LoadGPRHyprms(void)
{
    ifstream fin;
    fin.open(Options.GetOptLoadHyprms());
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << Options.GetOptLoadHyprms();
        RUNTIME_ERROR(error);
    }

    int ncvs = Accumulator.GetNumberOfCoords();
    NCorr.CreateVector(ncvs);
    WFac.CreateVector(ncvs);

    string line;
    while( getline(fin,line) ){
        // is it comment?
        if( (line.size() > 0) && (line[0] == '#') ) continue;

        // parse line
        stringstream str(line);
        string key, buf;
        double value;
        str >> key >> buf >> value;
        if( ! str ){
            CSmallString error;
            error << "GPR hyperparameters file, unable to decode line: " << line.c_str();
            RUNTIME_ERROR(error);
        }
        if( key == "SigmaF2" ){
            SigmaF2 = value;
        } else if( key == "NCorr" ){
            NCorr[0] = value;
            if( Options.GetOptSplitNCorr() ){
                vout << "   >> WARNING: Expanding NCorr due to --splitncorr" << endl;
                for(int i=0; i < ncvs;i++){
                    NCorr[1] = value;
                }
            }
        } else if( key.find("NCorr#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string sncorr;
            int    cvind;
            kstr >> sncorr >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode ncorr key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( (cvind < 0) || (cvind >= ncvs) ){
                CSmallString error;
                error << "ncorr index " << (cvind+1) << " out-of-range 1-" << ncvs;
                RUNTIME_ERROR(error);
            }
            if( Options.GetOptSplitNCorr() ){
                NCorr[cvind] = value;
            } else {
                if( cvind > 0 ){
                    vout << "   >> WARNING: NCorr with index higher than one ignored (no --splitncorr)!" << endl;
                } else {
                    NCorr[cvind] = value;
                }
            }

        } else if( key.find("WFac#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode wfac key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( (cvind < 0) || (cvind >= ncvs) ){
                CSmallString error;
                error << "wfac index " << (cvind+1) << " out-of-range 1-" << ncvs;
                RUNTIME_ERROR(error);
            }
            WFac[cvind] = value;
        } else {
            CSmallString error;
            error << "GPR hyperparameters file, unrecognized key: " << key.c_str();
            RUNTIME_ERROR(error);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFOptGPRHyprms::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    CSmallTime dur;
    dur = dt - StartTime;

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-optgprhyprms terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

