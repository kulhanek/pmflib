// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "OptGPRHyprms.hpp"
#include <iomanip>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <SciLapack.hpp>
// ----------
#include <IntegratorGPR.hpp>
#include <SmootherGPR.hpp>
#include <GHSIntegratorGPR0A.hpp>
#include <GHSIntegratorGPRcA.hpp>
#include <GHSIntegratorGPR0B.hpp>
// ----------
#include <ABFProxy_dG.hpp>
#include <ABFProxy_mTdS.hpp>
#include <CSTProxy_dG.hpp>
#include <CSTProxy_mTdS.hpp>
#include <PMFProxy_dH.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(COptGPRHyprms)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

COptGPRHyprms::COptGPRHyprms(void)
{
    NCVs = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int COptGPRHyprms::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CPMFIntOpts
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
    vout << "# gpr-opthyprms (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

        vout << "# PMF Accumulator       : " << Options.GetArgAccuFile() << endl;
        vout << "# Processed realm       : " << Options.GetArgRealm() << endl;
    if( Options.GetArgHyprmsFile() != "-") {
        vout << "# Optimized GPR hyprms  : " << Options.GetArgHyprmsFile() << endl;
    } else {
        vout << "# Optimized GPR hyprms  : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
    if( Options.IsOptLoadHyprmsSet() ){
        vout << "# GPR hyperprms file    : " << Options.GetOptLoadHyprms() << endl;
    } else {
        vout << "# SigmaF2               : " << (const char*)Options.GetOptSigmaF2() << endl;
        vout << "# CoVar                 : " << (const char*)Options.GetOptCoVar() << endl;
        vout << "# Width factor wfac     : " << (const char*)Options.GetOptWFac() << endl;
        vout << "# NCorr                 : " << (const char*)Options.GetOptNCorr() << endl;
        vout << "# SigmaN2               : " << (const char*)Options.GetOptSigmaN2() << endl;
    }
    vout << "# ------------------------------------------------" << endl;
        vout << "# Linear algebra        : " << Options.GetOptLAMethod() << endl;
    if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2") || (Options.GetOptLAMethod() == "default") ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
    }
        vout << "# Optimized target      : " << Options.GetOptTarget() << endl;

    // here we assume that kernels in IntegratorGPR and SmootherGPR are the same
    CSmootherGPR gpr;
    gpr.SetKernel(Options.GetOptGPRKernel());
        vout << "# GPR Kernel            : " << gpr.GetKernelName() << endl;
    vout << "# ------------------------------------------------" << endl;

    if( Options.GetOptTarget() == "logml" ){
        Target = EGOT_LOGML;
    } else if ( Options.GetOptTarget() == "logpl" ){
        Target = EGOT_LOGPL;
    } else {
        RUNTIME_ERROR("unsupported target");
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool COptGPRHyprms::Run(void)
{
    State = 1;

// load accumulator
    vout << endl;
    vout << format("%02d:Loading PMF accumulator: %s")%State%string(Options.GetArgAccuFile()) << endl;
    State++;
    Accu = CPMFAccumulatorPtr(new CPMFAccumulator);
    try {
        Accu->Load(Options.GetArgAccuFile());
    } catch(...) {
        CSmallString error;
        error << "unable to load the input PMF accumulator file '" << Options.GetArgAccuFile() << "'";
        ES_ERROR(error);
        return(false);
    }
    FES = CEnergySurfacePtr(new CEnergySurface);
    FES->Allocate(Accu);
    HES = CEnergySurfacePtr(new CEnergySurface);
    HES->Allocate(Accu);
    SES = CEnergySurfacePtr(new CEnergySurface);
    SES->Allocate(Accu);

    if( Options.IsOptGlobalMinSet() ){
        FES->SetGlobalMin(Options.GetOptGlobalMin());
    }
    vout << "   Done" << endl;

// statistics
    vout << endl;
    vout << format("%02d:PMF accumulator statistics ...")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done" << endl;

// realms
    vout << endl;
    vout << format("%02d:Initializing realm: %s")%State%string(Options.GetArgRealm()) << endl;
    State++;
    InitGPREngine();
    vout << "   Done" << endl;

// run optimization
    bool result = true;
    InitOptimizer();
    if( Options.GetOptTest() ){
        Test();
    }
    result = Optimize();
    if( Options.GetOptTest() ){
        Test();
    }

    if( Options.GetOptSPType() ){
        CalcHessian();
    }

    if( Options.GetOptPrintStat() ){
        ShowGPRStat();
    }

    CSmallString output = Options.GetProgArg(Options.GetNumberOfProgArgs()-1);
    vout << endl;
    vout << format("%02d:Saving optimized hyperparameters: %s")%State%string(output) << endl;
    State++;
    if( OutputFile.Open(output,"w") == false ){
        ES_ERROR("unable to open output file");
        return(false);
    }
    if( WriteHyperPrms(OutputFile) == false ){
        ES_ERROR("unable to save optimized hyperparameters");
        return(false);
    }
    if( result == false ){
        // indicate that hyperparameters are not at stationary point
        fprintf(OutputFile,"# not stationary point\n");
    }

    return(result);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitOptimizer(void)
{
    vout << endl;
    vout << format("%02d:Optimization of GPR hyperparameters ...")%State << endl;
    State++;
    CIntegratorGPR gpr;
    gpr.PrintExecInfo(vout);

// what should be optimized?
    DecodeEList(Options.GetOptSigmaF2Enabled(),SigmaF2,SigmaF2Enabled,"--enablesigmaf2");
    DecodeEList(Options.GetOptCoVarEnabled(),CoVar,CoVarEnabled,"--enablecovar");
    DecodeEList(Options.GetOptWFacEnabled(),WFac,WFacEnabled,"--enablewfac");
    DecodeEList(Options.GetOptNCorrEnabled(),NCorr,NCorrEnabled,"--enablencorr");
    DecodeEList(Options.GetOptSigmaN2Enabled(),SigmaN2,SigmaN2Enabled,"--enablesigman2");

// number of optimized parameters
    NumOfOptPrms = 0;
    NumOfPrms = 0;
// sigmaf2
    for(size_t i=0; i < SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ) NumOfOptPrms++;
        NumOfPrms++;
    }
// covar
    for(size_t i=0; i < CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] ) NumOfOptPrms++;
        NumOfPrms++;
    }
// wfac
    for(size_t i=0; i < WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) NumOfOptPrms++;
        NumOfPrms++;
    }
// ncorr
    for(size_t i=0; i < NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] ) NumOfOptPrms++;
        NumOfPrms++;
    }
// sigman2
    for(size_t i=0; i < SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) NumOfOptPrms++;
        NumOfPrms++;
    }

    if( NumOfOptPrms == 0 ){
        RUNTIME_ERROR("no hyperparameter to optimize");
    }

    int     noptsteps   = Options.GetOptNOptSteps();
    double  termeps     = Options.GetOptTermEps();
    double  termval     = Options.GetOptTermVal();
    int     nlbfgscorr  = Options.GetOptNumOfLBFGSCorr();

// print header
    vout << "   Maximum number of L-BFGS steps  = " << noptsteps << endl;
    vout << "   Max number of optimizer resets  = " << Options.GetOptNumOfResets() << endl;
    vout << "   Termination criteria (L-BFGS)   = " << std::scientific << termeps << std::fixed << endl;
    vout << "   Termination criteria (value)    = " << std::scientific << termval << std::fixed << endl;
    vout << "   Number of L-BFGS corrections    = " << nlbfgscorr << endl;
    vout << "   Total number of hyprms          = " << NumOfPrms << endl;
    vout << "   Number of optimized hyprms      = " << NumOfOptPrms << endl;

// initial values
    if( Options.IsOptLoadHyprmsSet() ){
        LoadGPRHyprms();
    } else {
        DecodeVList(Options.GetOptSigmaF2(),SigmaF2,"--sigmaf2",15.0);
        for(size_t i=0; i < SigmaF2.GetLength(); i++){
            if( SigmaF2[i] < Options.GetOptMinSigmaF2() ){
                RUNTIME_ERROR("--sigmaf2 has to be greater than or equal to --minsigmaf2");
            }
        }
        DecodeVList(Options.GetOptCoVar(),CoVar,"--covar",0.0);
        for(size_t i=0; i < CoVar.GetLength(); i++){
            if( CoVar[i] < Options.GetOptMinCoVar() ){
                RUNTIME_ERROR("--covar has to be greater than or equal to --mincovar");
            }
        }
        DecodeVList(Options.GetOptWFac(),WFac,"--wfac",3.0);
        for(size_t i=0; i < WFac.GetLength(); i++){
            if( WFac[i] < Options.GetOptMinWFac() ){
                RUNTIME_ERROR("--wfac has to be greater than or equal to --minwfac");
            }
        }
        DecodeVList(Options.GetOptNCorr(),NCorr,"--ncorr",0.0);
        for(size_t i=0; i < NCorr.GetLength(); i++){
            if( NCorr[i] < Options.GetOptMinNCorr() ){
                RUNTIME_ERROR("--ncorr has to be greater than or equal to --minncorr");
            }
        }
        DecodeVList(Options.GetOptSigmaN2(),SigmaN2,"--sigman2",0.0);
        for(size_t i=0; i < SigmaN2.GetLength(); i++){
            if( SigmaN2[i] < Options.GetOptMinSigmaN2() ){
                RUNTIME_ERROR("--sigman2 has to be greater than or equal to --minsigman2");
            }
        }
    }

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
    for(int k=0; k < (int)SigmaF2.GetLength(); k++ ){
        vout << format("      SigmaF2#%-2d = %10.4f")%(k+1)%SigmaF2[k] << endl;
    }
    for(int k=0; k < (int)CoVar.GetLength(); k++ ){
        vout << format("      CoVar#%-2d   = %10.4f")%(k+1)%CoVar[k] << endl;
    }
    for(int k=0; k < (int)WFac.GetLength(); k++ ){
        vout << format("      WFac#%-2d    = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(int k=0; k < (int)NCorr.GetLength(); k++ ){
        vout << format("      NCorr#%-2d   = %10.4e")%(k+1)%NCorr[k] << endl;
    }
    for(int k=0; k < (int)SigmaN2.GetLength(); k++ ){
        vout << format("      SigmaN2#%-2d = %10.4e")%(k+1)%SigmaN2[k] << endl;
    }

    Hyprms.CreateVector(NumOfOptPrms);
    HyprmsGrd.CreateVector(NumOfOptPrms);

// set initial hyprms
    int ind = 0;
    for(int i=0; i < (int)SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ){
            Hyprms[ind] = sqrt(SigmaF2[i]-Options.GetOptMinSigmaF2());
            ind++;
        }
    }
    for(int i=0; i < (int)CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] ){
            Hyprms[ind] = sqrt(CoVar[i]-Options.GetOptMinCoVar());
            ind++;
        }
    }
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ){
            Hyprms[ind] = sqrt(WFac[i]-Options.GetOptMinWFac());
            ind++;
        }
    }
    for(int i=0; i < (int)NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] ){
            Hyprms[ind] = sqrt(NCorr[i]-Options.GetOptMinNCorr());
            ind++;
        }
    }
    for(int i=0; i < (int)SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ){
            // Hyprms[ind] = sqrt(SigmaN2[i]-Options.GetOptMinSigmaN2());
            Hyprms[ind] = log(SigmaN2[i]);
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

void COptGPRHyprms::DecodeEList(const CSmallString& spec, CSimpleVector<double>& vlist, std::vector<bool>& elist,const CSmallString& optionname)
{
    if( vlist.GetLength() == 0 ) return;

    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( slist.size() > vlist.GetLength() ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << vlist.GetLength() << ")";
        RUNTIME_ERROR(error);
    }

    elist.resize(vlist.GetLength());

    // parse values
    bool last_st = false;
    for(size_t i=0; i < slist.size(); i++){
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
    for(size_t i=slist.size(); i < vlist.GetLength(); i++){
        elist[i] = last_st;
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::DecodeVList(const CSmallString& spec, CSimpleVector<double>& vlist,const CSmallString& optionname,double defv)
{
    if( vlist.GetLength() == 0 ) return;

    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( slist.size() > vlist.GetLength() ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << vlist.GetLength() << ")";
        RUNTIME_ERROR(error);
    }

    // parse values
    double last_st = defv;
    for(size_t i=0; i < slist.size(); i++){
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
    for(size_t i=slist.size(); i < vlist.GetLength(); i++){
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

#ifdef HAVE_MKL_ILP64
// 64-bit integer - ILP64
typedef long int FT_INT;
#else
// standard 32-bit integer
typedef int FT_INT;
#endif

extern "C" void lbfgs_(FT_INT* N,FT_INT* M,double* X,double* F,double* G,FT_INT* DIAGCO,double* DIAG,
                       FT_INT* IPRINT,double* EPS,double* XTOL,double* W,FT_INT* IFLAG);

//------------------------------------------------------------------------------

bool COptGPRHyprms::Optimize(void)
{
    bool result = true;

    FT_INT iflag;
    FT_INT diacgo = 0;
    FT_INT iprint[2];
    iprint[0] = -1;
    iprint[1] = 1;

    vout << endl;
    switch(Target){
        case(EGOT_LOGML):
            vout << "# step         logML";
        break;
        case(EGOT_LOGPL):
            vout << "# step         logPL";
        break;
    }

    for(size_t i=0; i < SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ) vout << format(" SigmaF2#%-2d")%(i+1);
    }
    for(size_t i=0; i < CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] )   vout << format("   CoVar#%-2d")%(i+1);
    }
    for(size_t i=0; i < WFacEnabled.size(); i++){
        if( WFacEnabled[i] )    vout << format("    WFac#%-2d")%(i+1);
    }
    for(size_t i=0; i < NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] )   vout << format("   NCorr#%-2d")%(i+1);
    }
    for(size_t i=0; i < SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) vout << format(" SigmaN2#%-2d")%(i+1);
    }
    vout << endl;

    vout << "# ---- -------------";
    for(size_t i=0; i < SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ) vout << " ----------";
    }
    for(size_t i=0; i < CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] )   vout << " ----------";
    }
    for(size_t i=0; i < WFacEnabled.size(); i++){
        if( WFacEnabled[i] )    vout << " ----------";
    }
    for(size_t i=0; i < NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] )   vout << " ----------";
    }
    for(size_t i=0; i < SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) vout << " ----------";
    }
    vout << endl;


    double xtol = 1e-15;
    iflag = 0;

    int     noptsteps   = Options.GetOptNOptSteps();
    double  termeps     = Options.GetOptTermEps();
    FT_INT  nlbfgscorr  = Options.GetOptNumOfLBFGSCorr();
    int     numofreset  = 0;
    double  last_logtrg = 0;

    int istep = 0;

    CSimpleVector<double>   OldHyprms;
    OldHyprms.CreateVector(NumOfOptPrms);

    // input parameters cannot be zero or negative
    for(int i=0; i < NumOfOptPrms; i++){
        if( Hyprms[i] <= 0.0 ){
    //        Hyprms[i] = 1.0;
        }
    }

    int insupr = 0;
    int nochan = 0;

    for(int k=1; k <= noptsteps; k++ ){
        OldHyprms = Hyprms;

        if( Options.GetOptNumeric() ){
            RunGPRNumerical();
        } else {
            RunGPRAnalytical();
        }
        WriteResults(istep);

        // we need to maximize logML, thus reverse curvature
        double rv = -logTarget;
        for(int i=0; i < NumOfOptPrms; i++){
            HyprmsGrd[i] = -HyprmsGrd[i];
        }

        // run L-BFGS
        FT_INT nx = NumOfOptPrms;
        lbfgs_(&nx,&nlbfgscorr,Hyprms,&rv,HyprmsGrd,&diacgo,TmpXG,iprint,&termeps,&xtol,Work,&iflag);

        istep++;

        if( iflag == 0 ) break;
        if( iflag < 0 ) {
            insupr++;
            if( insupr > 3 ){
                vout << endl;
                vout << "<b><blue>>>> INFO: Insufficient progress, resetting ...</blue></b>" << endl;
                if( ResetOpt(numofreset) == false ){
                    result = false;
                    break;
                }
                vout << endl;
                iflag = 0;
                last_logtrg = logTarget - 10; // be sure that the number is somehow different
                continue;
            }
        } else {
            insupr = 0; // reset counter
        }

        if( istep > 1 ){
            if( fabs(logTarget - last_logtrg) < Options.GetOptTermVal() ){
                nochan++;
                if( nochan > 3 ) {
                    double gnorm = GetGNorm();
                    if( gnorm > termeps*1000 ){
                        vout << endl;
                        vout << "<b><blue>>>> INFO: No significant change, but gradient is large - resetting ...</blue></b>" << endl;
                        vout <<   format("          gnorm = %14.6e")%gnorm << endl;
                        if( ResetOpt(numofreset) == false ){
                            result = false;
                            break;
                        }
                        vout << endl;
                        iflag = 0;
                        last_logtrg = logTarget - 10; // be sure that the number is somehow different
                        continue;
                    } else {
                        vout << endl;
                        vout << "<b><blue>>>> INFO: No significant change - terminating ...</blue></b>" << endl;
                        vout << endl;
                        break;
                    }
                }
            } else {
                nochan = 0;
            }
        }
        last_logtrg = logTarget;
    }

    if( noptsteps > 0 ) {
        vout << "# ---- -------------";
        for(size_t i=0; i < SigmaF2Enabled.size(); i++){
            if( SigmaF2Enabled[i] ) vout << " ----------";
        }
        for(size_t i=0; i < CoVarEnabled.size(); i++){
            if( CoVarEnabled[i] ) vout << " ----------";
        }
        for(size_t i=0; i < WFacEnabled.size(); i++){
            if( WFacEnabled[i] ) vout << " ----------";
        }
        for(size_t i=0; i < NCorrEnabled.size(); i++){
            if( NCorrEnabled[i] ) vout << " ----------";
        }
        for(size_t i=0; i < SigmaN2Enabled.size(); i++){
            if( SigmaN2Enabled[i] ) vout << " ----------";
        }
        vout << endl;
    }

    // final gradients - for the last values of hyprms or SP type calculations
    if( Options.GetOptNumeric() ){
        RunGPRNumerical();
    } else {
        RunGPRAnalytical();
    }
    WriteResults(istep);

    if( (istep >= noptsteps) && (noptsteps > 0) ){
        vout << endl;
        vout << "<red><b>>>> ERROR: Insufficient number of optimization steps ...</b></red>" << endl;
        result = false;
    }

    PrintGradientSummary();

    double gnorm = GetGNorm();
    if( gnorm > termeps*1000 ){
        vout << endl;
        vout << "<b><red>>>> ERROR: Gradient is too large, stationary point was not likely found ...</red></b>" << endl;
        vout <<  format("           gnorm = %14.6e")%gnorm << endl;
        result = false;
    }

    return(result);
}

//------------------------------------------------------------------------------

bool COptGPRHyprms::ResetOpt(int& numofreset)
{
    numofreset++;
    if( numofreset > Options.GetOptNumOfResets() ){
        vout <<  "<b><red>          Too many resets, exiting ...</red></b>" << endl;
        vout << endl;
        return(false);
    }
    if( numofreset > 1 ){
        vout <<          "          Perturbing parameters ..." << endl;
        // perturbing parameters
        for(int i=0; i < NumOfOptPrms; i++){
            float r = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
            Hyprms[i] += 0.5*r + 0.1;    // perturb only up due to down limits
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

double COptGPRHyprms::GetGNorm(void)
{
    double gnorm = 0.0;
    int ind = 0;
    for(size_t prm=0; prm < HyprmsEnabled.size(); prm++){
        if( HyprmsEnabled[prm] ){
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            ind++;
        }
    }
    if( NumOfOptPrms > 0 ){
        gnorm = sqrt(gnorm/(double)NumOfOptPrms);
    }
    return(gnorm);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::Test(void)
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
            vout << setw(10) << left << GetPrmName(prm) << " " << right;
            vout << format("%14.6e %14.6e %14.6e")%grd1[ind]%grd2[ind]%diff << rel << endl;
            if( fabs(diff) > trh ){
                vout << "</bold></red>";
            }
            ind++;
        }
    }
}

//------------------------------------------------------------------------------

std::string COptGPRHyprms::GetPrmName(int prm)
{
    int ind = 0;
    for(size_t i=0; i < SigmaF2Enabled.size(); i++){
        stringstream str;
        str << format("SigmaF2#%-2d")%(i+1);
        if( ind == prm) return(str.str());
        ind++;
    }
    for(size_t i=0; i < CoVarEnabled.size(); i++){
        stringstream str;
        str << format("CoVar#%-2d")%(i+1);
        if( ind == prm) return(str.str());
        ind++;
    }
    for(size_t i=0; i < WFacEnabled.size(); i++){
        stringstream str;
        str << format("WFac#%-2d")%(i+1);
        if( ind == prm) return(str.str());
        ind++;
    }
    for(size_t i=0; i < NCorrEnabled.size(); i++){
        stringstream str;
        str << format("NCorr#%-2d")%(i+1);
        if( ind == prm) return(str.str());
        ind++;
    }
    for(size_t i=0; i < SigmaN2Enabled.size(); i++){
        stringstream str;
        str << format("SigmaN2#%-2d")%(i+1);
        if( ind == prm) return(str.str());
        ind++;
    }
    return("");
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CalcHessian(void)
{
    vout << endl;
    vout << format("%02d:Calculating numerical hessian by central differences ... ")%State << endl;
    State++;

    Hessian.CreateMatrix(NumOfOptPrms,NumOfOptPrms);
    Hessian.SetZero();
    EigenValues.CreateVector(NumOfOptPrms);
    EigenValues.SetZero();

    CSimpleVector<double>   tmp_prms;
    tmp_prms.CreateVector(NumOfOptPrms);

    CSimpleVector< CSimpleVector<double> >   grd1;
    grd1.CreateVector(NumOfOptPrms);
    for(int i=0; i < NumOfOptPrms; i++){
        grd1[i].CreateVector(NumOfOptPrms);
    }

    CSimpleVector< CSimpleVector<double> >   grd2;
    grd2.CreateVector(NumOfOptPrms);
    for(int i=0; i < NumOfOptPrms; i++){
        grd2[i].CreateVector(NumOfOptPrms);
    }

// derivatives
    double dh = 1e-5;

    for(int i=0; i < NumOfOptPrms; i++){

        vout << "   + perturbation for hyprm: " << setw(3) << (i+1) << endl;
        grd1[i].SetZero();
        tmp_prms = Hyprms;
        tmp_prms[i] += dh*Hyprms[i];
        ScatterHyprms(tmp_prms);

        // calculate gradient and logML
        vout << high;
        CreateGPREngine();
        GetTarget();
        GetTargetDerivatives(grd1[i]);
        vout << low;

        vout << "   - perturbation for hyprm: " << setw(3) << (i+1) << endl;
        grd2[i].SetZero();
        tmp_prms = Hyprms;
        tmp_prms[i] -= dh*Hyprms[i];
        ScatterHyprms(tmp_prms);


        // calculate gradient and logML
        vout << high;
        CreateGPREngine();
        GetTarget();
        GetTargetDerivatives(grd2[i]);
        vout << low;
    }

// hessian by central differences from gradients
// https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
    for(int i=0; i < NumOfOptPrms; i++){
        for(int j=0; j < NumOfOptPrms; j++){
            Hessian[i][j] = (grd1[j][i] - grd2[j][i])/(4.0*dh*Hyprms[i]) + (grd1[i][j] - grd2[i][j])/(4.0*dh*Hyprms[i]);
        }
    }

// find eigenvalues
    CSciLapack::syev('V','U',Hessian,EigenValues);

    vout << endl;
    vout << "# idx    hess eigval";

    int ind = 0;
    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( ! HyprmsEnabled[prm] ) continue;
        vout << " " << setw(10) << left << GetPrmName(prm) << right;
    }
    vout << endl;

    vout << "# --- --------------";
    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( HyprmsEnabled[prm] ){
            vout << " ----------";
        }
    }
    vout << endl;

    for(int i=0; i < NumOfOptPrms; i++){
        vout << format("%5d %14.6e")%(i+1)%EigenValues[i];
        ind = 0;
        for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
            if( HyprmsEnabled[prm] ){
                vout << format(" %10.4f")%Hessian[ind][i];
                ind++;
            }
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::PrintGradientSummary(void)
{
    vout << endl;
    vout << "# Final results ..." << endl;
    vout << "# idx       prms          value       gradient" << endl;
    vout << "# --- ---------- -------------- --------------" << endl;

    int ind = 0;

    double gnorm = 0.0;

    for(size_t i=0; i < SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ) {
            vout << format("%5d ")%(ind+1);
            vout << format("SigmaF2#%-2d ")%(i+1);
            double value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinSigmaF2();
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            vout << format("%14.6e %14.6e")%(value)%HyprmsGrd[ind] << endl;
            ind++;
        }
    }
    for(size_t i=0; i < CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] ) {
            vout << format("%5d ")%(ind+1);
            vout << format("CoVar#%-2d   ")%(i+1);
            double value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinCoVar();
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            vout << format("%14.6e %14.6e")%(value)%HyprmsGrd[ind] << endl;
            ind++;
        }
    }
    for(size_t i=0; i < WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) {
            vout << format("%5d ")%(ind+1);
            vout << format("WFac#%-2d    ")%(i+1);
            double value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinWFac();
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            vout << format("%14.6e %14.6e")%(value)%HyprmsGrd[ind] << endl;
            ind++;
        }
    }
    for(size_t i=0; i < NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] ) {
            vout << format("%5d ")%(ind+1);
            vout << format("NCorr#%-2d   ")%(i+1);
            double value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinNCorr();
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            vout << format("%14.6e %14.6e")%(value)%HyprmsGrd[ind] << endl;
            ind++;
        }
    }
    for(size_t i=0; i < SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) {
            vout << format("%5d ")%(ind+1);
            vout << format("SigmaN2#%-2d ")%(i+1);
            double value = exp(Hyprms[ind]);
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            vout << format("%14.6e %14.6e")%(value)%HyprmsGrd[ind] << endl;
            ind++;
        }
    }

    if( NumOfOptPrms > 0 ){
        gnorm = sqrt(gnorm/(double)NumOfOptPrms);
    }

    vout << "# --- ---------- -------------- --------------" << endl;
    vout << format("# Gradient norm                 %14.6e")%gnorm << endl;
    vout << format("# Sum of targets                %14.6e")%logTarget << endl;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::ShowGPRStat(void)
{
    vout << endl;
    vout << format("%02d:Statistics per PMF accumulator ...")%State << endl;
    State++;

// setup parameters
    ScatterHyprms(Hyprms);

    vout << low;
    CreateGPREngine();
}

//------------------------------------------------------------------------------

void COptGPRHyprms::RunGPRAnalytical(void)
{
// setup parameters
    ScatterHyprms(Hyprms);
    HyprmsGrd.SetZero();

// get data
    logTarget = GetTarget();
    GetTargetDerivatives(HyprmsGrd);

// transform gradients
//    for(int ind=0; ind < NumOfOptPrms; ind++){
//        HyprmsGrd[ind] = 2.0*HyprmsGrd[ind]*Hyprms[ind];
//    }


    int ind = 0;
    for(int i=0; i < (int)SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ){
            HyprmsGrd[ind] = 2.0*HyprmsGrd[ind]*Hyprms[ind];
            ind++;
        }
    }
    for(int i=0; i < (int)CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] ){
            HyprmsGrd[ind] = 2.0*HyprmsGrd[ind]*Hyprms[ind];
            ind++;
        }
    }
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ){
            HyprmsGrd[ind] = 2.0*HyprmsGrd[ind]*Hyprms[ind];
            ind++;
        }
    }
    for(int i=0; i < (int)NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] ){
            HyprmsGrd[ind] = 2.0*HyprmsGrd[ind]*Hyprms[ind];
            ind++;
        }
    }
    for(int i=0; i < (int)SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ){
            // Hyprms[ind] = sqrt(SigmaN2[i]-Options.GetOptMinSigmaN2());
            HyprmsGrd[ind] = HyprmsGrd[ind]*exp(Hyprms[ind]);
            ind++;
        }
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::RunGPRNumerical(void)
{
// setup parameters
    HyprmsGrd.SetZero();
    logTarget = 0;

    CSimpleVector<double>   tmp_prms;
    tmp_prms.CreateVector(NumOfOptPrms);

// value
    tmp_prms = Hyprms;
    ScatterHyprms(tmp_prms);
    logTarget = GetTarget();

// derivatives
    double dh  = 1e-4;
    double dh2 = 1e-6;

    for(int i=0; i < NumOfOptPrms; i++){
        if( Options.GetOptCD5() ){
            double v1,v2,v3,v4;
            tmp_prms = Hyprms;
            tmp_prms[i] -= 2.0*(dh*Hyprms[i]+dh2);
            ScatterHyprms(tmp_prms);
            v1 = GetTarget();
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh*Hyprms[i]+dh2;
            ScatterHyprms(tmp_prms);
            v2 = GetTarget();
            tmp_prms = Hyprms;
            tmp_prms[i] += dh*Hyprms[i]+dh2;
            ScatterHyprms(tmp_prms);
            v3 = GetTarget();
            tmp_prms = Hyprms;
            tmp_prms[i] += 2.0*(dh*Hyprms[i]+dh2);
            ScatterHyprms(tmp_prms);
            v4 = GetTarget();
            HyprmsGrd[i] += (v1-8.0*v2+8.0*v3-v4)/(12.0*(dh*Hyprms[i]+dh2));
        } else {
            double lv,rv;
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh*Hyprms[i]+dh2;
            ScatterHyprms(tmp_prms);
            lv = GetTarget();
            tmp_prms = Hyprms;
            tmp_prms[i] += dh*Hyprms[i]+dh2;
            ScatterHyprms(tmp_prms);
            rv = GetTarget();
            HyprmsGrd[i] += (rv-lv)/(2.0*(dh*Hyprms[i]+dh2));
        }
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::ScatterHyprms(CSimpleVector<double>& hyprsm)
{
    // update parameters
    HyprmsEnabled.resize(NumOfPrms);

    int ind = 0;
    int i = 0;
// ---------------
    for(int k=0; k < (int)SigmaF2Enabled.size(); k++){
        if( SigmaF2Enabled[k] ){
            SigmaF2[k] = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinSigmaF2();
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
// ---------------
    for(int k=0; k < (int)CoVarEnabled.size(); k++){
        if( CoVarEnabled[k] ){
            CoVar[k] = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinCoVar();
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
            WFac[k] = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinWFac();
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
// ---------------
    for(int k=0; k < (int)NCorrEnabled.size(); k++){
        if( NCorrEnabled[k] ){
            NCorr[k] = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinNCorr();
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
// ---------------
    for(int k=0; k < (int)SigmaN2Enabled.size(); k++){
        if( SigmaN2Enabled[k] ){
            SigmaN2[k] = exp(hyprsm[ind]);
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitGPREngine(void)
{
    if( Options.GetArgRealm() == "dG/dx" ) {
        InitGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "dH/dx" ) {
        InitGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "-TdS/dx" ) {
        InitGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "mTdS/dx" ) {
        InitGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "dH" ) {
        InitGPREngine_dF();
    } else if( Options.GetArgRealm() == "GHS_dH_A" ) {
        InitGPREngine_GHS_dH_A();
    } else if( Options.GetArgRealm() == "cGHS_dH_A" ) {
        InitGPREngine_cGHS_dH_A();
    } else if( Options.GetArgRealm() == "GHS_dH_B" ) {
        InitGPREngine_GHS_dH_B();
    } else {
        CSmallString error;
        error << "unsupported realm: " <<  Options.GetArgRealm();
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitGPREngine_dF_dx(void)
{
    SigmaF2.CreateVector(1);
    NCorr.CreateVector(Accu->GetNumOfCVs());
    WFac.CreateVector(Accu->GetNumOfCVs());
    SigmaN2.CreateVector(Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitGPREngine_dF(void)
{
    SigmaF2.CreateVector(1);
    NCorr.CreateVector(Accu->GetNumOfCVs());
    WFac.CreateVector(Accu->GetNumOfCVs());
    SigmaN2.CreateVector(Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitGPREngine_GHS_dH_A(void)
{
    SigmaF2.CreateVector(3);
    WFac.CreateVector(Accu->GetNumOfCVs());
    SigmaN2.CreateVector(3*Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitGPREngine_cGHS_dH_A(void)
{
    SigmaF2.CreateVector(3);
    WFac.CreateVector(Accu->GetNumOfCVs());
    SigmaN2.CreateVector(3*Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitGPREngine_GHS_dH_B(void)
{
    SigmaF2.CreateVector(3);
    WFac.CreateVector(Accu->GetNumOfCVs());
    SigmaN2.CreateVector(3*Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CreateGPREngine(void)
{
    if( Options.GetArgRealm() == "dG/dx" ) {
        CreateGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "dH/dx" ) {
        CreateGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "-TdS/dx" ) {
        CreateGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "mTdS/dx" ) {
        CreateGPREngine_dF_dx();
    } else if( Options.GetArgRealm() == "dH" ) {
        CreateGPREngine_dF();
    } else if( Options.GetArgRealm() == "GHS_dH_A" ) {
        CreateGPREngine_GHS_dH_A();
    } else if( Options.GetArgRealm() == "cGHS_dH_A" ) {
        CreateGPREngine_cGHS_dH_A();
    } else if( Options.GetArgRealm() == "GHS_dH_B" ) {
        CreateGPREngine_GHS_dH_B();
    } else {
        CSmallString error;
        error << "unsupported realm: " <<  Options.GetArgRealm();
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CreateGPREngine_dF_dx(void)
{
    CEnergyDerProxyPtr proxy;

    if( Options.GetArgRealm() == "dG/dx" ) {
        if( CABFProxy_dG::IsCompatible(Accu) ){
           proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        } else if (CCSTProxy_dG::IsCompatible(Accu) ) {
            proxy    = CCSTProxy_dG_Ptr(new CCSTProxy_dG);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
    } else if( (Options.GetArgRealm() == "-TdS/dx") || (Options.GetArgRealm() == "mTdS/dx") ) {
        if( CABFProxy_mTdS::IsCompatible(Accu) ){
            proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        } else if (CCSTProxy_mTdS::IsCompatible(Accu) ) {
            proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
    } else {
            CSmallString error;
            error << "unsupported realm: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
    }

    proxy->Init(Accu);

    CIntegratorGPRPtr gpr = CIntegratorGPRPtr(new CIntegratorGPR);

    gpr->SetOutputES(FES);
    gpr->AddInputEnergyDerProxy(proxy);

    gpr->SetRCond(Options.GetOptRCond());

    gpr->SetIncludeError(false);
    gpr->SetNoEnergy(false);
    gpr->IncludeGluedAreas(false);

    gpr->SetLAMethod(Options.GetOptLAMethod());
    gpr->SetKernel(Options.GetOptGPRKernel());
    gpr->SetUseInv(Options.GetOptGPRUseInv());
    gpr->SetCalcLogPL(Options.GetOptGPRCalcLogPL() || (Target == EGOT_LOGPL));

// set parameters
    gpr->SetSigmaF2(SigmaF2);
    gpr->SetWFac(WFac);
    gpr->SetNCorr(NCorr);
    gpr->SetSigmaN2(SigmaN2);

// run integrator
    gpr->PrepForHyprmsGrd(true);
    gpr->Integrate(vout,false);

    GPREngine = gpr;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CreateGPREngine_dF(void)
{
    CEnergyProxyPtr proxy;

    if( Options.GetArgRealm() == "dH" ) {
        proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
    } else {
        CSmallString error;
        error << "unsupported realm: " <<  Options.GetArgRealm();
        RUNTIME_ERROR(error);
    }
    proxy->Init(Accu);

    CSmootherGPRPtr gpr = CSmootherGPRPtr(new CSmootherGPR);

    gpr->SetOutputES(FES);
    gpr->AddInputEnergyProxy(proxy);

    gpr->SetRCond(Options.GetOptRCond());

    gpr->SetIncludeError(false);

    gpr->SetLAMethod(Options.GetOptLAMethod());

    gpr->SetKernel(Options.GetOptGPRKernel());
    gpr->SetUseInv(Options.GetOptGPRUseInv());
    gpr->SetCalcLogPL(Options.GetOptGPRCalcLogPL() || (Target == EGOT_LOGPL));

// set parameters
    gpr->SetSigmaF2(SigmaF2);
    gpr->SetWFac(WFac);
    gpr->SetNCorr(NCorr);
    gpr->SetSigmaN2(SigmaN2);

// run interpolator
    gpr->PrepForHyprmsGrd(true);
    gpr->Interpolate(vout,false);

    GPREngine = gpr;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CreateGPREngine_GHS_dH_A(void)
{
    CEnergyDerProxyPtr proxy_dg;
    CEnergyProxyPtr    proxy_dh;
    CEnergyDerProxyPtr proxy_ds;

    if( Options.GetArgRealm() == "GHS_dH_A" ) {
        if( CABFProxy_dG::IsCompatible(Accu) ){
           proxy_dg = CABFProxy_dG_Ptr(new CABFProxy_dG);
           proxy_dg->Init(Accu);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm for dG/dx: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
        proxy_dh = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy_dh->Init(Accu);
        if( CABFProxy_mTdS::IsCompatible(Accu) ){
            proxy_ds    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy_ds->Init(Accu);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm for -TdS/dx: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
    } else {
        CSmallString error;
        error << "unsupported realm: " <<  Options.GetArgRealm();
        RUNTIME_ERROR(error);
    }

    CGHSIntegratorGPR0APtr gpr = CGHSIntegratorGPR0APtr(new CGHSIntegratorGPR0A);

    gpr->SetAccumulator(Accu);

    gpr->SetOutputFES(FES);
    gpr->SetOutputHES(HES);
    gpr->SetOutputSES(SES);

    gpr->SetGDerProxy(proxy_dg);
    gpr->SetHEneProxy(proxy_dh);
    gpr->SetSDerProxy(proxy_ds);

    gpr->SetRCond(Options.GetOptRCond());

    gpr->SetIncludeError(false);
    gpr->SetNoEnergy(false);

    gpr->SetLAMethod(Options.GetOptLAMethod());
    gpr->SetKernel(Options.GetOptGPRKernel());
    gpr->SetUseInv(Options.GetOptGPRUseInv());
    gpr->SetCalcLogPL(Options.GetOptGPRCalcLogPL() || (Target == EGOT_LOGPL));

// set parameters
    gpr->SetSigmaF2(SigmaF2);
    gpr->SetWFac(WFac);
    gpr->SetSigmaN2(SigmaN2);

// run integrator
    gpr->PrepForHyprmsGrd(true);
    gpr->Integrate(vout,false);

    GPREngine = gpr;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CreateGPREngine_cGHS_dH_A(void)
{
    CEnergyDerProxyPtr proxy_dg;
    CEnergyProxyPtr    proxy_dh;
    CEnergyDerProxyPtr proxy_ds;

    if( Options.GetArgRealm() == "cGHS_dH_A" ) {
        if( CABFProxy_dG::IsCompatible(Accu) ){
           proxy_dg = CABFProxy_dG_Ptr(new CABFProxy_dG);
           proxy_dg->Init(Accu);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm for dG/dx: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
        proxy_dh = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy_dh->Init(Accu);
        if( CABFProxy_mTdS::IsCompatible(Accu) ){
            proxy_ds    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy_ds->Init(Accu);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm for -TdS/dx: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
    } else {
        CSmallString error;
        error << "unsupported realm: " <<  Options.GetArgRealm();
        RUNTIME_ERROR(error);
    }

    CGHSIntegratorGPRcAPtr gpr = CGHSIntegratorGPRcAPtr(new CGHSIntegratorGPRcA);

    gpr->SetAccumulator(Accu);

    gpr->SetOutputFES(FES);
    gpr->SetOutputHES(HES);
    gpr->SetOutputSES(SES);

    gpr->SetGDerProxy(proxy_dg);
    gpr->SetHEneProxy(proxy_dh);
    gpr->SetSDerProxy(proxy_ds);

    gpr->SetRCond(Options.GetOptRCond());

    gpr->SetIncludeError(false);
    gpr->SetNoEnergy(false);

    gpr->SetLAMethod(Options.GetOptLAMethod());
    gpr->SetKernel(Options.GetOptGPRKernel());
    gpr->SetUseInv(Options.GetOptGPRUseInv());
    gpr->SetCalcLogPL(Options.GetOptGPRCalcLogPL() || (Target == EGOT_LOGPL));

    // FIXME
   // gpr->SetUseNumDiff(true);

// set parameters
    gpr->SetSigmaF2(SigmaF2);
    gpr->SetWFac(WFac);
    gpr->SetSigmaN2(SigmaN2);

// run integrator
    gpr->PrepForHyprmsGrd(true);
    gpr->Integrate(vout,false);

    GPREngine = gpr;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::CreateGPREngine_GHS_dH_B(void)
{
    CEnergyDerProxyPtr proxy_dg;
    CEnergyProxyPtr    proxy_dh;
    CEnergyDerProxyPtr proxy_ds;

    if( Options.GetArgRealm() == "GHS_dH_B" ) {
        if( CABFProxy_dG::IsCompatible(Accu) ){
           proxy_dg = CABFProxy_dG_Ptr(new CABFProxy_dG);
           proxy_dg->Init(Accu);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm for dG/dx: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
        proxy_dh = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
        proxy_dh->Init(Accu);
        if( CABFProxy_mTdS::IsCompatible(Accu) ){
            proxy_ds    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy_ds->Init(Accu);
        } else {
            CSmallString error;
            error << "incompatible method: " << Accu->GetMethod() << " with requested realm for -TdS/dx: " <<  Options.GetArgRealm();
            RUNTIME_ERROR(error);
        }
    } else {
        CSmallString error;
        error << "unsupported realm: " <<  Options.GetArgRealm();
        RUNTIME_ERROR(error);
    }

    CGHSIntegratorGPR0BPtr gpr = CGHSIntegratorGPR0BPtr(new CGHSIntegratorGPR0B);

    gpr->SetAccumulator(Accu);

    gpr->SetOutputFES(FES);
    gpr->SetOutputHES(HES);
    gpr->SetOutputSES(SES);

    gpr->SetGDerProxy(proxy_dg);
    gpr->SetHEneProxy(proxy_dh);
    gpr->SetSDerProxy(proxy_ds);

    gpr->SetRCond(Options.GetOptRCond());

    gpr->SetIncludeError(false);
    gpr->SetNoEnergy(false);

    gpr->SetLAMethod(Options.GetOptLAMethod());
    gpr->SetKernel(Options.GetOptGPRKernel());
    gpr->SetUseInv(Options.GetOptGPRUseInv());
    gpr->SetCalcLogPL(Options.GetOptGPRCalcLogPL() || (Target == EGOT_LOGPL));

// set parameters
    gpr->SetSigmaF2(SigmaF2);
    gpr->SetWFac(WFac);
    gpr->SetSigmaN2(SigmaN2);

// run integrator
    gpr->PrepForHyprmsGrd(true);
    gpr->Integrate(vout,false);

    GPREngine = gpr;
}

//------------------------------------------------------------------------------

double  COptGPRHyprms::GetTarget(void)
{
    vout << high;

    CreateGPREngine();

    double target = 0.0;
    switch(Target){
        case(EGOT_LOGML):
            // calculate logML
            target = GPREngine->GetLogML();
            vout << "      logML     = " << setprecision(5) << target << endl;
        break;
        case(EGOT_LOGPL):
            // calculate logLOO
            target = GPREngine->GetLogPL();
            vout << "      logPL     = " << setprecision(5) << target << endl;
        break;
    }

    vout << low;

    return(target);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::GetTargetDerivatives(CSimpleVector<double>& der)
{
    switch(Target){
        case(EGOT_LOGML):
            GPREngine->GetLogMLDerivatives(HyprmsEnabled,der);
        break;
        case(EGOT_LOGPL):
            GPREngine->GetLogPLDerivatives(HyprmsEnabled,der);
        break;
    }

}

//------------------------------------------------------------------------------

void COptGPRHyprms::WriteResults(int istep)
{
    vout << format("%6d %13e")%istep%logTarget;

    for(int i=0; i < (int)SigmaF2Enabled.size(); i++){
        if( SigmaF2Enabled[i] ) vout << format(" %10.4e")%SigmaF2[i];
    }
    for(int i=0; i < (int)CoVarEnabled.size(); i++){
        if( CoVarEnabled[i] ) vout << format(" %10.4e")%CoVar[i];
    }
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) vout << format(" %10.3f")%WFac[i];
    }
    for(int i=0; i < (int)NCorrEnabled.size(); i++){
        if( NCorrEnabled[i] ) vout << format(" %10.4e")%NCorr[i];
    }
    for(int i=0; i < (int)SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) vout << format(" %10.4e")%SigmaN2[i];
    }

    vout << endl;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = Accu->GetNumOfBins();
    int    sampled = 0;
    int    limit = 0;
    for(int ibin=0; ibin < Accu->GetNumOfBins(); ibin++) {
        if( Accu->GetNumOfSamples(ibin) > 0 ) {
            sampled++;
        }
        if( Accu->GetNumOfSamples(ibin) > Options.GetOptLimit() ) {
            limit++;
        } else {
            Accu->SetNumOfSamples(ibin,0);
        }
    }
    if( maxbins > 0 ){
        vout << " Sampled area: "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" ;
        vout << " ... Within limit: "
             << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%";
    }
    vout << endl;
}

//------------------------------------------------------------------------------

bool COptGPRHyprms::WriteHyperPrms(FILE* p_fout)
{
    if( p_fout == NULL ){
        return(false);
    }

    ScatterHyprms(Hyprms);

    if( fprintf(p_fout,"# GPR hyper-parameters\n") <= 0 ) return(false);

    for(size_t i=0; i < SigmaF2.GetLength(); i++ ){
        if( fprintf(p_fout,"SigmaF2#%-2ld = %16.10e\n",i+1,SigmaF2[i]) <= 0 ) return(false);
    }
    for(size_t i=0; i < CoVar.GetLength(); i++ ){
        if( fprintf(p_fout,"CoVar#%-2ld   = %16.10e\n",i+1,CoVar[i]) <= 0 ) return(false);
    }
    for(size_t i=0; i < WFac.GetLength(); i++ ){
        if( fprintf(p_fout,"WFac#%-2ld    = %16.10e\n",i+1,WFac[i]) <= 0 ) return(false);
    }
    for(size_t i=0; i < NCorr.GetLength(); i++ ){
        if( fprintf(p_fout,"NCorr#%-2ld   = %16.10e\n",i+1,NCorr[i]) <= 0 ) return(false);
    }
    for(size_t i=0; i < SigmaN2.GetLength(); i++ ){
        if( fprintf(p_fout,"SigmaN2#%-2ld = %16.10e\n",i+1,SigmaN2[i]) <= 0 ) return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::LoadGPRHyprms(void)
{
    ifstream fin;
    fin.open(Options.GetOptLoadHyprms());
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << Options.GetOptLoadHyprms();
        RUNTIME_ERROR(error);
    }

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
        if( key.find("SigmaF2#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode sigmaf2 key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( (cvind < 0) || (cvind >= (int)SigmaF2.GetLength()) ){
                CSmallString error;
                error << "sigmaf2 index " << (cvind+1) << " out-of-range 1-" << SigmaF2.GetLength();
                RUNTIME_ERROR(error);
            }
            SigmaF2[cvind] = value;
        } else if( key.find("CoVar#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode covar key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( (cvind < 0) || (cvind >= (int)CoVar.GetLength()) ){
                CSmallString error;
                error << "covar index " << (cvind+1) << " out-of-range 1-" << CoVar.GetLength();
                RUNTIME_ERROR(error);
            }
            CoVar[cvind] = value;
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
            if( (cvind < 0) || (cvind >= (int)WFac.GetLength()) ){
                CSmallString error;
                error << "wfac index " << (cvind+1) << " out-of-range 1-" << WFac.GetLength();
                RUNTIME_ERROR(error);
            }
            WFac[cvind] = value;
        } else if( key.find("NCorr#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode ncorr key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( (cvind < 0) || (cvind >= (int)NCorr.GetLength()) ){
                CSmallString error;
                error << "ncorr index " << (cvind+1) << " out-of-range 1-" << NCorr.GetLength();
                RUNTIME_ERROR(error);
            }
            NCorr[cvind] = value;
        } else if( key.find("SigmaN2#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string ssigman2;
            int    cvind;
            kstr >> ssigman2 >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode sigman2 key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            if( (cvind < 0) || (cvind >= (int)SigmaN2.GetLength()) ){
                CSmallString error;
                error << "sigman2 index " << (cvind+1) << " out-of-range 1-" << SigmaN2.GetLength();
                RUNTIME_ERROR(error);
            }
            SigmaN2[cvind] = value;
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

void COptGPRHyprms::Finalize(void)
{
    // close files if they are own by program
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    CSmallTime dur;
    dur = dt - StartTime;

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# gpr-opthyprms terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

