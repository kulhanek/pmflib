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

    CSmallString output = Options.GetProgArg(Options.GetNumberOfProgArgs()-1);
    if( output != "-") {
        vout << "# Optimized GPR hyprms (out): " << output << endl;
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
    CProxyRealmPtr curr;
    vout << format("%02d:Loading PMF accumulators ...")%State << endl;
    State++;
    int naccu = 1;
    for(int i=0; i < Options.GetNumberOfProgArgs()-1; i++){
        CSmallString name = Options.GetProgArg(i);
        if( IsRealm(name) ){
            if( curr == NULL ){
                CSmallString error;
                error << "no accumulators defined for realm'" << name << "'";
                ES_ERROR(error);
                return(false);
            }
            curr->Name = name;
            curr = CProxyRealmPtr();
        } else {
            vout << format("** PMFAccumulator #%05d: %s")%(naccu)%string(name) << endl;
            CPMFAccumulatorPtr p_accu(new CPMFAccumulator);
            try {
                p_accu->Load(name);
            } catch(...) {
                CSmallString error;
                error << "unable to load the input PMF accumulator file '" << name << "'";
                ES_ERROR(error);
                return(false);
            }
            if( curr == NULL ){
                curr = CProxyRealmPtr(new CProxyRealm);
                RealmProxies.push_back(curr);
            }
            curr->Accumulators.push_back(p_accu);
            naccu++;
        }
    }
    vout << "   Done" << endl;

// realms
    vout << endl;
    vout << format("%02d:Initializing realms ...")%State << endl;
    State++;
    for(size_t i=0; i < RealmProxies.size(); i++){
        CProxyRealmPtr realm = RealmProxies[i];
        vout << format("   ** ----------> Realm #%02d: %s")%(i+1)%realm->Name << endl;
        vout << format("      # of PMF accumulators: %d")%realm->Accumulators.size() << endl;
        InitRealm(realm);
    }

// statistics
    vout << endl;
    vout << format("%02d:PMF accumulator statistics ...")%State << endl;
    State++;
    PrintSampledStat();
    vout << "   Done" << endl;

// run optimization
    bool result = true;
    InitOptimizer();
    if( ! Options.GetOptTest() ){
        result = Optimize();
    } else {
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

bool COptGPRHyprms::IsRealm(const CSmallString& name)
{
    if( name == "dG" )      return(true);
    if( name == "-TdS" )    return(true);
    if( name == "mTdS" )    return(true);
    if( name == "dH" )      return(true);
    return(false);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::InitRealm(CProxyRealmPtr realm)
{
// init energyder proxy
    if( realm->Name == "dG" ){
        for(size_t i=0; i < realm->Accumulators.size(); i++){
            CPMFAccumulatorPtr accu = realm->Accumulators[i];
            CEnergyDerProxyPtr proxy;
            if( CABFProxy_dG::IsCompatible(accu) ){
               proxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
            } else if (CCSTProxy_dG::IsCompatible(accu) ) {
                proxy    = CCSTProxy_dG_Ptr(new CCSTProxy_dG);
            } else {
                CSmallString error;
                error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm->Name;
                RUNTIME_ERROR(error);
            }
            proxy->Init(accu);
            realm->DerProxies.push_back(proxy);
        }
// -----------------------------------------------
    } else if ( (realm->Name == "-TdS") || (realm->Name == "mTdS") ) {
        for(size_t i=0; i < realm->Accumulators.size(); i++){
            CPMFAccumulatorPtr accu = realm->Accumulators[i];
            CEnergyDerProxyPtr proxy;
            if( CABFProxy_mTdS::IsCompatible(accu) ){
                proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            } else if (CCSTProxy_mTdS::IsCompatible(accu) ) {
                proxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
            } else {
                CSmallString error;
                error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm->Name;
                RUNTIME_ERROR(error);
            }
            proxy->Init(accu);
            realm->DerProxies.push_back(proxy);
        }
    } else if ( realm->Name == "dH" ) {
        for(size_t i=0; i < realm->Accumulators.size(); i++){
            CPMFAccumulatorPtr accu = realm->Accumulators[i];
            CEnergyProxyPtr proxy    = CPMFProxy_dH_Ptr(new CPMFProxy_dH);
            proxy->Init(accu);
            realm->EnergyProxies.push_back(proxy);
        }
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << realm->Name;
        RUNTIME_ERROR(error);
    }
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
    SigmaF2Enabled  = Options.GetOptSigmaF2Enabled();
    NCorrEnabled    = Options.GetOptNCorrEnabled();
    DecodeEList(Options.GetOptWFacEnabled(),WFacEnabled,"--enablewfac");
    DecodeEList(Options.GetOptSigmaN2Enabled(),SigmaN2Enabled,"--enablesigman2");

// number of optimized parameters
    NumOfOptPrms = 0;
    NumOfPrms = 0;
// sigmaf2
    if( SigmaF2Enabled ) NumOfOptPrms++;
    NumOfPrms++;
// ncorr
    if( NCorrEnabled ) NumOfOptPrms++;
    NumOfPrms++;
// wfac
    for(size_t i=0; i < WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) NumOfOptPrms++;
    }
    NumOfPrms += NCVs;
// sigman2
    for(size_t i=0; i < SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) NumOfOptPrms++;
    }
    NumOfPrms += NCVs;

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
        SigmaF2 = Options.GetOptSigmaF2();
        NCorr   = Options.GetOptNCorr();
        DecodeVList(Options.GetOptWFac(),WFac,"--wfac",3.0);
        for(size_t i=0; i < WFac.GetLength(); i++){
            if( WFac[i] < Options.GetOptMinWFac() ){
                RUNTIME_ERROR("--wfac has to be greater than --minwfac");
            }
        }
        DecodeVList(Options.GetOptSigmaN2(),SigmaN2,"--sigman2",0.1);
        for(size_t i=0; i < SigmaN2.GetLength(); i++){
            if( SigmaN2[i] < Options.GetOptMinSigmaN2() ){
                RUNTIME_ERROR("--sigman2 has to be greater than --minsigman2");
            }
        }
    }

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
        vout << format("      SigmaF2    = %10.4f")%SigmaF2 << endl;
        vout << format("      NCorr      = %10.4f")%NCorr   << endl;

    for(int k=0; k < NCVs; k++ ){
        vout << format("      WFac#%-2d    = %10.4f")%(k+1)%WFac[k] << endl;
    }
    for(int k=0; k < NCVs; k++ ){
        vout << format("      SigmaN2#%-2d = %10.4e")%(k+1)%SigmaN2[k] << endl;
    }

    if( ! Options.GetOptTest() ){
        vout << endl;
        switch(Target){
            case(EGOT_LOGML):
                vout << "# step         logML";
            break;
            case(EGOT_LOGPL):
                vout << "# step         logPL";
            break;
        }

        if( SigmaF2Enabled )    vout << "    SigmaF2";
        if( NCorrEnabled )      vout << "      NCorr";

        for(size_t i=0; i < WFacEnabled.size(); i++){
            if( WFacEnabled[i] ) vout << format("     WFac#%-2d")%(i+1);
        }
        for(size_t i=0; i < SigmaN2Enabled.size(); i++){
            if( SigmaN2Enabled[i] ) vout << format(" SigmaN2#%-2d")%(i+1);
        }
        vout << endl;

        vout << "# ---- -------------";
        if( SigmaF2Enabled )    vout << " ----------";
        if( NCorrEnabled )      vout << " ----------";
        for(size_t i=0; i < WFacEnabled.size(); i++){
            if( WFacEnabled[i] ) vout << " ----------";
        }
        for(size_t i=0; i < SigmaN2Enabled.size(); i++){
            if( SigmaN2Enabled[i] ) vout << " ----------";
        }
        vout << endl;
    }

    Hyprms.CreateVector(NumOfOptPrms);
    HyprmsGrd.CreateVector(NumOfOptPrms);

// set initial hyprms
    int ind = 0;
    if( SigmaF2Enabled ){
        Hyprms[ind] = sqrt(SigmaF2-Options.GetOptMinSigmaF2());
        ind++;
    }
    if( NCorrEnabled ){
        Hyprms[ind] = sqrt(NCorr-Options.GetOptMinNCorr());
        ind++;
    }
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ){
            Hyprms[ind] = sqrt(WFac[i]-Options.GetOptMinWFac());
            ind++;
        }
    }
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ){
            Hyprms[ind] = sqrt(SigmaN2[i]-Options.GetOptMinSigmaN2());
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

void COptGPRHyprms::DecodeEList(const CSmallString& spec, std::vector<bool>& elist,const CSmallString& optionname)
{
    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( (int)slist.size() > NCVs ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << NCVs << ")";
        RUNTIME_ERROR(error);
    }

    elist.resize(NCVs);

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
    for(int i=slist.size(); i < NCVs; i++){
        elist[i] = last_st;
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::DecodeVList(const CSmallString& spec, CSimpleVector<double>& vlist,const CSmallString& optionname,double defv)
{
    string          sspecen(spec);
    vector<string>  slist;

    split(slist,sspecen,is_any_of("x"),token_compress_on);

    if( (int)slist.size() > NCVs ){
        CSmallString error;
        error << "too many flags (" << slist.size() << ") for " << optionname << " than required (" << NCVs << ")";
        RUNTIME_ERROR(error);
    }

    vlist.CreateVector(NCVs);

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
    for(int i=slist.size(); i < NCVs; i++){
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
            Hyprms[i] = 1.0;
        }
    }

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

        if( istep > 1 ){
            if( fabs(logTarget - last_logtrg) < Options.GetOptTermVal() ){
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
        }
        last_logtrg = logTarget;
    }

    if( noptsteps > 0 ) {
        vout << "# ---- -------------";
        if( SigmaF2Enabled )    vout << " ----------";
        if( NCorrEnabled )      vout << " ----------";
        for(size_t i=0; i < WFacEnabled.size(); i++){
            if( WFacEnabled[i] ) vout << " ----------";
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
            if( prm == 0 ){
                vout << "SigmaF2    ";
            } else if( prm == 1 ){
                vout << "NCorr      ";
            } else if( (prm >= 2) && (prm < NCVs+2) ) {
                int cv = prm - 2;
                vout << format("WFac#%-2d    ")%(cv+1);
            } else {
                int cv = prm - (NCVs+2);
                vout << format("SigmaN2#%-2d ")%(cv+1);
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
    double dh = 1e-3;

    for(int i=0; i < NumOfOptPrms; i++){

        vout << "   + perturbation for hyprm: " << setw(3) << (i+1) << endl;
        grd1[i].SetZero();
        tmp_prms = Hyprms;
        tmp_prms[i] += dh;
        ScatterHyprms(tmp_prms);

        // calculate gradient and logML
        for(size_t j=0; j < RealmProxies.size(); j++){
            CProxyRealmPtr proxy = RealmProxies[j];

        // IntegratorGPR
            if( proxy->DerProxies.size() > 0 ) {
                CIntegratorGPR gpr;
                gpr.PrepForHyprmsGrd(true);
                GetTargetFromIntegrator(gpr,proxy);
                switch(Target){
                    case(EGOT_LOGML):
                        gpr.GetLogMLDerivatives(HyprmsEnabled,grd1[i]);
                    break;
                    case(EGOT_LOGPL):
                        gpr.GetLogPLDerivatives(HyprmsEnabled,grd1[i]);
                    break;
                }
        // SmootherGPR
            } else if( proxy->EnergyProxies.size() > 0 ) {
                CSmootherGPR gpr;
                gpr.PrepForHyprmsGrd(true);
                GetTargetFromSmoother(gpr,proxy);
                switch(Target){
                    case(EGOT_LOGML):
                        gpr.GetLogMLDerivatives(HyprmsEnabled,grd1[i]);
                    break;
                    case(EGOT_LOGPL):
                        gpr.GetLogPLDerivatives(HyprmsEnabled,grd1[i]);
                    break;
                }
            } else {
                RUNTIME_ERROR("undefined proxy")
            }
        }

        vout << "   - perturbation for hyprm: " << setw(3) << (i+1) << endl;
        grd2[i].SetZero();
        tmp_prms = Hyprms;
        tmp_prms[i] -= dh;
        ScatterHyprms(tmp_prms);


// calculate gradient and logML
        for(size_t j=0; j < RealmProxies.size(); j++){
            CProxyRealmPtr proxy = RealmProxies[j];

        // IntegratorGPR
            if( proxy->DerProxies.size() > 0  ) {
                CIntegratorGPR gpr;
                gpr.PrepForHyprmsGrd(true);
                GetTargetFromIntegrator(gpr,proxy);
                switch(Target){
                    case(EGOT_LOGML):
                        gpr.GetLogMLDerivatives(HyprmsEnabled,grd2[i]);
                    break;
                    case(EGOT_LOGPL):
                        gpr.GetLogPLDerivatives(HyprmsEnabled,grd2[i]);
                    break;
                }
        // SmootherGPR
            } else if( proxy->EnergyProxies.size() > 0 ) {
                CSmootherGPR gpr;
                gpr.PrepForHyprmsGrd(true);
                GetTargetFromSmoother(gpr,proxy);
                switch(Target){
                    case(EGOT_LOGML):
                        gpr.GetLogMLDerivatives(HyprmsEnabled,grd2[i]);
                    break;
                    case(EGOT_LOGPL):
                        gpr.GetLogPLDerivatives(HyprmsEnabled,grd2[i]);
                    break;
                }
            } else {
                RUNTIME_ERROR("undefined proxy")
            }
        }
    }

// hessian by central differences from gradients
// https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
    for(int i=0; i < NumOfOptPrms; i++){
        for(int j=0; j < NumOfOptPrms; j++){
            Hessian[i][j] = (grd1[j][i] - grd2[j][i])/(4.0*dh) + (grd1[i][j] - grd2[i][j])/(4.0*dh);
        }
    }

// find eigenvalues
    CSciLapack::syev('V','U',Hessian,EigenValues);

    vout << endl;
    vout << "# idx    hess eigval";

    int ind = 0;
    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( ! HyprmsEnabled[prm] ) continue;
        if( prm == 0 ){
            vout << " SigmaF2   ";
        } else if( prm == 1 ){
            vout << " NCorr     ";
        } else if( (prm >= 2) && (prm < NCVs+2) ) {
            int cv = prm - 2;
            vout << format(" WFac#%-2d   ")%(cv+1);
        } else {
            int cv = prm - (NCVs+2);
            vout << format(" SigmaN2#%-2d")%(cv+1);
        }
        ind++;
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

    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( ! HyprmsEnabled[prm] ) continue;
        vout << format("%5d ")%(ind+1);
        double value = 0.0;
        if( prm == 0 ){
            vout << "SigmaF2    ";
            value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinSigmaF2();
        } else if( prm == 1 ){
            vout << "NCorr      ";
            value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinNCorr();
        } else if( (prm >= 2) && (prm < NCVs+2)) {
            int cv = prm - 2;
            vout << format("WFac#%-2d    ")%(cv+1);
            value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinWFac();
        } else {
            int cv = prm - (NCVs+2);
            vout << format("SigmaN2#%-2d ")%(cv+1);
            value = Hyprms[ind]*Hyprms[ind] + Options.GetOptMinSigmaN2();
        }
        gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
        vout << format("%14.6e %14.6e")%(value)%HyprmsGrd[ind] << endl;
        ind++;
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

// run GPR integration
    for(size_t i=0; i < RealmProxies.size(); i++){
        CProxyRealmPtr proxy = RealmProxies[i];

        vout << format("** PMF Realm Set #%02d ...")%(i+1) << endl;

    // IntegratorGPR
        if( proxy->DerProxies.size() > 0 ) {
            CIntegratorGPR   gpr;

            CEnergySurfacePtr fes(new CEnergySurface);
            fes->Allocate(proxy->Accumulators.front());

            gpr.SetOutputES(fes);
            for(size_t i=0; i < proxy->DerProxies.size(); i++){
                gpr.AddInputEnergyDerProxy(proxy->DerProxies[i]);
            }

            gpr.SetRCond(Options.GetOptRCond());

            gpr.SetIncludeError(false);
            gpr.SetNoEnergy(false);
            gpr.IncludeGluedAreas(false);

            gpr.SetLAMethod(Options.GetOptLAMethod());
            gpr.SetKernel(Options.GetOptGPRKernel());
            gpr.SetUseInv(Options.GetOptGPRUseInv());
            gpr.SetCalcLogPL(Options.GetOptGPRCalcLogPL() || Target == EGOT_LOGPL);

            if( Options.IsOptGlobalMinSet() ){
                gpr.SetGlobalMin(Options.GetOptGlobalMin());
            }

        // run integrator
            gpr.SetSigmaF2(SigmaF2);
            gpr.SetNCorr(NCorr);
            gpr.SetWFac(WFac);
            gpr.SetSigmaN2(SigmaN2);
            gpr.Integrate(vout,false);
    // SmootherGPR
        } else if( proxy->EnergyProxies.size() > 0 ) {
            CSmootherGPR   gpr;

            CEnergySurfacePtr fes(new CEnergySurface);
            fes->Allocate(proxy->Accumulators.front());

            gpr.SetOutputES(fes);
            for(size_t i=0; i < proxy->EnergyProxies.size(); i++){
                gpr.AddInputEnergyProxy(proxy->EnergyProxies[i]);
            }

            gpr.SetRCond(Options.GetOptRCond());

            gpr.SetIncludeError(false);

            gpr.SetLAMethod(Options.GetOptLAMethod());

            gpr.SetKernel(Options.GetOptGPRKernel());
            gpr.SetUseInv(Options.GetOptGPRUseInv());
            gpr.SetCalcLogPL(Options.GetOptGPRCalcLogPL() || Target == EGOT_LOGPL);

            if( Options.IsOptGlobalMinSet() ){
                gpr.SetGlobalMin(Options.GetOptGlobalMin());
            }

        // run integrator
            gpr.SetSigmaF2(SigmaF2);
            gpr.SetNCorr(NCorr);
            gpr.SetWFac(WFac);
            gpr.SetSigmaN2(SigmaN2);
            gpr.Interpolate(vout,false);
        } else {
            RUNTIME_ERROR("undefined proxy")
        }

   }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::RunGPRAnalytical(void)
{
// setup parameters
    ScatterHyprms(Hyprms);

    logTarget = 0;
    HyprmsGrd.SetZero();

// calculate gradient and logML
    for(size_t i=0; i < RealmProxies.size(); i++){
        CProxyRealmPtr proxy = RealmProxies[i];

    // IntegratorGPR
        if( proxy->DerProxies.size() > 0 ) {
            CIntegratorGPR gpr;
            gpr.PrepForHyprmsGrd(true);
            logTarget += GetTargetFromIntegrator(gpr,proxy);
            switch(Target){
                case(EGOT_LOGML):
                    gpr.GetLogMLDerivatives(HyprmsEnabled,HyprmsGrd);
                break;
                case(EGOT_LOGPL):
                    gpr.GetLogPLDerivatives(HyprmsEnabled,HyprmsGrd);
                break;
            }
    // SmootherGPR
        } else if( proxy->EnergyProxies.size() > 0 ) {
            CSmootherGPR gpr;
            gpr.PrepForHyprmsGrd(true);
            logTarget += GetTargetFromSmoother(gpr,proxy);
            switch(Target){
                case(EGOT_LOGML):
                    gpr.GetLogMLDerivatives(HyprmsEnabled,HyprmsGrd);
                break;
                case(EGOT_LOGPL):
                    gpr.GetLogPLDerivatives(HyprmsEnabled,HyprmsGrd);
                break;
            }
        } else {
            RUNTIME_ERROR("undefined proxy")
        }
    }

    // transform gradients
    for(int ind=0; ind < NumOfOptPrms; ind++){
        HyprmsGrd[ind] = 2.0*HyprmsGrd[ind]*Hyprms[ind];
    }
}

//------------------------------------------------------------------------------

void COptGPRHyprms::RunGPRNumerical(void)
{
    logTarget = 0;
    HyprmsGrd.SetZero();
    for(size_t i=0; i < RealmProxies.size(); i++){
        CProxyRealmPtr proxy = RealmProxies[i];

    // IntegratorGPR
        if( proxy->DerProxies.size() > 0 ) {
            logTarget += RunGPRNumericalIntegrator(proxy,HyprmsGrd);
    // SmootherGPR
        } else if( proxy->EnergyProxies.size() > 0 ) {
            logTarget += RunGPRNumericalSmoother(proxy,HyprmsGrd);
        } else {
            RUNTIME_ERROR("undefined proxy")
        }
    }
}

//------------------------------------------------------------------------------

double COptGPRHyprms::RunGPRNumericalIntegrator(CProxyRealmPtr derproxy,CSimpleVector<double>& der)
{
    CSimpleVector<double>   tmp_prms;
    tmp_prms.CreateVector(NumOfOptPrms);

    double lml = 0.0;

// value
    tmp_prms = Hyprms;
    ScatterHyprms(tmp_prms);
    {
        CIntegratorGPR gpr;
        lml = GetTargetFromIntegrator(gpr,derproxy);
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
                CIntegratorGPR gpr;
                v1 = GetTargetFromIntegrator(gpr,derproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CIntegratorGPR gpr;
                v2 = GetTargetFromIntegrator(gpr,derproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CIntegratorGPR gpr;
                v3 = GetTargetFromIntegrator(gpr,derproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += 2.0*dh;
            ScatterHyprms(tmp_prms);
            {
                CIntegratorGPR gpr;
                v4 = GetTargetFromIntegrator(gpr,derproxy);
            }
            der[i] += (v1-8.0*v2+8.0*v3-v4)/(12.0*dh);
        } else {
            double lv,rv;
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CIntegratorGPR gpr;
                lv = GetTargetFromIntegrator(gpr,derproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CIntegratorGPR gpr;
                rv = GetTargetFromIntegrator(gpr,derproxy);
            }
            der[i] += (rv-lv)/(2.0*dh);
        }

    }

    return(lml);
}

//------------------------------------------------------------------------------

double COptGPRHyprms::RunGPRNumericalSmoother(CProxyRealmPtr eneproxy,CSimpleVector<double>& der)
{
    CSimpleVector<double>   tmp_prms;
    tmp_prms.CreateVector(NumOfOptPrms);

    double lml = 0.0;

// value
    tmp_prms = Hyprms;
    ScatterHyprms(tmp_prms);
    {
        CSmootherGPR gpr;
        lml = GetTargetFromSmoother(gpr,eneproxy);
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
                CSmootherGPR gpr;
                v1 = GetTargetFromSmoother(gpr,eneproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CSmootherGPR gpr;
                v2 = GetTargetFromSmoother(gpr,eneproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CSmootherGPR gpr;
                v3 = GetTargetFromSmoother(gpr,eneproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += 2.0*dh;
            ScatterHyprms(tmp_prms);
            {
                CSmootherGPR gpr;
                v4 = GetTargetFromSmoother(gpr,eneproxy);
            }
            der[i] += (v1-8.0*v2+8.0*v3-v4)/(12.0*dh);
        } else {
            double lv,rv;
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CSmootherGPR gpr;
                lv = GetTargetFromSmoother(gpr,eneproxy);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CSmootherGPR gpr;
                rv = GetTargetFromSmoother(gpr,eneproxy);
            }
            der[i] += (rv-lv)/(2.0*dh);
        }

    }

    return(lml);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::ScatterHyprms(CSimpleVector<double>& hyprsm)
{
    // update parameters
    HyprmsEnabled.resize(NumOfPrms);

    int ind = 0;
    int i = 0;
// ---------------
    if( SigmaF2Enabled ){
        SigmaF2 = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinSigmaF2();
        ind++;
        HyprmsEnabled[i] = true;
    } else {
        HyprmsEnabled[i] = false;
    }
    i++;
// ---------------
    if( NCorrEnabled ){
        NCorr = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinNCorr();
        ind++;
        HyprmsEnabled[i] = true;
    } else {
        HyprmsEnabled[i] = false;
    }
    i++;
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
    for(int k=0; k < (int)SigmaN2Enabled.size(); k++){
        if( SigmaN2Enabled[k] ){
            SigmaN2[k] = hyprsm[ind]*hyprsm[ind] + Options.GetOptMinSigmaN2();
            ind++;
            HyprmsEnabled[i] = true;
        } else {
            HyprmsEnabled[i] = false;
        }
        i++;
    }
}

//------------------------------------------------------------------------------

double COptGPRHyprms::GetTargetFromIntegrator(CIntegratorGPR& gpr,CProxyRealmPtr proxy)
{
    CEnergySurfacePtr fes(new CEnergySurface);
    fes->Allocate(proxy->Accumulators.front());

    gpr.SetOutputES(fes);
    for(size_t i=0; i < proxy->DerProxies.size(); i++){
        gpr.AddInputEnergyDerProxy(proxy->DerProxies[i]);
    }

    gpr.SetRCond(Options.GetOptRCond());

    gpr.SetIncludeError(false);
    gpr.SetNoEnergy(true);
    gpr.IncludeGluedAreas(false);

    gpr.SetLAMethod(Options.GetOptLAMethod());
    gpr.SetUseInv(Options.GetOptGPRUseInv());
    gpr.SetKernel(Options.GetOptGPRKernel());

    if( Options.IsOptGlobalMinSet() ){
        gpr.SetGlobalMin(Options.GetOptGlobalMin());
    }

    if( Target == EGOT_LOGPL) gpr.SetCalcLogPL(true);

// setup integrator
    gpr.SetSigmaF2(SigmaF2);
    gpr.SetNCorr(NCorr);
    gpr.SetWFac(WFac);
    gpr.SetSigmaN2(SigmaN2);
    vout << high;
    gpr.Integrate(vout,true);

    double target = 0.0;
    if( Options.GetOptTarget() == "logml" ){
        // calculate logML
        target = gpr.GetLogML();
        vout << "      logML     = " << setprecision(5) << target << endl;
    } else if ( Options.GetOptTarget() == "logpl" ) {
        // calculate logLOO
        target = gpr.GetLogPL();
        vout << "      logPL     = " << setprecision(5) << target << endl;
    }

    vout << low;

    return(target);
}

//------------------------------------------------------------------------------

double  COptGPRHyprms::GetTargetFromSmoother(CSmootherGPR& gpr,CProxyRealmPtr proxy)
{
    CEnergySurfacePtr fes(new CEnergySurface);
    fes->Allocate(proxy->Accumulators.front());

    gpr.SetOutputES(fes);
    for(size_t i=0; i < proxy->EnergyProxies.size(); i++){
        gpr.AddInputEnergyProxy(proxy->EnergyProxies[i]);
    }

    gpr.SetRCond(Options.GetOptRCond());

    gpr.SetIncludeError(false);

    gpr.SetLAMethod(Options.GetOptLAMethod());
    gpr.SetUseInv(Options.GetOptGPRUseInv());

    gpr.SetKernel(Options.GetOptGPRKernel());

    if( Options.IsOptGlobalMinSet() ){
        gpr.SetGlobalMin(Options.GetOptGlobalMin());
    }

    if( Target == EGOT_LOGPL) gpr.SetCalcLogPL(true);

// setup integrator
    gpr.SetSigmaF2(SigmaF2);
    gpr.SetNCorr(NCorr);
    gpr.SetWFac(WFac);
    gpr.SetSigmaN2(SigmaN2);
    vout << high;

    gpr.Interpolate(vout,true);

    double target = 0.0;
    if( Options.GetOptTarget() == "logml" ){
        // calculate logML
        target = gpr.GetLogML();
        vout << "      logML     = " << setprecision(5) << target << endl;
    } else if ( Options.GetOptTarget() == "logpl" ) {
        // calculate logLOO
        target = gpr.GetLogPL();
        vout << "      logPL     = " << setprecision(5) << target << endl;
    }

    vout << low;

    return(target);
}

//------------------------------------------------------------------------------

void COptGPRHyprms::WriteResults(int istep)
{
    vout << format("%6d %13e")%istep%logTarget;
    if( SigmaF2Enabled ) vout << format(" %10.3f")%SigmaF2;
    if( NCorrEnabled )   vout << format(" %10.3f")%NCorr;
    for(int i=0; i < (int)WFacEnabled.size(); i++){
        if( WFacEnabled[i] ) vout << format(" %10.3f")%WFac[i];
    }
    for(int i=0; i < (int)SigmaN2Enabled.size(); i++){
        if( SigmaN2Enabled[i] ) vout << format(" %10.4e")%SigmaN2[i];
    }
    vout << endl;
}

//------------------------------------------------------------------------------

void COptGPRHyprms::PrintSampledStat(void)
{
    for(size_t r=0; r < RealmProxies.size(); r++){
        vout << format("** -------> Realm #%02d: %s")%(r+1)%RealmProxies[r]->Name << endl;
        for(size_t i=0; i < RealmProxies[r]->Accumulators.size(); i++){
            CPMFAccumulatorPtr accu = RealmProxies[r]->Accumulators[i];
            vout << format("   ** PMF Accumulator #%05d ...")%(i+1);
            // calculate sampled area
            double maxbins = accu->GetNumOfBins();
            int    sampled = 0;
            int    limit = 0;
            for(int ibin=0; ibin < accu->GetNumOfBins(); ibin++) {
                if( accu->GetNumOfSamples(ibin) > 0 ) {
                    sampled++;
                }
                if( accu->GetNumOfSamples(ibin) > Options.GetOptLimit() ) {
                    limit++;
                } else {
                    accu->SetNumOfSamples(ibin,0);
                }
            }
            if( maxbins > 0 ){
                vout << " Sampled area: "
                     << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" ;
                vout << " ... Within limit: "
                     << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%";
            }
            vout << endl;
            int ncvs = accu->GetNumOfCVs();
            if( i == 0 ){
                NCVs = ncvs;
            }
            if( NCVs != ncvs ){
                CSmallString error;
                error << "inconsistent dimensions (NCVs) of PMF accumulator: " << ncvs << "; the first accu: " << NCVs;
                RUNTIME_ERROR(error);
            }
       }
    }
}

//------------------------------------------------------------------------------

bool COptGPRHyprms::WriteHyperPrms(FILE* p_fout)
{
    if( p_fout == NULL ){
        return(false);
    }

    ScatterHyprms(Hyprms);

    if( fprintf(p_fout,"# GPR hyper-parameters\n") <= 0 ) return(false);
        if( fprintf(p_fout,"SigmaF2    = %16.10e\n",SigmaF2) <= 0 ) return(false);
        if( fprintf(p_fout,"NCorr      = %16.10e\n",NCorr) <= 0 ) return(false);

    for(size_t i=0; i < WFac.GetLength(); i++ ){
        if( fprintf(p_fout,"WFac#%-2ld    = %16.10e\n",i+1,WFac[i]) <= 0 ) return(false);
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

    WFac.CreateVector(NCVs);
    SigmaN2.CreateVector(NCVs);

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
            NCorr = value;
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
            if( (cvind < 0) || (cvind >= NCVs) ){
                CSmallString error;
                error << "wfac index " << (cvind+1) << " out-of-range 1-" << NCVs;
                RUNTIME_ERROR(error);
            }
            WFac[cvind] = value;
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
            if( (cvind < 0) || (cvind >= NCVs) ){
                CSmallString error;
                error << "sigman2 index " << (cvind+1) << " out-of-range 1-" << NCVs;
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

