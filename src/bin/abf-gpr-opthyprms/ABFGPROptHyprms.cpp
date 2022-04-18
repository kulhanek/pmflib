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
#include "ABFGPROptHyprms.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iomanip>
#include <Vector.hpp>
#include <SciLapack.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFGPROptHyprms)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFGPROptHyprms::CABFGPROptHyprms(void)
{
    State       = 1;
    NSTLimit    = 0;
    NumOfBins   = 0;
    NumOfCVs    = 0;
    KRank       = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFGPROptHyprms::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CIntOpts
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

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-gpr-opthyprms (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgInAccuName() != "-") {
        vout << "# Input PMF accumulator (in)      : " << Options.GetArgInAccuName() << endl;
    } else {
        vout << "# Input PMF accumulator (in)      : - (standard input)" << endl;
    }
    if(Options.GetArgHyprmsName() != "-") {
        vout << "# Resulting hyperparameters (out) : " << Options.GetArgHyprmsName() << endl;
    } else {
        vout << "# Resulting hyperparameters (out) : - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFGPROptHyprms::Run(void)
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
    TimeStep    = InAccu->GetTimeStep();

    vout << "   Done." << endl;

    InAccu->PrintInfo(vout);

// -----------------------------------------------------------------------------


    if( Options.GetOptOptTarget() == "logml" ){
        OptTarget = EGOT_LOGML;
    } else if( Options.GetOptOptTarget() == "logpl" ){
        OptTarget = EGOT_LOGPL;
    } else {
        RUNTIME_ERROR("unsupported opt target");
    }

    if( Options.GetOptOptRealm() == "CVS" ){
        OptRealm = EGOR_CVS;
    } else if( Options.GetOptOptRealm() == "ICF" ){
        OptRealm = EGOR_ICF;
    } else if( Options.GetOptOptRealm() == "KIN" ){
        OptRealm = EGOR_KIN;
    } else if( Options.GetOptOptRealm() == "TOT" ){
        OptRealm = EGOR_TOT;
    } else {
        RUNTIME_ERROR("unsupported opt realm");
    }

    switch(OptRealm){
        case(EGOR_CVS):
            InData = InAccu->GetSectionData("TCVS")->GetDataBlob(InCVS);
            break;
        case(EGOR_ICF):
            RUNTIME_ERROR("not implemented");
            break;
        case(EGOR_KIN):
            InData = InAccu->GetSectionData("TEKIN")->GetDataBlob();
            break;
        case(EGOR_TOT):
            GetEtot();
            InData = InAccu->GetSectionData("TETOT")->GetDataBlob();
            break;
    }

    GPRLen = Options.GetOptGPRLen();

// -----------------------------------------------------------------------------
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFGPROptHyprms::GetEtot(void)
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

void CABFGPROptHyprms::InitOptimizer(void)
{
    vout << endl;
    vout << format("%02d:Optimization of ABF GPR hyperparameters ...")%State << endl;
    State++;
    CGPFilter gpr;
    gpr.PrintExecInfo(vout);

// what should be optimized?
    WidthEnabled = Options.GetOptWidthEnabled();
    NoiseEnabled = Options.GetOptNoiseEnabled();

// number of optimized parameters
    NumOfOptPrms = 0;
    NumOfPrms = 0;
// width
    if( WidthEnabled ) NumOfOptPrms++;
    NumOfPrms++;
// noise
    if( NoiseEnabled ) NumOfOptPrms++;
    NumOfPrms++;

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
    vout << "   Optimized realm                 = " << Options.GetOptOptRealm() << endl;
    vout << "   GPR kernel                      = " << Options.GetOptGPRKernel() << endl;

// initial values
    Width = Options.GetOptWidth();
    Noise = Options.GetOptNoise();

    // print hyperparameters
    vout << "   Hyperparameters ..." << endl;
    vout << format("      Width  = %16.4f")%Width << endl;
    vout << format("      Noise  = %16.4f")%Noise << endl;

    if( ! Options.GetOptTest() ){
        vout << endl;
        switch(OptTarget){
            case(EGOT_LOGML):
                vout << "# step         logML";
            break;
            case(EGOT_LOGPL):
                vout << "# step         logPL";
            break;
        }

        if( WidthEnabled )    vout << "      Width";
        if( NoiseEnabled )    vout << "      Noise";
        vout << "      KRank";
        vout << endl;

        vout << "# ---- -------------";
        if( WidthEnabled )    vout << " ----------";
        if( NoiseEnabled )    vout << " ----------";
        vout << " ----------";
        vout << endl;
    }

    Hyprms.CreateVector(NumOfOptPrms);
    HyprmsGrd.CreateVector(NumOfOptPrms);


// set initial hyprms
    int ind = 0;
    if( WidthEnabled ){
        Hyprms[ind] = Width;
        ind++;
    }
    if( NoiseEnabled ){
        Hyprms[ind] = Noise;
        ind++;
    }

// init L-BFGS
    // N*(2*M+1)+2*M
    Work.CreateVector(NumOfOptPrms*(2*nlbfgscorr+1)+2*nlbfgscorr);
    Work.SetZero();
    TmpXG.CreateVector(NumOfOptPrms);
    TmpXG.SetZero();
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

bool CABFGPROptHyprms::Optimize(void)
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

//    // input parameters cannot be zero or negative
//    for(int i=0; i < NumOfOptPrms; i++){
//        if( Hyprms[i] <= 0.0 ){
//            Hyprms[i] = 1.0;
//        }
//    }

    for(int k=1; k <= noptsteps; k++ ){
        // parameters cannot be zero or negative
        bool reset = false;
//        for(int i=0; i < NumOfOptPrms; i++){
//            if( Hyprms[i] <= 0.0 ){
//                reset = true;
//            }
//        }
        if( reset ){
            vout << endl;
            vout << "<b><blue>>>> INFO: Parameters out-of-range, resetting ...</blue></b>" << endl;
            Hyprms = OldHyprms;
            if( ResetOpt(numofreset) == false ){
                result = false;
                break;
            }
            vout << endl;
            iflag = 0;  // reset optimizer
            last_logtrg = last_logtrg - 10; // be sure that the number is somehow different
        }

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
        if( WidthEnabled )    vout << " ----------";
        if( NoiseEnabled )    vout << " ----------";
        vout << " ----------";
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

bool CABFGPROptHyprms::ResetOpt(int& numofreset)
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

double CABFGPROptHyprms::GetGNorm(void)
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

void CABFGPROptHyprms::Test(void)
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

    noffset = soffset+1;
    woffset = noffset+1;

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
                vout << "NCorr      ";
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

void CABFGPROptHyprms::CalcHessian(void)
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
    int    krank;

    for(int i=0; i < NumOfOptPrms; i++){

        vout << "   + perturbation for hyprm: " << setw(3) << (i+1) << endl;
        grd1[i].SetZero();
        tmp_prms = Hyprms;
        tmp_prms[i] += dh;
        ScatterHyprms(tmp_prms);

        {
            // calculate gradient and logML
            CGPFilter gpr;
            GetTargetFromGPFilter(gpr,HyprmsEnabled,grd1[i],krank);
        }

        vout << "   - perturbation for hyprm: " << setw(3) << (i+1) << endl;
        grd2[i].SetZero();
        tmp_prms = Hyprms;
        tmp_prms[i] -= dh;
        ScatterHyprms(tmp_prms);

        {
            // calculate gradient and logPL
            CGPFilter gpr;
            GetTargetFromGPFilter(gpr,HyprmsEnabled,grd2[i],krank);
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

    int soffset = 0;
    int noffset;
    int woffset;

    noffset = soffset+1;
    woffset = noffset+1;

    int ind = 0;
    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( HyprmsEnabled[prm] ){
            if( prm == 0 ){
                vout << " SigmaF2   ";
            } else if( (prm >= noffset) && (prm < woffset) ){
                vout << " NCorr     ";
            } else {
                int cv = prm - woffset;
                vout << format(" WFac#%-2d   ")%(cv+1);
            }
            ind++;
        }
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

void CABFGPROptHyprms::PrintGradientSummary(void)
{
    vout << endl;
    vout << "# Final results ..." << endl;
    vout << "# idx       prms          value       gradient" << endl;
    vout << "# --- ---------- -------------- --------------" << endl;

    double gnorm = 0.0;

    int ind = 0;
    for(int prm=0; prm < (int)HyprmsEnabled.size(); prm++){
        if( HyprmsEnabled[prm] ){
            vout << format("%5d ")%(ind+1);
            if( prm == 0 ){
                vout << "Width      ";
            } else {
                vout << "Noise      ";
            }
            gnorm += HyprmsGrd[ind]*HyprmsGrd[ind];
            vout << format("%14.6e %14.6e")%Hyprms[ind]%HyprmsGrd[ind] << endl;
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

void CABFGPROptHyprms::RunGPRAnalytical(void)
{
// setup parameters
    ScatterHyprms(Hyprms);
    HyprmsGrd.SetZero();

    // calculate gradient and logML/PL
    CGPFilter gpr;
    logTarget = GetTargetFromGPFilter(gpr,HyprmsEnabled,HyprmsGrd,KRank);
}

//------------------------------------------------------------------------------

void CABFGPROptHyprms::RunGPRNumerical(void)
{
    CSimpleVector<double>   tmp_prms;
    tmp_prms.CreateVector(NumOfOptPrms);

    CSimpleVector<double>   tmp_grd;
    tmp_grd.CreateVector(NumOfOptPrms);

    logTarget = 0.0;
    HyprmsGrd.SetZero();

// value
    tmp_prms = Hyprms;
    ScatterHyprms(tmp_prms);
    {
        CGPFilter gpr;
        logTarget = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,KRank);
    }

// derivatives
    double dh = 1e-3;
    int    krank;

    for(int i=0; i < NumOfOptPrms; i++){
        if( Options.GetOptCD5() ){
            double v1,v2,v3,v4;
            tmp_prms = Hyprms;
            tmp_prms[i] -= 2.0*dh;
            ScatterHyprms(tmp_prms);
            {
                CGPFilter gpr;
                v1 = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,krank);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CGPFilter gpr;
                v2 = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,krank);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CGPFilter gpr;
                v3 = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,krank);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += 2.0*dh;
            ScatterHyprms(tmp_prms);
            {
                CGPFilter gpr;
                v4 = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,krank);
            }
            HyprmsGrd[i] += (v1-8.0*v2+8.0*v3-v4)/(12.0*dh);
        } else {
            double lv,rv;
            tmp_prms = Hyprms;
            tmp_prms[i] -= dh;
            ScatterHyprms(tmp_prms);
            {
                CGPFilter gpr;
                lv = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,krank);
            }
            tmp_prms = Hyprms;
            tmp_prms[i] += dh;
            ScatterHyprms(tmp_prms);
            {
                CGPFilter gpr;
                rv = GetTargetFromGPFilter(gpr,HyprmsEnabled,tmp_grd,krank);
            }
            HyprmsGrd[i] += (rv-lv)/(2.0*dh);
        }

    }

    ScatterHyprms(Hyprms);
}

//------------------------------------------------------------------------------

void CABFGPROptHyprms::ScatterHyprms(CSimpleVector<double>& hyprsm)
{
    // update parameters
    HyprmsEnabled.resize(NumOfPrms);

    int ind = 0;
    int i = 0;
// ---------------
    if( WidthEnabled ){
        Width = hyprsm[ind];
        ind++;
        HyprmsEnabled[i] = true;
    } else {
        HyprmsEnabled[i] = false;
    }
    i++;
// ---------------
    if( NoiseEnabled ){
        Noise = hyprsm[ind];
        ind++;
        HyprmsEnabled[i] = true;
    } else {
        HyprmsEnabled[i] = false;
    }
}

//------------------------------------------------------------------------------

double CABFGPROptHyprms::GetTargetFromGPFilter(CGPFilter& gpr,std::vector<bool>& hypen,CSimpleVector<double>& grd,int& irank)
{
// setup filter
    gpr.SetFilter(TimeStep,GPRLen);
    gpr.SetRCond(Options.GetOptRCond());
    gpr.SetKernel(Options.GetOptGPRKernel());
    gpr.SetWidth(Width);
    gpr.SetNoise(Noise); // pow(10.0,-Noise));

    vout << high;
    gpr.PrepareProcess(vout);
    irank = gpr.GetKRank();

    double target = 0.0;
    double nblocks = 0.0;
    grd.SetZero();

    for(size_t t=GPRLen; t < NSTLimit-GPRLen; t += GPRLen){

        // train model
        gpr.TrainProcess(InData,t);

        // get data
        switch(OptTarget){
            case(EGOT_LOGML):
                target += gpr.GetLogML();
                gpr.GetLogMLDerivatives(hypen,grd);
                nblocks++;
            break;
            case(EGOT_LOGPL):
                target += gpr.GetLogPL();
                gpr.GetLogPLDerivatives(hypen,grd);
                nblocks++;
            break;
        }
    }

    if( nblocks > 0 ){
        target /= nblocks;
        for(int i=0; i < NumOfOptPrms; i++){
            grd[i] /= nblocks;
        }
    }

    switch(OptTarget){
        case(EGOT_LOGML):
            vout << "      logML     = " << setprecision(5) << target << endl;
        break;
        case(EGOT_LOGPL):
            vout << "      logPL     = " << setprecision(5) << target << endl;
        break;
    }

    vout << low;

    return(target);
}

//------------------------------------------------------------------------------

void CABFGPROptHyprms::WriteResults(int istep)
{
    vout << format("%6d %13e")%istep%logTarget;
    if( WidthEnabled ) vout << format(" %10.3f")%Width;
    if( NoiseEnabled ) vout << format(" %10.3f")%Noise;
    vout << format(" %10d")%KRank;
    vout << endl;
}

//------------------------------------------------------------------------------

bool CABFGPROptHyprms::WriteHyperPrms(FILE* p_fout)
{
    if( p_fout == NULL ){
        return(false);
    }

    ScatterHyprms(Hyprms);

    if( fprintf(p_fout,"# GPR hyper-parameters\n") <= 0 ) return(false);
    if( fprintf(p_fout,"Width  = %10.4f\n",Width) <= 0 ) return(false);
    if( fprintf(p_fout,"Noise  = %16.4f\n",Noise) <= 0 ) return(false);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFGPROptHyprms::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-gpr-opthyprms terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

