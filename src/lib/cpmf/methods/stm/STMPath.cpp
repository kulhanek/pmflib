// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2025,2026 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
// ===============================================================================

#include <sstream>
#include <STMPath.hpp>
#include <TemplIterator.hpp>
#include <ErrorSystem.hpp>
#include <iomanip>
#include <string>
#include <iterator>
#include <math.h>
#include <XMLIterator.hpp>
#include <PrmUtils.hpp>
#include <CVSplineInterpolatingCubic.hpp>
#include <CVSplineSmoothingCubic.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

/*
Methods:
* GD            - gradient descent (gradient)
* NGD           - normalized gradient descent (normalized gradient)
* NGD-AUTO      - gradient descent (switch between GD and NGD)
* ADAM          - Adaptive Moment Estimation
* AMSGrad       - AMSGrad
* AMSGradBC     - AMSGrad + bias corrected estimates
* AdaBelief     - Adam-Belief
*/

//------------------------------------------------------------------------------

// https://en.wikipedia.org/wiki/Stochastic_gradient_descent << ADAM
// https://www.ruder.io/optimizing-gradient-descent/

// https://en.wikipedia.org/wiki/Barzilai-Borwein_method
// Barzilai-Borwein method does not work
// because gradients are too noisy

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSTMPath::CSTMPath(void)
{
    // files
    InputPath = "{PATHS}";
    OutputPath = "_stm.path";
    PathTrajectory = "_stm.traj";

    OutputPathSummary = "_stm.results";
    OptimizationLog = "_stm.log";

    // intervals
    SmoothInterval  = 0;
    ReparamInterval = 1;

    TrajInterval    = 0;
    OutInterval     = 100;

    // stm
    InitPeriod = 10000;         // initialization period
    SoloInitPeriod = true;

    EquiPeriod =  1000;         // equilibration period
    SoloEquiPeriod = true;

    AccuPeriod =  5000;         // accumulation period
    ProdPeriod = 50000;         // final production period

    MaxSTMSteps         = 100;
    OptMethod           = "adabelief";
    ShifGlobalMinA2Zero = true;

    MaxGNormForGD       = 5.0;
    MinGNormEps         = 1e-7;
    StepSize            = 0.003;
    AdamB1              = 0.7;
    AdamB2              = 0.99;
    ResetAdamAlg        = 0;
    MemoryLength        = 0;

    SmoothingFac        = 0.0;

    AsynchronousMode = false;    // update per bead or path

    STMStep         = 0;

    FinalPLenChange     = 0.001;
    FinalMaxBeadMove    = 0.005;
    FinalAveBeadMove    = 0.005;
    FinalMaxpMFSize     = 10.00;
    FinalAvepMFSize     = 2.00;

    PLenChange = 0;
    MaxBeadMove = 0;
    MaxBeadMoveID = 0;
    AveBeadMove = 0;
    MaxpMFSize = 0;
    MaxpMFSizeID = 0;
    AvepMFSize = 0;

    MABufLength = 3;    // buffers are allocated later

    // control
    NumOfRendezvousBeads = 0;
    STMStatus = ESTMS_INITIALIZED;
    HeaderPrinted = false;
    Terminate = false;

    CVSplineType = "smoothing-cubic";

    // how many points are used to calculate path segment length
    SegmentDiscretization = 10;

    ClearPath();
}

//------------------------------------------------------------------------------

void CSTMPath::AllocatePath(void)
{
    if( NumOfCVs < 2 ){
        RUNTIME_ERROR("number of CVs must be larger than or equal to 2");
    }
    if( NumOfBeads < 3 ){
        RUNTIME_ERROR("number of beads must be larger than or equal to 3");
    }

// CVs
    for(int i=0; i < NumOfCVs; i++){
        CColVariablePtr cv = CColVariablePtr(new CColVariable);
        CVs.push_back(cv);
    }

// path beads
    for(int b=0; b < NumOfBeads; b++){
        CBeadPtr bead = CBeadPtr(new CBead);
        Beads.push_back(bead);
        bead->InitBead(this,NumOfCVs);
    }
}

//------------------------------------------------------------------------------

void CSTMPath::ClearPath(void)
{
    NumOfBeads = 0;
    NumOfCVs = 0;
    CVs.clear();
    Beads.clear();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMPath::AttachVerboseStream(std::ostream& str,bool verbose)
{
    vout.Attach(str);
    if( verbose ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessFilesControl(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [files] ====================================================================" << endl;
    if(prmfile.OpenSection("files") == false) {
        vout << "Input path (input)                             = " << setw(20) << InputPath
             << "  (default)" << endl;
        vout << "Output path (output)                           = " << setw(20) << OutputPath
             << "  (default)" << endl;
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory
             << "  (default)" << endl;
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary
             << "  (default)" << endl;
        vout << "Optimization journal (optlog)                  = " << setw(20) << OptimizationLog
             << "  (default)" << endl;
        return(true);
    }

    if(prmfile.GetStringByKey("input",InputPath) == true) {
        vout << "Input path (input)                             = " << setw(20) << InputPath << endl;
    } else {
        vout << "Input path (input)                             = " << setw(20) << InputPath
             << "  (default)" << endl;
    }

    if(prmfile.GetStringByKey("output",OutputPath) == true) {
        vout << "Output path (output)                           = " << setw(20) << OutputPath << endl;
    } else {
        vout << "Output path (output)                           = " << setw(20) << OutputPath
             << "  (default)" << endl;
    }

    if(prmfile.GetStringByKey("trajectory",PathTrajectory) == true) {
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory << endl;
    } else {
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory
             << "  (default)" << endl;
    }

    if(prmfile.GetStringByKey("summary",OutputPathSummary) == true) {
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary << endl;
    } else {
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary
             << "  (default)" << endl;
    }

    if(prmfile.GetStringByKey("optlog",OptimizationLog) == true) {
        vout << "Optimization journal (optlog)                  = " << setw(20) << OptimizationLog << endl;
    } else {
        vout << "Optimization journal (optlog)                  = " << setw(20) << OptimizationLog
             << "  (default)" << endl;
    }


    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessSTMControl(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [stm] ======================================================================" << endl;
    if(prmfile.OpenSection("stm") == false) {
        vout << "Max number of STM steps (steps)                = " << setw(12) << right << MaxSTMSteps
            << "          (default)" << endl;
        vout << "Optimization method (optmethod)                = " << setw(12) << right << OptMethod
            << "          (default)" << endl;
        vout << "Shift global A min to zero (shift2zero)        = " << setw(12) << right << PrmFileOnOff(ShifGlobalMinA2Zero)
            << "          (default)" << endl;

        vout << "Initialization period (init)                   = " << setw(12) << right << InitPeriod
            << "          (default)" << endl;
        vout << "Solo initialization period (soloinit)          = " << setw(12) << right << PrmFileOnOff(SoloInitPeriod)
            << "          (default)" << endl;
        vout << "Accumulation period (accu)                     = " << setw(12) << right << AccuPeriod
            << "          (default)" << endl;
        vout << "Equilibration period (equi)                    = " << setw(12) << right << EquiPeriod
            << "          (default)" << endl;
        vout << "Solo equilibration period (soloequi)           = " << setw(12) << right << PrmFileOnOff(SoloEquiPeriod)
            << "          (default)" << endl;
        vout << "Final production period (prod)                 = " << setw(12) << right << ProdPeriod
            << "          (default)" << endl;

        vout << "Path smoothing factor (sfac)                   = " << setw(12) << right << SmoothingFac
            << "          (default)" << endl;
        vout << "Asynchronous mode (async)                      = " << setw(12) << right << PrmFileOnOff(AsynchronousMode)
            << "          (default)" << endl;

        // opt method
        bool result  = ProcessGDOptMethodSetup(prmfile);
        return(result);
    }

    if(prmfile.GetIntegerByKey("steps",MaxSTMSteps) == true) {
        vout << "Max number of STM steps (steps)                = " << setw(12) << right << MaxSTMSteps << endl;
    } else {
        vout << "Max number of STM steps (steps)                = " << setw(12) << right << MaxSTMSteps
            << "          (default)" << endl;
    }

    if(prmfile.GetStringByKey("optmethod",OptMethod) == true) {
        vout << "Optimization method (optmethod)                = " << setw(12) << right << OptMethod << endl;
    } else {
        vout << "Optimization method (optmethod)                = " << setw(12) << right << OptMethod
            << "          (default)" << endl;
    }

    if(prmfile.GetLogicalByKey("shift2zero",ShifGlobalMinA2Zero) == true) {
        vout << "Shift global A min to zero (shift2zero)        = " << setw(12) << right << PrmFileOnOff(ShifGlobalMinA2Zero) << left << endl;
    } else {
        vout << "Shift global A min to zero (shift2zero)        = " << setw(12) << right << PrmFileOnOff(ShifGlobalMinA2Zero)
            << "          (default)" << endl;
    }

    
    if(prmfile.GetIntegerByKey("init",InitPeriod) == true) {
        vout << "Initialization period (init)                   = " << setw(12) << right << InitPeriod << endl;
    } else {
        vout << "Initialization period (init)                   = " << setw(12) << right << InitPeriod
            << "          (default)" << endl;
    }

    if(prmfile.GetLogicalByKey("soloinit",SoloInitPeriod) == true) {
        vout << "Solo initialization period (soloinit)          = " << setw(12) << right << PrmFileOnOff(SoloInitPeriod) << left << endl;
    } else {
        vout << "Solo initialization period (soloinit)          = " << setw(12) << right << PrmFileOnOff(SoloInitPeriod)
            << "          (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("accu",AccuPeriod) == true) {
        vout << "Accumulation period (accu)                     = " << setw(12) << right << AccuPeriod << endl;
    } else {
        vout << "Accumulation period (accu)                     = " << setw(12) << right << AccuPeriod
            << "          (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("equi",EquiPeriod) == true) {
        vout << "Equilibration period (equi)                    = " << setw(12) << right << EquiPeriod << endl;
    } else {
        vout << "Equilibration period (equi)                    = " << setw(12) << right << EquiPeriod
            << "          (default)" << endl;
    }
    if(prmfile.GetLogicalByKey("soloequi",SoloEquiPeriod) == true) {
        vout << "Solo equilibration period (soloequi)           = " << setw(12) << right << PrmFileOnOff(SoloEquiPeriod) << left << endl;
    } else {
        vout << "Solo equilibration period (soloequi)           = " << setw(12) << right << PrmFileOnOff(SoloEquiPeriod)
            << "          (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("prod",ProdPeriod) == true) {
        vout << "Final production period (prod)                 = " << setw(12) << right << ProdPeriod << endl;
    } else {
        vout << "Final production period (prod)                 = " << setw(12) << right << ProdPeriod
            << "          (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("sfac",SmoothingFac) == true) {
        vout << "Path smoothing factor (sfac)                   = " << setw(12) << right << SmoothingFac << left << endl;
    } else {
        vout << "Path smoothing factor (sfac)                   = " << setw(12) << right << SmoothingFac
            << "          (default)" << endl;
    }

    if(prmfile.GetLogicalByKey("async",AsynchronousMode) == true) {
        vout << "Asynchronous mode (async)                      = " << setw(12) << right << PrmFileOnOff(AsynchronousMode) << left << endl;
    } else {
        vout << "Asynchronous mode (async)                      = " << setw(12) << right << PrmFileOnOff(AsynchronousMode)
            << "          (default)" << endl;
    }

// optimization method setup
    OptMethod.ToLowerCase();

    bool result = true;
    if( OptMethod == "gd" ){
        result = ProcessGDOptMethodSetup(prmfile);
    } else if( OptMethod == "ngd" ){
        result = ProcessNGDOptMethodSetup(prmfile);
    } else if( OptMethod == "ngd-auto" ){
        result = ProcessNGDAutoOptMethodSetup(prmfile);
    } else if( OptMethod == "adam" ){
        result = ProcessAdamOptMethodSetup(prmfile);
    } else if( OptMethod == "adabelief" ){
        result = ProcessADABeliefOptMethodSetup(prmfile);
    } else if( OptMethod == "amsgrad" ){
        result = ProcessAMSGradOptMethodSetup(prmfile);
    } else if( OptMethod == "amsgradbc" ){
        result = ProcessAMSGradBCOptMethodSetup(prmfile);
    } else {
        RUNTIME_ERROR("not implemented opt method");
    }

    return(result);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessSTMTerminationControl(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [termination] ==============================================================" << endl;
    if(prmfile.OpenSection("termination") == false) {

        vout << "Final path length change (fplch)               = " << setw(9) << FinalPLenChange
             << left << "             (default)" << endl;
        vout << "Max final bead movement (maxbmove)             = " << setw(9) << FinalMaxBeadMove
             << left << "             (default)" << endl;
        vout << "Average final bead movement (avebmove)         = " << setw(9) << FinalAveBeadMove
             << left << "             (default)" << endl;

        vout << "Max perpendicular mean force (maxppmf)         = " << setw(9) << FinalMaxpMFSize
             << left << "             (default)" << endl;
        vout << "Average perpendicular mean force (aveppmf)     = " << setw(9) << FinalAvepMFSize
             << left << "             (default)" << endl;

        vout << "Term buffer length (tbuflen)                   = " << setw(9) << MABufLength
             << left << "             (default)" << endl;
        
        return(true);
    }

    if(prmfile.GetDoubleByKey("fplch",FinalPLenChange) == true) {
        vout << "Final path length change (fplch)               = " << setw(9) << FinalPLenChange << left << endl;
    } else {
        vout << "Final path length change (fplch)               = " << setw(9) << FinalPLenChange
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("maxbmove",FinalMaxBeadMove) == true) {
        vout << "Max final bead movement (maxbmove)             = " << setw(9) << FinalMaxBeadMove << left << endl;
    } else {
        vout << "Max final bead movement (maxbmove)             = " << setw(9) << FinalMaxBeadMove
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("avebmove",FinalAveBeadMove) == true) {
        vout << "Average final bead movement (avebmove)         = " << setw(9) << FinalAveBeadMove << left << endl;
    } else {
        vout << "Average final bead movement (avebmove)         = " << setw(9) << FinalAveBeadMove
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("",FinalMaxpMFSize) == true) {
        vout << "Max perpendicular mean force (maxppmf)         = " << setw(9) << FinalMaxpMFSize << left << endl;
    } else {
        vout << "Max perpendicular mean force (maxppmf)         = " << setw(9) << FinalMaxpMFSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("aveppmf",FinalAvepMFSize) == true) {
        vout << "Average perpendicular mean force (aveppmf)     = " << setw(9) << FinalAvepMFSize << left << endl;
    } else {
        vout << "Average perpendicular mean force (aveppmf)     = " << setw(9) << FinalAvepMFSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("tbuflen",MABufLength) == true) {
        vout << "Term buffer length (tbuflen)                   = " << setw(9) << MABufLength << endl;
    } else {
        vout << "Term buffer length (tbuflen)                   = " << setw(9) << MABufLength
             << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessGDOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [gd] =======================================================================" << endl;
    if(prmfile.OpenSection("gd") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessNGDOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [ngd] ======================================================================" << endl;
    if(prmfile.OpenSection("ngd") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("mingnormeps",MinGNormEps) == true) {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps << left << endl;
    } else {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessNGDAutoOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [ngd] ======================================================================" << endl;
    if(prmfile.OpenSection("ngd") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "Max gnorm to switch to GD (maxgnormforgd)      = " << setw(9) << MaxGNormForGD
             << left << "             (default)" << endl;
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("maxgnormforgd",MaxGNormForGD) == true) {
        vout << "Max gnorm to switch to GD (maxgnormforgd)      = " << setw(9) << MaxGNormForGD << left << endl;
    } else {
        vout << "Max gnorm to switch to GD (maxgnormforgd)      = " << setw(9) << MaxGNormForGD
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("mingnormeps",MinGNormEps) == true) {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps << left << endl;
    } else {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessAdamOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [adam] =====================================================================" << endl;
    if(prmfile.OpenSection("adam") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta1",AdamB1) == true) {
        vout << "beta1                                          = " << setw(9) << AdamB1 << left << endl;
    } else {
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta2",AdamB2) == true) {
        vout << "beta2                                          = " << setw(9) << AdamB2 << left << endl;
    } else {
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("mingnormeps",MinGNormEps) == true) {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps << left << endl;
    } else {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("maxresets",ResetAdamAlg) == true) {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg << left << endl;
    } else {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("memsteps",MemoryLength) == true) {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength << left << endl;
    } else {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessADABeliefOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [adabelief] ================================================================" << endl;

    if(prmfile.OpenSection("adabelief") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta1",AdamB1) == true) {
        vout << "beta1                                          = " << setw(9) << AdamB1 << left << endl;
    } else {
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta2",AdamB2) == true) {
        vout << "beta2                                          = " << setw(9) << AdamB2 << left << endl;
    } else {
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("mingnormeps",MinGNormEps) == true) {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps << left << endl;
    } else {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("maxresets",ResetAdamAlg) == true) {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg << left << endl;
    } else {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("memsteps",MemoryLength) == true) {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength << left << endl;
    } else {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessAMSGradOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [amsgrad] ==================================================================" << endl;
    if(prmfile.OpenSection("amsgrad") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg
             << left << "             (default)" << endl;
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta1",AdamB1) == true) {
        vout << "beta1                                          = " << setw(9) << AdamB1 << left << endl;
    } else {
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta2",AdamB2) == true) {
        vout << "beta2                                          = " << setw(9) << AdamB2 << left << endl;
    } else {
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("mingnormeps",MinGNormEps) == true) {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps << left << endl;
    } else {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("maxresets",ResetAdamAlg) == true) {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg << left << endl;
    } else {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("memsteps",MemoryLength) == true) {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength << left << endl;
    } else {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessAMSGradBCOptMethodSetup(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [amsgradbc] ================================================================" << endl;
    if(prmfile.OpenSection("amsgradbc") == false) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg
             << left << "             (default)" << endl;
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
        return(true);
    }

    if(prmfile.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta1",AdamB1) == true) {
        vout << "beta1                                          = " << setw(9) << AdamB1 << left << endl;
    } else {
        vout << "beta1                                          = " << setw(9) << AdamB1
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("beta2",AdamB2) == true) {
        vout << "beta2                                          = " << setw(9) << AdamB2 << left << endl;
    } else {
        vout << "beta2                                          = " << setw(9) << AdamB2
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("mingnormeps",MinGNormEps) == true) {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps << left << endl;
    } else {
        vout << "Min gnorm value (mingnormesp)                  = " << setw(9) << MinGNormEps
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("maxresets",ResetAdamAlg) == true) {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg << left << endl;
    } else {
        vout << "Max resets (maxresets)                         = " << setw(9) << ResetAdamAlg
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("memsteps",MemoryLength) == true) {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength << left << endl;
    } else {
        vout << "Reset memory every (memsteps)                  = " << setw(9) << MemoryLength
             << left << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessIntervalsControl(CPrmFile& prmfile)
{
    vout << endl;
    vout << "=== [intervals] ================================================================" << endl;
    if(prmfile.OpenSection("intervals") == false) {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval
             << "             (default)" << endl;
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval
             << "             (default)" << endl;
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval
             << "             (default)" << endl;
        vout << "Output path update (output)                    = " << setw(9) << OutInterval
             << "             (default)" << endl;

        return(true);
    }

    if(prmfile.GetIntegerByKey("reparam",ReparamInterval) == true) {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval << endl;
    } else {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("smooth",SmoothInterval) == true) {
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval << endl;
    } else {
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("trajectory",TrajInterval) == true) {
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval << endl;
    } else {
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("output",OutInterval) == true) {
        vout << "Output path update (output)                    = " << setw(9) << OutInterval << endl;
    } else {
        vout << "Output path update (output)                    = " << setw(9) << OutInterval
             << "             (default)" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::ProcessPathControl(CPrmFile& prmfile)
{
//    ! [PATH]
//    ! nbeads     number_of_beads
//    ! ncvs       number_of_cvs
//    ! names      cv1 cv2 cv3 ... cvn
//    ! types      cv1 cv2 cv3 ... cvn
//    ! min        cv1 cv2 cv3 ... cvn
//    ! max        cv1 cv2 cv3 ... cvn
//    ! maxmov     cv1 cv2 cv3 ... cvn
//    ! normal     cv1 cv2 cv3 ... cvn
//    ! permanent  cv1 cv2 cv3 ... cvn
//    ! free       cv1 cv2 cv3 ... cvn

    // clear path
    ClearPath();

    vout << endl;
    vout << "=== [PATH] =====================================================================" << endl;
    if(prmfile.OpenSection("PATH") == false) {
        RUNTIME_ERROR("[PATH] section is mandatory for a path specification");
    }

    if(prmfile.GetStringByKey("name",PathName) == true) {
        vout << "Path name                         = " << PathName << endl;
    } else {
        RUNTIME_ERROR("path name (name) is not specified");
    }

    // read number of beads and CVs
    if(prmfile.GetIntegerByKey("nbeads",NumOfBeads) == true) {
        vout << "Number of beads (nbeads)          = " << NumOfBeads << endl;
    } else {
        RUNTIME_ERROR("number of beads (nbeads) is not specified");
    }
    if(prmfile.GetIntegerByKey("ncvs",NumOfCVs) == true) {
        vout << "Number of CVS (ncvs)              = " << NumOfCVs << endl;
    } else {
        RUNTIME_ERROR("number of CVs (ncvs) is not specified");
    }

// allocate path -------------
    AllocatePath();

// print header --------------
    // legends
    vout << endl;
    vout << " ID   Type ";
    for(int i=0; i < NumOfCVs; i++){
        vout << "     CV" << left << setw(2) << i+1 << "    ";
    }
    vout << endl;

    // delimiters
    vout << "---- ------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

    // load cvs,types,min and max items
    ReadPathControls(prmfile);

    // read temporary path - only points specified by user
    int num_of_user_beads = ReadPathNumberOfUserBeads(prmfile);

    vout << debug << "Number of user provided beads: " << num_of_user_beads << endl;
    vout << high;

    for(int i=0; i < num_of_user_beads; i++){
        CBeadPtr bead = CBeadPtr(new CBead);
        InputBeads.push_back(bead);
    }

    ReadPathUserBeads(prmfile,InputBeads);

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::LoadCVSplines(CPrmFile& prmfile)
{
    if( NumOfCVs < 2 ){
        RUNTIME_ERROR("number of CVs must be larger than or equal to 2");
    }

    if( prmfile.OpenGroup("CVSPLINES") == false ) {
        vout << ">> Info: No {CVSPLINES} group is specified - using the default interpolating cubic splines ..." << endl;
        CVSplineType = "interpolating-cubic";
    } else {
        vout << endl;
        vout << "=== [setup] ====================================================================" << endl;
        if( prmfile.OpenSection("setup") == false ) {
            vout << "CV spline type (type)                          = " << left << setw(20) << CVSplineType << "  (default)" << endl;
        } else {
            if( prmfile.GetStringByKey("type",CVSplineType) == true  ) {
                vout << "CV spline type (type)                          = " << left << setw(20) << CVSplineType << endl;
            } else {
                vout << "CV spline type (type)                          = " << left << setw(20) << CVSplineType << "  (default)" << endl;
            }
        }
    }

// path splines
    for(int i=0; i < NumOfCVs; i++){
        vout << endl;
        vout << "# ### CV: " << CVs[i]->GetName() << endl;
        vout << "# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

        CCVSplinePtr cvspline;
        if( CVSplineType ==  "interpolating-cubic" ){
            cvspline = CCVSplinePtr(new CCVSplineInterpolatingCubic);
        } else if( CVSplineType ==  "smoothing-cubic" ){
            cvspline = CCVSplinePtr(new CCVSplineSmoothingCubic);
        } else {
            RUNTIME_ERROR("not implemented CV spline");
        }

        // search CV spline setup
        bool found = false;
        bool res = prmfile.FirstSection();
        while( res == true ){
            if( prmfile.GetSectionName() == "CV" ){
                CSmallString cv_name;
                if( prmfile.GetStringByKey("name",cv_name) && (cv_name == CVs[i]->GetName()) ){
                    found = true;
                    cvspline->LoadSetup(prmfile,vout);
                }
            }
            res = prmfile.NextSection();
        }
        if( found == false ){
            if( prmfile.OpenSection("default") == true ) {
                cvspline->LoadSetup(prmfile,vout);
            }
        }

        cvspline->PrintSetup(vout);

        CVSplines.push_back(cvspline);
    }

    vout << endl;
    vout << ":::::::::::::::::::::::::::::::::: | FULL PATH | :::::::::::::::::::::::::::::::" << endl;

// helper vector - spline segment length
    SPos.CreateVector(NumOfCVs);

    int num_of_user_beads = InputBeads.size();

    if( num_of_user_beads != NumOfBeads ){
        // incomplete path

        // optimize path
        for(int b=0; b < num_of_user_beads; b++){
            if( (b != 0) && (b != num_of_user_beads - 1) ){
                if( InputBeads[b]->BeadType != BTY_NORMAL ){
                    RUNTIME_ERROR("the internal beads for incomplete path must be all of \"normal\" type!");
                }
            }
            InputBeads[b]->PPos = InputBeads[b]->Pos;
        }
        OptimizePath(InputBeads);

        // generate missing points or re-optimize path
        Beads[0]->Alpha = 0.0;
        Beads[0]->BeadType = InputBeads[0]->BeadType;
        Beads[0]->BeadID = 1;
        Beads[NumOfBeads-1]->Alpha = 1.0;
        Beads[NumOfBeads-1]->BeadType = InputBeads[num_of_user_beads-1]->BeadType;
        Beads[NumOfBeads-1]->BeadID = NumOfBeads;
        for(int i=0; i < NumOfCVs; i++){
            Beads[0]->Pos[i] = CVSplines[i]->GetCV(0.0);
            Beads[NumOfBeads-1]->Pos[i] = CVSplines[i]->GetCV(1.0);
            for(int b=1; b < NumOfBeads-1; b++){
                double alpha = (double)b / ((double)NumOfBeads-1.0);
                Beads[b]->Pos[i] = CVSplines[i]->GetCV(alpha);
                Beads[b]->Alpha = alpha;
                Beads[b]->BeadID = b + 1;
                Beads[b]->BeadType = BTY_NORMAL;
            }
        }

    } else {
        OptimizePath(InputBeads);

        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->Pos = InputBeads[i]->Pos;
            Beads[i]->Alpha = InputBeads[i]->Alpha;
            Beads[i]->BeadID = InputBeads[i]->BeadID;
            Beads[i]->BeadType = InputBeads[i]->BeadType;
        }
    }

    // check boundaries
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->FPos = Beads[b]->Pos;
    }
    CheckBoundaries();

    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->PPos = Beads[b]->FPos;
        Beads[b]->Pos  = Beads[b]->FPos;
    }
    CurrentPathLength = OptimizePath(Beads);

    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::ReadPathControls(CPrmFile& file)
{
    // load names,types,min and max items

// cvs ---------------------
    CSmallString tmp;
    if(file.GetStringByKey("names",tmp) == true) {
        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( (int)tokens.size() != NumOfCVs ){
            RUNTIME_ERROR("incorrect number of cvs");
        }

        vout << left << "     names " << right;
        for(int i=0; i < NumOfCVs; i++){
            CVs[i]->ID = i;
            CVs[i]->SetName(tokens[i]);
            vout << " " << setw(12) << tokens[i];
        }
        vout << endl;
    }

// types ---------------------
    if(file.GetStringByKey("types",tmp) == true) {
        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( (int)tokens.size() != NumOfCVs ){
            RUNTIME_ERROR("incorrect number of types");
        }

        vout << left << "     types " << right;
        for(int i=0; i < NumOfCVs; i++){
            CVs[i]->SetType(tokens[i]);
            vout << " " << setw(12) << tokens[i];
        }
        vout << endl;
    }

// min values ----------------
    if(file.GetStringByKey("min",tmp) == true) {
        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( (int)tokens.size() != NumOfCVs ){
            RUNTIME_ERROR("incorrect number of min values");
        }

        vout << left << "     min   " << right << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            double min = CSmallString(tokens[i]).ToDouble();
            CVs[i]->SetMinValue(min);
            vout << " " << setw(12) << min;
        }
        vout << endl;
    }

// max values ----------------
    if(file.GetStringByKey("max",tmp) == true) {
        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( (int)tokens.size() != NumOfCVs ){
            RUNTIME_ERROR("incorrect number of max values");
        }

        vout << left << "     max   " << right << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            double max = CSmallString(tokens[i]).ToDouble();
            CVs[i]->SetMaxValue(max);
            vout << " " << setw(12) << max;
        }
        vout << endl;
    }

// max values ----------------
    if(file.GetStringByKey("maxmov",tmp) == true) {
        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( (int)tokens.size() != NumOfCVs ){
            RUNTIME_ERROR("incorrect number of maxmov values");
        }

        vout << left << "     maxmov" << right << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            double maxmov = CSmallString(tokens[i]).ToDouble();
            vout << " " << setw(12) << maxmov;
            maxmov = maxmov / (CVs[i]->GetMaxValue() - CVs[i]->GetMinValue());
            CVs[i]->SetMaxMovement(maxmov); // scaled
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

int CSTMPath::ReadPathNumberOfUserBeads(CPrmFile& file)
{
    CSmallString tmp;

    file.FirstLine();
    int nbeads = 0;

    while( file.GetLine(tmp) ){
        file.NextLine();

        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( tokens.size() >= 2 ){
            if( tokens[0] == "name" ) continue;
            if( tokens[0] == "ncvs" ) continue;
            if( tokens[0] == "nbeads" ) continue;
        }

        if( (int)tokens.size() != (NumOfCVs+1) ){
            CSmallString error;
            error << "incorrect number of items for '" << tokens[0] << "' key";
            RUNTIME_ERROR(error);
        }
        // skip already processed keys
        if( tokens[0] == "names" ) continue;
        if( tokens[0] == "types" ) continue;
        if( tokens[0] == "min" ) continue;
        if( tokens[0] == "max" ) continue;
        if( tokens[0] == "maxmov" ) continue;

        if( (tokens[0] != "normal") && (tokens[0] != "permanent") && (tokens[0] != "free") ){
            CSmallString error;
            error << "unsupported key '" << tokens[0] << "'";
            RUNTIME_ERROR(error)
        }
        if( nbeads > NumOfBeads ){
            RUNTIME_ERROR("more beads specification than nbeads");
        }
        nbeads++;
    }

    return(nbeads);
}

//------------------------------------------------------------------------------

void CSTMPath::ReadPathUserBeads(CPrmFile& file,std::vector<CBeadPtr>& beads)
{
    CSmallString tmp;

    file.FirstLine();
    int beadid = 0;

    while( file.GetLine(tmp) ){
        file.SetCurrentLineProcessed();
        file.NextLine();

        vector<string> tokens;
        SplitString(string(tmp),tokens);

        if( tokens.size() >= 2 ){
            if( tokens[0] == "name" ) continue;
            if( tokens[0] == "ncvs" ) continue;
            if( tokens[0] == "nbeads" ) continue;
        }

        if( (int)tokens.size() != (NumOfCVs+1) ){
            CSmallString error;
            error << "incorrect number of items for '" << tokens[0] << "' key";
            RUNTIME_ERROR(error);
        }
        // skip already processed keys
        if( tokens[0] == "names" ) continue;
        if( tokens[0] == "types" ) continue;
        if( tokens[0] == "min" ) continue;
        if( tokens[0] == "max" ) continue;
        if( tokens[0] == "maxmov" ) continue;

        if( (tokens[0] != "normal") && (tokens[0] != "permanent") && (tokens[0] != "free")  ){
            CSmallString error;
            error << "unsupported key '" << tokens[0] << "'";
            RUNTIME_ERROR(error)
        }
        if( beadid > NumOfBeads ){
            RUNTIME_ERROR("more beads specification than nbeads");
        }

        // process point definition
        beads[beadid]->InitBead(this,NumOfCVs);
        for(int i=0; i < NumOfCVs; i++){
            double unscaled = CSmallString(tokens[i+1]).ToDouble();
            double scaled = CVs[i]->GetScaledValue(unscaled);
            beads[beadid]->Pos[i] = scaled;
        }

        if( tokens[0] == "free" ){
            vout << setw(4) << beadid+1 << " F     " << scientific << setprecision(5);
            beads[beadid]->BeadType = BTY_FREE;
        } else if( tokens[0] == "permanent" ){
            vout << setw(4) << beadid+1 << " P     " << scientific << setprecision(5);
            beads[beadid]->BeadType = BTY_PERMANENT;
        } else if( tokens[0] == "normal" ){
            vout << setw(4) << beadid+1 << " N     " << scientific << setprecision(5);
            beads[beadid]->BeadType = BTY_NORMAL;
        } else {
            CSmallString error;
            error << "unsupported bead type '" << tokens[0] << "'";
            RUNTIME_ERROR(error) 
        }

        for(int i=0; i < NumOfCVs; i++){
            double unscaled = CVs[i]->GetUnscaledValue(beads[beadid]->Pos[i]);
            vout << " " << setw(12) << unscaled;
        }
        vout << endl;
        beadid++;
    }
}

//------------------------------------------------------------------------------

bool CSTMPath::GetSoloInitPeriod(void) const
{
    return(SoloInitPeriod);
}

//------------------------------------------------------------------------------

void CSTMPath::SetSoloInitPeriod(bool set) 
{
    SoloInitPeriod = set;
}

//------------------------------------------------------------------------------

bool CSTMPath::GetSoloEquiPeriod(void) const
{
    return(SoloEquiPeriod);
}

//------------------------------------------------------------------------------

void CSTMPath::SetSoloEquiPeriod(bool set)
{
    SoloEquiPeriod = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSTMPath::LoadPath(CPrmFile& file)
{
    // load path ------------------------------
    if( (InputPath.GetLength() > 0) && (InputPath[0] == '{') ){
        CSmallString grpname;
        grpname = InputPath.GetSubStringFromTo(1,InputPath.GetLength()-2);
        vout << ">> Info: The path is taken from the server control file." << endl;
        if( file.OpenGroup(grpname) == false ){
            CSmallString error;
            error << "unable to open path group '" << InputPath << "'";
            ES_ERROR(error);
            return(false);
        }
        try{
            ProcessPathControl(file);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to load path",e);
            return(false);
        }
    } else {
        vout << ">> Info: The path is taken from the file: " << InputPath << endl;
        CPrmFile path_file;

        if(path_file.Read(InputPath) == false) {
            ES_ERROR("unable to read path file");
            return(false);
        }
        try{
            ProcessPathControl(path_file);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to load path",e);
            return(false);
        }

        if( path_file.CountULines() > 0 ){
            vout << endl;
            ES_ERROR("unprocessed items found in path file");
            path_file.Dump(stderr,true);
            return(false);
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CSTMPath::SavePath(void)
{
    vout << "Output STM path:         " << OutputPath <<  endl;
    try{
        SavePath(OutputPath);
    } catch(...){
        return(false);
    }
    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::SavePath(const CSmallString& name)
{
    ofstream ofs;
    ofs.open(name);
    if( ! ofs ){
        CSmallString error;
        error << "unable to open file '" << name << "'";
        RUNTIME_ERROR(error);
    }
    PrintPath(ofs);
    if( ! ofs ){
        CSmallString error;
        error << "unable to write into file '" << name << "'";
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

bool CSTMPath::SavePathSummary(void)
{
    vout << "Output STM path summary: " << OutputPathSummary <<  endl;
    try{
        SavePathSummary(OutputPathSummary);
    } catch(...){
        return(false);
    }
    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::SavePathSummary(const CSmallString& name)
{
    ofstream ofs;
    ofs.open(name);
    if( ! ofs ){
        CSmallString error;
        error << "unable to open file '" << name << "'";
        RUNTIME_ERROR(error);
    }
    PrintPathSummary(ofs);
    if( ! ofs ){
        CSmallString error;
        error << "unable to write into file '" << name << "'";
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

void CSTMPath::FlushPath(void)
{
    try{
        ProcessingMutex.Lock();
        SavePath(OutputPath);
        SavePathSummary(OutputPathSummary);
        ProcessingMutex.Unlock();
    } catch(...) {
        ProcessingMutex.Unlock();
        throw;
    }
}

//------------------------------------------------------------------------------

bool CSTMPath::OpenTrajectory(void)
{
    if( TrajInterval <= 0 ) return(true);

    TrajectoryFOut.open(PathTrajectory);
    if( ! TrajectoryFOut ){
        ES_ERROR("unable to open path trajectory file");
        return(false);
    }

    // write header
    TrajectoryFOut << "# STMTRAJ " << NumOfCVs << " " << NumOfBeads << endl;
    PrintPathSummaryHeader(TrajectoryFOut);

    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::SaveTrajectorySnapshot(void)
{
    if( TrajInterval <= 0 ) return;
    TrajectoryFOut << "# STMSNAP " << STMStep / TrajInterval << endl;
    PrintPathSummaryData(TrajectoryFOut);
    TrajectoryFOut << endl; // necessary for gnuplot
}

//------------------------------------------------------------------------------

void CSTMPath::CloseTrajectory(void)
{
    TrajectoryFOut.close();
}

//------------------------------------------------------------------------------

bool CSTMPath::OpenOptLog(void)
{
    if( TrajInterval <= 0 ) return(true);

    OptLogFOut.open(OptimizationLog);
    if( ! OptLogFOut ){
        ES_ERROR("unable to open optimization journal");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::CloseOptLog(void)
{
    OptLogFOut.close();
}

//------------------------------------------------------------------------------

void CSTMPath::ForceTermination(void)
{
    ProcessingMutex.Lock();
        Terminate = true;
        vout << ">>> INFO: Received soft termination signal." << endl;
    ProcessingMutex.Unlock();
}

//------------------------------------------------------------------------------

void CSTMPath::SetAsynchronousMode(bool set)
{
    AsynchronousMode = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSTMPath::CheckClient(CXMLElement* p_cele)
{
    try {
        ProcessingMutex.Lock();

    // check CVs
        CXMLElement* p_cvsele = p_cele->GetFirstChildElement("CVS");
        if(p_cvsele == NULL) {
            LOGIC_ERROR("unable to open CVS element");
        }
        if( CheckCoords(p_cvsele) == false) {
            RUNTIME_ERROR("unable to check coordinates");
        }

        ProcessingMutex.Unlock();
    } catch(...) {
        ProcessingMutex.Unlock();
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::RegisterBead(int bead_id,int client_id)
{
    try {
        ProcessingMutex.Lock();

        CBeadPtr p_bead = GetBead(bead_id);
        if( p_bead == NULL ){
            CSmallString error;
            error << "bead with ID=" << bead_id << " not found";
            RUNTIME_ERROR(error);
        }

        if( p_bead->GetClientID() > 0 ){
            CSmallString error;
            error << "bead ID=" << bead_id << " is already registered to client ID=" << p_bead->GetClientID();
            RUNTIME_ERROR(error);
        }

        // assign client
        p_bead->SetClientID(client_id);

        if( HeaderPrinted == false ){
            vout << endl;
            vout << "::::::::::::::::::::::::::::::: Path optimization ::::::::::::::::::::::::::::::" << endl;
            if( AsynchronousMode ){
                vout << "# Entering asynchronous mode ..." << endl;
            } else {
                vout << "# Entering synchronous mode ..." << endl;
                // move to first mode
                for(int i=0; i < NumOfBeads; i++){
                    Beads[i]->MoveToNextMode();
                }
            }
            
            PrintSTMHeader();

            STMStatus = ESTMS_OPTIMIZING;
        }

        ProcessingMutex.Unlock();
    } catch(...) {
        ProcessingMutex.Unlock();
        throw;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::BeginAsynchronousMode(void)
{
    // move to first mode
    for(int i=0; i < NumOfBeads; i++){
        Beads[i]->MoveToNextMode();
    }
}

//------------------------------------------------------------------------------

void CSTMPath::ExchangeData(CXMLElement* p_cele,CXMLElement* p_rele)
{
    if( p_cele == NULL ){
        INVALID_ARGUMENT("p_cele is NULL");
    }
    if( p_rele == NULL ){
        INVALID_ARGUMENT("p_rele is NULL");
    }

    // no exchange
    if( (STMStatus == ESTMS_MAX_STEPS_REACHED) || (STMStatus == ESTMS_COMPLETED) ) return;

    if( AsynchronousMode == true ){
        ExchangeDataAsynchronously(p_cele,p_rele);
    } else {
        ExchangeDataSynchronously(p_cele,p_rele);
    }

}

//------------------------------------------------------------------------------

void CSTMPath::ExchangeDataSynchronously(CXMLElement* p_cele,CXMLElement* p_rele)
{
    if( p_cele == NULL ){
        INVALID_ARGUMENT("p_cele is NULL");
    }
    if( p_rele == NULL ){
        INVALID_ARGUMENT("p_rele is NULL");
    }

    int bead_id =  -1;
    if(p_cele->GetAttribute("bead_id",bead_id) == false) {
        RUNTIME_ERROR("unable to get bead_id");
    }

    // get bead ----------------------------------
    CBeadPtr p_bead = GetBead(bead_id);
    if( p_bead == NULL ){
        CSmallString error;
        error << "unable to find bead ID=" << bead_id;
        RUNTIME_ERROR(error);
    }

    // update bead data --------------------------
    int mode = BMO_UNKNOWN;
    p_cele->GetAttribute("mode",mode);

    if( mode != BMO_UNKNOWN ){
        // do we have correct bead mode?
        if( mode != p_bead->GetMode() ){
            CSmallString error;
            error << "client is in inconsistent mode";
            RUNTIME_ERROR(error);
        }
    }

    if( p_bead->GetModeStatus() == BMS_PREPARED ){
        // first step or bead was released
        p_bead->SetNextStepData(p_rele);
        return;
    }

    if( p_bead->GetModeStatus() != BMS_RUNNING  ){
        CSmallString error;
        error << "bead ID=" << bead_id << " is not in running state - unable to exchange data";
        RUNTIME_ERROR(error);
    }

    if( (p_bead->GetMode() == BMO_ACCUMULATION) || (p_bead->GetMode() == BMO_PRODUCTION) ){
        // only if it is in production state
        p_bead->GetProductionData(p_cele);
        ProcessProductionData(p_bead);
    } else {
        p_bead->SkipProductionData();
    }

    if( (STMStatus == ESTMS_MAX_STEPS_REACHED) || (STMStatus == ESTMS_COMPLETED) ){
        // sent termination status to client
        TerminateClient(p_rele);
        return;
    }

    // update program ----------------------------
    p_bead->MoveToNextMode();

    // set data for client -----------------------
    p_bead->SetNextStepData(p_rele);
}

//------------------------------------------------------------------------------

void CSTMPath::ProcessProductionData(CBeadPtr p_bead)
{
    // how many beads are waiting
    try{
        RendezvousMutex.Lock();

        NumOfRendezvousBeads++;
        if( NumOfRendezvousBeads != NumOfBeads ){
            p_bead->Mode = BMO_WAITFORRENDEZVOUS;
            RendezvousCond.WaitForSignal(RendezvousMutex);
        } else {
            // do all operations on the whole path
            if( STMStatus == ESTMS_PATH_FOUND ){
                CompletePathData();
                IntegratePath();
                STMStep++;
                SavePathAndTraj();

                vout << endl;
                PrintSTMStepInfo();

                // this will happen if the path was found and final production runs are required
                STMStatus = ESTMS_COMPLETED;

                vout << endl;
                vout << ">> INFO: The server is terminated since all data were acquired.*" <<  endl;
            } else {
                try {
                    ProcessingMutex.Lock();
                    switch( p_bead->GetMode() ){
                        case BMO_ACCUMULATION:
                            p_bead->Mode = BMO_WAITFORRENDEZVOUS;

                            // regular data acquisition
                            CompletePathData();
                            IntegratePath();
                            SavePathAndTraj();

                                UpdateAllPositions();
                                SmoothAllPositions();
                                ReparametrizeAllPositions();
                                CheckBoundaries();

                            UpdateAllPositionsFinalize();  // call STMStep++;
                            PrintSTMStepInfo();

                            if( STMStatus == ESTMS_PATH_FOUND ){
                                if( ProdPeriod <= 0 ){
                                    for(int b=0; b < NumOfBeads; b++){
                                        Beads[b]->Mode = BMO_ACCUMULATION;
                                    }
                                    STMStatus = ESTMS_COMPLETED;
                                    vout << ">> INFO: The server is terminated since all data were acquired.**" <<  endl;
                                }
                            } else {
                                if( MaxSTMSteps <= STMStep ){
                                    vout << ">> Max number of optimization steps reached, but requested convergence not reached. Server is terminating." << endl;
                                    STMStatus = ESTMS_MAX_STEPS_REACHED;
                                }
                            }
                            break;
                        case BMO_PRODUCTION:
                            CompletePathData();
                            IntegratePath();
                            STMStep++;
                            SavePathAndTraj();

                            vout << endl;

                            PrintSTMStepInfo();

                            for(int b=0; b < NumOfBeads; b++){
                                Beads[b]->Mode = BMO_PRODUCTION;
                            }
                            // this can happen only when STM with production period is run
                            // e.g. init, equi, accu periods are zero
                            STMStatus = ESTMS_COMPLETED;

                            vout << endl;
                            vout << ">> INFO: The server is terminated since all data were acquired.***" <<  endl;
                            break;
                    }
                    if( Terminate ){
                        if( STMStatus == ESTMS_COMPLETED ){
                            vout << ">> INFO: The server is terminated since it was requested by stm-admin." <<  endl;
                            STMStatus = ESTMS_COMPLETED;
                        }
                    }
                    ProcessingMutex.Unlock();
                } catch(...){
                    ProcessingMutex.Unlock();
                    RendezvousCond.BroadcastSignal();
                    NumOfRendezvousBeads = 0;
                    throw;
                }
            }
            // unblock all other beads
            RendezvousCond.BroadcastSignal();
            NumOfRendezvousBeads = 0;
        }
        RendezvousMutex.Unlock();
    } catch(...){
        RendezvousMutex.Unlock();
        throw;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::ProcessPathAsynchronously(void)
{
    try {
        ProcessingMutex.Lock();
        if( STMStatus == ESTMS_OPTIMIZING ){

            CompletePathData();
            IntegratePath();
            SavePathAndTraj();

            UpdateAllPositions();
            SmoothAllPositions();
            ReparametrizeAllPositions();
            CheckBoundaries();

            UpdateAllPositionsFinalize();  // call STMStep++;
            PrintSTMStepInfo();

            if( STMStatus == ESTMS_PATH_FOUND ){
                if( ProdPeriod <= 0 ){
                    for(int b=0; b < NumOfBeads; b++){
                        Beads[b]->Mode = BMO_ACCUMULATION;
                    }
                    STMStatus = ESTMS_COMPLETED;
                    vout << ">> INFO: The server is terminated since all data were acquired." <<  endl;
                }
            } else {
                if( MaxSTMSteps <= STMStep ){
                    vout << ">> Max number of optimization steps reached, but requested convergence not reached. Server is terminating." << endl;
                    STMStatus = ESTMS_MAX_STEPS_REACHED;
                }
            }

            if( (STMStatus == ESTMS_OPTIMIZING) || (STMStatus == ESTMS_PATH_FOUND) ){
                // move to next mode
                for(int i=0; i < NumOfBeads; i++) {
                    Beads[i]->MoveToNextMode();
                }
            }

        } else if( STMStatus == ESTMS_PATH_FOUND ){
            
            CompletePathData();
            IntegratePath();
            SavePathAndTraj();
            STMStep++;
            PrintSTMStepInfo();
            
            STMStatus = ESTMS_COMPLETED;
            vout << ">> INFO: The server is terminated since all data were acquired." <<  endl;
        } else {
            RUNTIME_ERROR("should not happen");
        }

        if( Terminate ){
            STMStatus = ESTMS_COMPLETED;
        }

        ProcessingMutex.Unlock();
    } catch(...){
        ProcessingMutex.Unlock();
        throw;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::ExchangeDataAsynchronously(CXMLElement* p_cele,CXMLElement* p_rele)
{
    if( p_cele == NULL ){
        INVALID_ARGUMENT("p_cele is NULL");
    }
    if( p_rele == NULL ){
        INVALID_ARGUMENT("p_rele is NULL");
    }

    int bead_id =  -1;
    if(p_cele->GetAttribute("bead_id",bead_id) == false) {
        RUNTIME_ERROR("unable to get bead_id");
    }

    // get bead ----------------------------------
    CBeadPtr p_bead = GetBead(bead_id);
    if( p_bead == NULL ){
        CSmallString error;
        error << "unable to find bead ID=" << bead_id;
        RUNTIME_ERROR(error);
    }

    // update bead data --------------------------
    int mode = BMO_UNKNOWN;
    p_cele->GetAttribute("mode",mode);

    if( mode != BMO_UNKNOWN ){
        // do we have correct bead mode?
        if( mode != p_bead->GetMode() ){
            CSmallString error;
            error << "client is in inconsistent mode";
            RUNTIME_ERROR(error);
        }
    }

    if( p_bead->GetModeStatus() == BMS_PREPARED ){
        // first step or bead was released
        p_bead->SetNextStepData(p_rele);
        return;
    }

    if( p_bead->GetModeStatus() != BMS_RUNNING  ){
        CSmallString error;
        error << "bead ID=" << bead_id << " is not in running state - unable to exchange data";
        RUNTIME_ERROR(error);
    }

    if( (p_bead->GetMode() == BMO_ACCUMULATION) || (p_bead->GetMode() == BMO_PRODUCTION) ){
        p_bead->GetProductionData(p_cele);
        p_bead->SetWaitForRendezvous();
        // shift to next mode is processed in Launcher
        // terminate client
        TerminateClient(p_rele);
        return;
    } else {
        p_bead->SkipProductionData();
    }

    if( (p_bead->GetMode() == BMO_INITIALIZATION) && (SoloInitPeriod == true) ){
        // terminate client
        TerminateClient(p_rele);
        return;
    }

    if( (p_bead->GetMode() == BMO_EQUILIBRATION) && (SoloEquiPeriod == true) ){
        // terminate client
        TerminateClient(p_rele);
        return;
    }

    // move to the next step

    // update program ----------------------------
    p_bead->MoveToNextMode();

    // set data for client -----------------------
    p_bead->SetNextStepData(p_rele);
}

//------------------------------------------------------------------------------

void CSTMPath::TerminateClient(CXMLElement* p_rele)
{
    if( p_rele == NULL ){
        INVALID_ARGUMENT("p_rele is NULL");
    }
    p_rele->SetAttribute("mode",BMO_TERMINATE);
    p_rele->SetAttribute("steps",0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMPath::LoadInfo(CXMLElement* p_ele)
{
    if( p_ele == NULL ){
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CXMLElement* p_mele = p_ele->GetFirstChildElement("BEADS");
    if( p_mele == NULL ){
        LOGIC_ERROR("no BEADS element");
    }

    bool result = true;
    result &= p_mele->GetAttribute("name",PathName);
    result &= p_mele->GetAttribute("ncvs",NumOfCVs);
    result &= p_mele->GetAttribute("nbeads",NumOfBeads);

    if( result == false ){
        LOGIC_ERROR("ncvs and/or nbeads is missing");
    }

    AllocatePath();

    CXMLElement* p_iele;
    p_iele = p_mele->GetFirstChildElement("COORD");
    for(int i=0; i < NumOfCVs; i++) {
        CVs[i]->LoadInfo(p_iele);
        p_iele = p_iele->GetNextSiblingElement("COORD");
    }

    CVSplines.clear();

    p_iele = p_mele->GetFirstChildElement("CVSPLINE");
    for(int i=0; i < NumOfCVs; i++) {
        if( p_iele != NULL ) {
            CSmallString cvstype = "none";
            p_iele->GetAttribute("type",cvstype);
            CCVSplinePtr cvspline;
            if( cvstype ==  "interpolating-cubic" ){
                cvspline = CCVSplinePtr(new CCVSplineInterpolatingCubic);
            } else if( cvstype ==  "smoothing-cubic" ){
                cvspline = CCVSplinePtr(new CCVSplineSmoothingCubic);
            } else {
                RUNTIME_ERROR("not implemented CV spline");
            }
            cvspline->LoadInfo(p_iele);
            CVSplines.push_back(cvspline);
        }
        p_iele = p_iele->GetNextSiblingElement("CVSPLINE");
    }

    if( CVSplines.size() != (size_t)NumOfCVs ){
        CVSplines.clear();
        for(int i=0; i < NumOfCVs; i++) {
            CCVSplinePtr cvspline;
            cvspline = CCVSplinePtr(new CCVSplineInterpolatingCubic);
            CVSplines.push_back(cvspline);
        }
    }

    SPos.CreateVector(NumOfCVs);

    p_iele = p_mele->GetFirstChildElement("BEAD");
    for(int b=0; b < NumOfBeads; b++) {
        Beads[b]->LoadInfo(p_iele);
        p_iele = p_iele->GetNextSiblingElement("BEAD");
    }
}

//------------------------------------------------------------------------------

void CSTMPath::SaveInfo(CXMLElement* p_ele)
{
    if( p_ele == NULL ){
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CXMLElement* p_mele = p_ele->CreateChildElement("BEADS");

    p_mele->SetAttribute("name",PathName);
    p_mele->SetAttribute("ncvs",NumOfCVs);
    p_mele->SetAttribute("nbeads",NumOfBeads);

    for(int i=0; i < NumOfCVs; i++) {
        CXMLElement* p_iele = p_mele->CreateChildElement("COORD");
        CVs[i]->SaveInfo(p_iele);
    }

    for(int i=0; i < NumOfCVs; i++) {
        CXMLElement* p_iele = p_mele->CreateChildElement("CVSPLINE");
        CVSplines[i]->SaveInfo(p_iele);
    }

    for(int b=0; b < NumOfBeads; b++) {
        CXMLElement* p_iele = p_mele->CreateChildElement("BEAD");
        Beads[b]->SaveInfo(p_iele);
    }
}

//------------------------------------------------------------------------------

bool CSTMPath::CheckCoords(CXMLElement* p_ele)
{
    if( p_ele == NULL ){
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CXMLIterator I(p_ele);
    int cnumofcvs = I.GetNumberOfChildElements("COORD");

    if( cnumofcvs != NumOfCVs ){
        CSmallString error;
        error << "inconsistent number of CVs declared by server (";
        error << NumOfCVs << ") and client (" << cnumofcvs << ")";
        ES_ERROR(error);
        return(false);
    }

    CXMLElement* p_nele = p_ele->GetFirstChildElement("COORD");
    int id = 0;
    while( p_nele != NULL ) {
        if( CVs[id]->CheckInfo(p_nele) == false ){
            CSmallString error;
            error << "CV" << id << " does not match server setup";
            ES_ERROR(error);
            return(false);
        }
        id++;
        p_nele = p_nele->GetNextSiblingElement("COORD");
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMPath::PrintPathSummary(std::ostream& vout)
{
    PrintPathSummaryHeader(vout);

    // optimize path and alphas for current position
    CompletePathData();
    IntegratePath();

    PrintPathSummaryData(vout);
}

//------------------------------------------------------------------------------

void CSTMPath::PrintPathSummaryHeader(std::ostream& vout)
{
    vout << "# === [PATH] ===================================================================" << endl;
    vout << "# Path name       = " << PathName << endl;
    vout << "# Number of CVs   = " << NumOfCVs << endl;
    vout << "# Number of beads = " << NumOfBeads << endl;

// header --------------------
    // legends
    vout << "#   ID   Type MO ST KinkA      α        dA/dα         Aint     CID Updates";
    for(int i=0; i < NumOfCVs; i++){
        vout << "         CV" << right << setw(2) << setfill('0') << i+1;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << "        sCV" << right << setw(2) << setfill('0') << i+1;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << "    dA/dsCV" << right << setw(2) << setfill('0') << i+1;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << "    dsCV" << right << setw(2) << setfill('0') << i+1 << "/dα";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << "     -|F" << right << setw(2) << setfill('0') << i+1 << "/dα";
    }
    vout << setfill(' ');
    vout << endl;

    // delimiters
    vout << "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

// data ----------------------
    vout << left << "#      names                                                              " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    vout << endl;
    vout << left << "#      types                                                              " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetType();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetType();
    }
    vout << endl; 
    vout << left << "#      min                                                                " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMinValue();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << 0.0;
    }
    vout << endl;
    vout << left << "#      max                                                                " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMaxValue();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << 1.0;
    }
    vout << endl;
    vout << left << "#      maxmov                                                             " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        if( CVs[i]->GetMaxMovement() > 0 ){
            vout << " " << setw(12) << CVs[i]->GetMaxMovement() * (CVs[i]->GetMaxValue() - CVs[i]->GetMinValue());
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    for(int i=0; i < NumOfCVs; i++){
        if( CVs[i]->GetMaxMovement() > 0 ){
            vout << " " << setw(12) << CVs[i]->GetMaxMovement();
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    vout << endl;

    vout << "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

    vout << "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " uuuuuuuuuuuu";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ssssssssssss";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ssssssssssss";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ssssssssssss";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ssssssssssss";
    }
    vout << endl;

    vout << "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

    vout << "#    1      2  3  4     5      6            7            8       9      10";
    int id = 11;
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    vout << endl;
    vout << "# ---- ------ -- -- ----- ------ ------------ ------------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

}

//------------------------------------------------------------------------------

void CSTMPath::PrintPathSummaryData(std::ostream& vout)
{
    for(int b=0; b < NumOfBeads; b++){
        vout << right;
        switch(Beads[b]->BeadType ){
            case(BTY_FREE):
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " F     ";
            break;
            case(BTY_PERMANENT):
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " P     ";
            break;
            case(BTY_NORMAL):
            default:
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " N     ";
            break;
        }

        switch(Beads[b]->GetMode()){
            case BMO_INITIALIZATION:
                vout << " I ";
                break;
            case BMO_ACCUMULATION:
                vout << " A ";
                break;
            case BMO_EQUILIBRATION:
                vout << " E ";
                break;
            case BMO_PRODUCTION:
                vout << " P ";
                break;
            case BMO_WAITFORRENDEZVOUS:
                vout << " W ";
                break;
            default:
                vout << " UN";
                break;
        }

        switch(Beads[b]->GetModeStatus()){
            case BMS_PREPARED:
                vout << " P ";
                break;
            case BMS_RUNNING:
                vout << " R ";
                break;
            case BMS_FINISHED:
                vout << " F ";
                break;
            default:
                vout << " UN";
                break;
        }

        vout << fixed << setprecision(1);
        vout << " " << setw(5) << Beads[b]->KinkA;

        vout << fixed << setprecision(4);
        vout << " " << setw(6) << Beads[b]->Alpha;

        vout << scientific << setprecision(5);
        vout << " " << setw(12) << Beads[b]->dAdAlpha;
        vout << " " << setw(12) << Beads[b]->A;
        if( Beads[b]->GetClientID() > 0 ){
            vout << " " << setw(7) << Beads[b]->GetClientID();
        } else {
            vout << " " << setw(7) << "--";
        }
        vout << setw(8) << Beads[b]->NumOfUpdates;
        vout << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            double unscaled = CVs[i]->GetUnscaledValue(Beads[b]->Pos[i]);
            vout << " " << setw(12) << unscaled;
        }
        for(int i=0; i < NumOfCVs; i++){
            double scaled = Beads[b]->Pos[i];
            vout << " " << setw(12) << scaled;
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b]->MF[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b]->dCVdAlpha[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b]->pMF[i];
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::PrintPathUpdate(std::ostream& vout)
{
    int num_of_updates = 0;
    for(int b=0; b < NumOfBeads; b++){
        num_of_updates += Beads[b]->NumOfUpdates;
    }

    if( num_of_updates > 0 ){

    } else {
        for(int b=0; b < NumOfBeads; b++){
            Beads[b]->Alpha = 0.0;
        }
    }

    vout << "# === [PATH UPDATE] ============================================================" << endl;
    vout << "# Path name       = " << PathName << endl;
    vout << "# Number of CVs   = " << NumOfCVs << endl;
    vout << "# Number of beads = " << NumOfBeads << endl;

// header --------------------
    // legends
    vout << "#  ID   Type  MO ST Nalpha     CID Updates  ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " old CV" << left << setw(2) << i+1 << "    ";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " new CV" << left << setw(2) << i+1 << "    ";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " diff  " << left << setw(2) << i+1 << "    ";
    }
    vout << endl;

    // delimiters
    vout << "# ---- ------ -- -- ------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

// data ----------------------
    vout << left << "#      names                           " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    vout << endl;
    vout << left << "#      types                           " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetType();
    }
    vout << endl;
    vout << left << "#      min                             " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMinValue();
    }
    vout << endl;
    vout << left << "#      max                             " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMaxValue();
    }
    vout << endl;
    vout << left << "#      maxmov                          " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        if( CVs[i]->GetMaxMovement() > 0 ){                     // unscalled
            vout << " " << setw(12) << CVs[i]->GetMaxMovement() * (CVs[i]->GetMaxValue() - CVs[i]->GetMinValue());
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    vout << endl;

    vout << "# ---- ------ -- -- ------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

    vout << "#    1      2  3  4      5       6       7";
    int id = 8;
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << right << setw(13) << id;
        id++;
    }
    vout << endl;
    vout << "# ---- ------ -- -- ------ ------- -------";
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " ------------";
    }
    vout << endl;

    for(int b=0; b < NumOfBeads; b++){
        vout << right;

        switch(Beads[b]->BeadType ){
            case(BTY_FREE):
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " F     ";
            break;
            case(BTY_PERMANENT):
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " P     ";
            break;
            case(BTY_NORMAL):
            default:
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " N     ";
            break;
        }
        
        switch(Beads[b]->GetMode()){
            case BMO_INITIALIZATION:
                vout << " I ";
                break;
            case BMO_ACCUMULATION:
                vout << " A ";
                break;
            case BMO_EQUILIBRATION:
                vout << " E ";
                break;
            case BMO_PRODUCTION:
                vout << " P ";
                break;
            case BMO_WAITFORRENDEZVOUS:
                vout << " W ";
                break;
            default:
                vout << " UN";
                break;
        }

        switch(Beads[b]->GetModeStatus()){
            case BMS_PREPARED:
                vout << " P ";
                break;
            case BMS_RUNNING:
                vout << " R ";
                break;
            case BMS_FINISHED:
                vout << " F ";
                break;
            default:
                vout << " UN";
                break;
        }

        vout << fixed << setprecision(4);
        vout << " " << setw(6) << Beads[b]->Alpha;

        if( Beads[b]->GetClientID() > 0 ){
            vout << " " << setw(7) << Beads[b]->GetClientID();
        } else {
            vout << " " << setw(7) << " --";
        }
        vout << setw(8) << Beads[b]->NumOfUpdates;
        vout << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            double unscaled = CVs[i]->GetUnscaledValue(Beads[b]->OPos[i]);
            vout << " " << setw(12) << unscaled;
        }
        if( Beads[b]->NumOfUpdates > 0 ){
            for(int i=0; i < NumOfCVs; i++){
                double unscaled = CVs[i]->GetUnscaledValue(Beads[b]->Pos[i]);
                vout << " " << setw(12) << unscaled;
            }
            for(int i=0; i < NumOfCVs; i++){
                double unscaled1 = CVs[i]->GetUnscaledValue(Beads[b]->Pos[i]);
                double unscaled2 = CVs[i]->GetUnscaledValue(Beads[b]->OPos[i]);
                vout << " " << setw(12) << unscaled1 - unscaled2;
            }
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::PrintPath(std::ostream& vout)
{
    vout << "[PATH]" << endl;
    vout << "name     " << PathName << endl;
    vout << "ncvs     " << NumOfCVs << endl;
    vout << "nbeads   " << NumOfBeads << endl;
    vout << "names    ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetName();
    }
    vout << endl;
    vout << "types    ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetType();
    }
    vout << endl;
    vout << "min      ";
    vout << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMinValue();
    }
    vout << endl;
    vout << "max      ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMaxValue();
    }
    vout << endl;
    vout << "maxmov   ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMaxMovement() * (CVs[i]->GetMaxValue() - CVs[i]->GetMinValue());;
    }
    vout << endl;

    for(int b=0; b < NumOfBeads; b++){
        vout << right;
        switch(Beads[b]->BeadType ){
            case(BTY_FREE):
            vout << "free     ";
            break;
            case(BTY_PERMANENT):
            vout << "permanent";
            break;
            case(BTY_NORMAL):
            default:
            vout << "normal   ";
            break;
        }

        for(int i=0; i < NumOfCVs; i++){
            double scaled = Beads[b]->Pos[i];
            double unscaled = CVs[i]->GetUnscaledValue(scaled);
            vout << " " << setw(12) << unscaled;
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::SplitString(string text,vector<string>& words)
{
    size_t i=0;
    char ch;
    string word;

    while( i < text.size() ) {
        ch = text[i++];
        if ( isspace(ch) )  {
          if (!word.empty()) {
            words.push_back(word);
          }
          word = "";
        }
        else {
          word += ch;
        }
    }
    if (!word.empty())  {
        words.push_back(word);
    }
}

//------------------------------------------------------------------------------

int CSTMPath::GetNumOfBeads(void)
{
    return(NumOfBeads);
}

//------------------------------------------------------------------------------

int CSTMPath::GetNumOfBeadsInRendezvousState(void)
{
    int count = 0;
    ProcessingMutex.Lock();

    for(int i=0; i < NumOfBeads; i++){
        if( Beads[i]->GetMode() == BMO_WAITFORRENDEZVOUS ) count++;
    }

    ProcessingMutex.Unlock();
    return(count);
}

//------------------------------------------------------------------------------

ESTMState CSTMPath::GetSTMStatus(void)
{
    return(STMStatus);
}

//------------------------------------------------------------------------------

bool CSTMPath::IsAsynchronous(void)
{
    return(AsynchronousMode);
}

//------------------------------------------------------------------------------

CBeadPtr CSTMPath::GetBead(int bead_id)
{
    if( (bead_id <= 0) || (bead_id > NumOfBeads) ){
        LOGIC_ERROR("bead_id out-of-legal range");
    }
    return(Beads[bead_id-1]);
}

//------------------------------------------------------------------------------

CBeadPtr CSTMPath::GetBeadByClientID(int client_id)
{
    for(int b=0; b < NumOfBeads; b++){
        if( Beads[b]->GetClientID() == client_id ){
            return(Beads[b]);
        }
    }
    return(NULL);
}

//------------------------------------------------------------------------------

CBeadPtr CSTMPath::GetNextFreeBead(void)
{
    for(int id=0; id < NumOfBeads; id++){
        if( Beads[id]->GetClientID() == -1 ){
            Beads[id]->SetClientID(0);
            return(Beads[id]);
        }
    }
    return(CBeadPtr());
}

//------------------------------------------------------------------------------

int CSTMPath::GetSTMStep(void)
{
    int step = 0;

    ProcessingMutex.Lock();
        step = STMStep;
    ProcessingMutex.Unlock();

    return(step);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMPath::PrintSTMHeader(void)
{
    PrintSTMHeaderF(vout);
    PrintSTMHeaderF(OptLogFOut);
    HeaderPrinted = true;

    MABufPLenChange.CreateVector(MABufLength);    
    MABufMaxBeadMove.CreateVector(MABufLength); 
    MABufAveBeadMove.CreateVector(MABufLength);   
    MABufMaxpMFSize.CreateVector(MABufLength);    
    MABufAvepMFSize.CreateVector(MABufLength); 
}

//------------------------------------------------------------------------------

void CSTMPath::PrintSTMStepInfo(void)
{
    CalculateSTMStepStat();
    PrintSTMStepInfoF(vout);
    PrintSTMStepInfoF(OptLogFOut);
}

//------------------------------------------------------------------------------

void CSTMPath::PrintSTMHeaderF(std::ostream& fout)
{
    fout << "#" << endl;
    fout << "# NOTE: All values are in scalled units." << endl;
    fout << "#" << endl;
    fout << "#     |                                         Current Values                                          |    |                             Running Averages                             |" << endl;
    fout << "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----|--------------|-------------- --------------|-------------- --------------|" << endl;
    fout << "# Step|  Path length  Length change | Max bead move BID Ave bead move | Max pMF size  BID  Ave pMF size |Term|Length change | Max bead move Ave bead move | Max pMF size   Ave pMF size |" << endl;
    fout << "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----|--------------|-------------- --------------|-------------- --------------|" << endl;
    fout << "#    1|             2              3|             4   5              6|             7   8              9|  10|            11|            12             13|            14             15|" << endl;
    fout << "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----|--------------|-------------- --------------|-------------- --------------|" << endl;
}

//------------------------------------------------------------------------------

void CSTMPath::CalculateSTMStepStat(void)
{

    PLenChange = UpdatedPathLength-CurrentPathLength;

    MaxBeadMove = 0;
    MaxBeadMoveID = 0;
    AveBeadMove = 0;

    MaxpMFSize = 0;
    MaxpMFSizeID = 0;
    AvepMFSize = 0;

    int bn = 0;
    for(int b=0; b < NumOfBeads; b++){
        if( Beads[b]->BeadType == BTY_PERMANENT ) continue; // skipt permanent beads
        double bmov = 0;
        double mfsize = 0.0;
        for(int i=0; i < NumOfCVs; i++){
            bmov += (Beads[b]->FPos[i]-Beads[b]->OPos[i])*(Beads[b]->FPos[i]-Beads[b]->OPos[i]);
            mfsize += (Beads[b]->uMF[i])*(Beads[b]->uMF[i]);
        }
        bmov = sqrt(bmov);
        mfsize = sqrt(mfsize);

        AveBeadMove += bmov;
        if( bmov > MaxBeadMove ){
            MaxBeadMove = bmov;
            MaxBeadMoveID = b+1;
        }

        AvepMFSize += mfsize;
        if( mfsize > MaxpMFSize ){
            MaxpMFSize = mfsize;
            MaxpMFSizeID = b+1;
        }
        bn++;
    }
    if( bn > 0 ){
        AveBeadMove = AveBeadMove / (double)bn;
        AvepMFSize = AvepMFSize / (double)bn;
    }

    // shft the moving average buffers
    for(int i=0; i < MABufLength-1; i++){
        MABufPLenChange[i]  = MABufPLenChange[i+1];
        MABufMaxBeadMove[i] = MABufMaxBeadMove[i+1];
        MABufAveBeadMove[i] = MABufAveBeadMove[i+1];
        MABufMaxpMFSize[i]  = MABufMaxpMFSize[i+1];
        MABufAvepMFSize[i]  = MABufAvepMFSize[i+1];
    }

    // add new values
    MABufPLenChange[MABufLength-1]  = fabs(PLenChange);
    MABufMaxBeadMove[MABufLength-1] = MaxBeadMove;
    MABufAveBeadMove[MABufLength-1] = AveBeadMove;
    MABufMaxpMFSize[MABufLength-1]  = MaxpMFSize;
    MABufAvepMFSize[MABufLength-1]  = AvepMFSize;

    // get moving averages
    MAPLenChange    = 0.0;    
    MAMaxBeadMove   = 0.0;
    MAAveBeadMove   = 0.0;  
    MAMaxpMFSize    = 0.0;    
    MAAvepMFSize    = 0.0;
    TermCrit        = 0;

    if( STMStep < MABufLength ) return;

    for(int i=0; i < MABufLength; i++){
        MAPLenChange    += MABufPLenChange[i];    
        MAMaxBeadMove   += MABufMaxBeadMove[i];
        MAAveBeadMove   += MABufAveBeadMove[i];  
        MAMaxpMFSize    += MABufMaxpMFSize[i];    
        MAAvepMFSize    += MABufAvepMFSize[i];
    }

    if( MABufLength > 0 ) {
        MAPLenChange    /= (double)MABufLength;    
        MAMaxBeadMove   /= (double)MABufLength;
        MAAveBeadMove   /= (double)MABufLength; 
        MAMaxpMFSize    /= (double)MABufLength;    
        MAAvepMFSize    /= (double)MABufLength;
    }

    // determine number of fullfilled termination criteria
    if( MAPLenChange < FinalPLenChange ){
        TermCrit++;
    }
    if( MAMaxBeadMove < FinalMaxBeadMove ){
        TermCrit++;
    }
    if( MAAveBeadMove < FinalAveBeadMove ){
        TermCrit++;
    }
    if( MAMaxpMFSize < FinalMaxpMFSize ){
        TermCrit++;
    }
    if( MAAvepMFSize < FinalAvepMFSize ){
        TermCrit++;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::PrintSTMStepInfoF(std::ostream& fout)
{
    fout << left << setw(6)  << STMStep << " " << right;

    fout << setw(14) << setprecision(7) << scientific << CurrentPathLength << " ";
    fout << setw(14) << setprecision(7) << scientific << PLenChange << " ";    

    fout << setw(14) << setprecision(7) << scientific << MaxBeadMove << " ";
    fout << setw(3) << MaxBeadMoveID << " ";
    fout << setw(14) << setprecision(7) << scientific << AveBeadMove << " ";

    fout << setw(14) << setprecision(7) << scientific << MaxpMFSize << " ";
    fout << setw(3) << MaxpMFSizeID << " ";
    fout << setw(14) << setprecision(7) << scientific << AvepMFSize << " ";

    fout << " " << setw(1) << TermCrit << "/" << "5 ";

    fout << setw(14) << setprecision(7) << scientific << MAPLenChange << " ";  

    fout << setw(14) << setprecision(7) << scientific << MAMaxBeadMove << " ";
    fout << setw(14) << setprecision(7) << scientific << MAAveBeadMove << " ";

    fout << setw(14) << setprecision(7) << scientific << MAMaxpMFSize << " ";
    fout << setw(14) << setprecision(7) << scientific << MAAvepMFSize << " ";

    fout << endl;

    if( (TermCrit == 5) && (STMStatus != ESTMS_PATH_FOUND) ){
        STMStatus = ESTMS_PATH_FOUND;
        fout << "#" << endl;
        fout << "# >> INFO: The path have converged." << endl;
        if( ProdPeriod > 0 ){
            fout << "# >> INFO: Entering production accumulation (" << ProdPeriod <<" steps)." <<  endl;
        }
    }
}

//------------------------------------------------------------------------------

void CSTMPath::SavePathAndTraj(void)
{
    // write output and trajectory
    if( (OutInterval > 0) && (STMStep % OutInterval == 0) ){
        SavePath(OutputPath);
        SavePathSummary(OutputPathSummary);
    }
    if( (TrajInterval > 0) && (STMStep % TrajInterval == 0) ){
        SaveTrajectorySnapshot();
    }
}

//------------------------------------------------------------------------------

void CSTMPath::CompletePathData(void)
{
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->ResetPosUpdates();
    }

// needed by OptimizePath()
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->PPos = Beads[b]->Pos;
        Beads[b]->SegLength = 0.0;
    }

// update segments
    InitPathSegments();

// per path segments
    std::vector< std::vector<CBeadPtr> >::iterator it = PathSegments.begin();
    std::vector< std::vector<CBeadPtr> >::iterator ie = PathSegments.end();

    while( it != ie ){
        std::vector<CBeadPtr>& bl = *it;
        double seglen = OptimizePath(bl);

        for(size_t b=0; b < bl.size(); b++){
            bl[b]->CalcBead();
            bl[b]->SegLength = seglen;
        }
        it++;
    }

    CurrentPathLength = OptimizePath(Beads);

    CSimpleVector<double> v1,v2;
    v1.CreateVector(NumOfCVs);
    v2.CreateVector(NumOfCVs); 

// calculate kink angle
    for(int b=0; b < NumOfBeads; b++){
        if( (b == 0) || (b == NumOfBeads - 1) ){
            Beads[b]->KinkA = 180.0;
            continue;
        }

        for(int i=0; i < NumOfCVs; i++){
            v1[i] = Beads[b+1]->Pos[i] - Beads[b]->Pos[i];
            v2[i] = Beads[b-1]->Pos[i] - Beads[b]->Pos[i];
        }

        double dot_product = 0.0;
        double norm2_v1 = 0.0;
        double norm2_v2 = 0.0;
        for(int i=0; i < NumOfCVs; i++){
            dot_product += v1[i] * v2[i];
            norm2_v1 += v1[i] * v1[i];
            norm2_v2 += v2[i] * v2[i];
        }

        if( (norm2_v1 > 0.0) && (norm2_v2 > 0.0)  ){
            double cos_theta = dot_product / (sqrt(norm2_v1) * sqrt(norm2_v2));
            if( cos_theta < -1.0 ) cos_theta = -1.0;
            if( cos_theta >  1.0 ) cos_theta =  1.0;

            Beads[b]->KinkA = acos(cos_theta) * 180.0 / M_PI;
        } else {
            Beads[b]->KinkA = 0.0;
        }
    }
}

//------------------------------------------------------------------------------

void CSTMPath::UpdateAllPositions(void)
{
    STMStep++;

    if( OptMethod == "gd" ){
        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->UpdatePositionGD(StepSize);
        }
    } else if ( OptMethod == "ngd" ){
        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->UpdatePositionNGD(StepSize,MinGNormEps);
        }
    } else if ( OptMethod == "ngd-auto" ){
        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->UpdatePositionNGDAuto(StepSize,MaxGNormForGD,MinGNormEps);
        }
    } else if ( OptMethod == "adam" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionADAM(StepSize,AdamB1,AdamB2,MinGNormEps);
        }
        if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ){
            vout << ">> INFO: Reset ADAM memory." << endl;
            ResetAdamAlg--;
        }
    } else if ( OptMethod == "adabelief" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionADABelief(StepSize,AdamB1,AdamB2,MinGNormEps);
        }
        if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ){
            vout << ">> INFO: Reset ADABelif memory." << endl;
            ResetAdamAlg--;
        }
    } else if ( OptMethod == "amsgrad" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionAMSGrad(StepSize,AdamB1,AdamB2,MinGNormEps);
        }
        if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ){
            vout << ">> INFO: Reset AMSGrad memory." << endl;
            ResetAdamAlg--;
        }
    } else if ( OptMethod == "amsgradbc" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionAMSGradBC(StepSize,AdamB1,AdamB2,MinGNormEps);
        }
        if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ){
            vout << ">> INFO: Reset AMSGradBC memory." << endl;
            ResetAdamAlg--;
        }
    } else {
        vout << endl;
        vout << ">> INFO: The optimization method '" << OptMethod << "' is not implemented (CSTMPath::UpdateAllPositions)!" << endl;
        STMStatus = ESTMS_MAX_STEPS_REACHED;
        return;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::UpdateAllPositionsFinalize(void)
{
    for(int b=0; b < NumOfBeads; b++) {
        Beads[b]->UpdatePositionFinalize();
    }
}

//------------------------------------------------------------------------------

void CSTMPath::SmoothAllPositions(void)
{
    if( (SmoothInterval == 0) || (STMStep % SmoothInterval != 0) ){
        for(int b=0; b < NumOfBeads; b++) {
            Beads[b]->SPos = Beads[b]->NPos;
        }
        return;
    }

    // vout << debug << "Smoothing positions ..." << endl << high;

    // smooth path
    for(int i=0; i < NumOfBeads; i++){
        if( (i != 0) && (i != NumOfBeads-1) && (Beads[i]->BeadType == BTY_NORMAL) ){
            for(int j=0; j < NumOfCVs; j++){
                Beads[i]->SPos[j] = (1.0-SmoothingFac)*Beads[i]->NPos[j]
                                + 0.5*SmoothingFac*(Beads[i-1]->NPos[j]+Beads[i+1]->NPos[j]);
            }
        } else {
            for(int j=0; j < NumOfCVs; j++){
                Beads[i]->SPos[j] = Beads[i]->NPos[j];
            }
        }
    }
}

//------------------------------------------------------------------------------

void CSTMPath::ReparametrizeAllPositions(void)
{
    if( (ReparamInterval == 0) || (STMStep % ReparamInterval != 0) ){
        for(int b=0; b < NumOfBeads; b++) {
            Beads[b]->FPos = Beads[b]->SPos;
        }
        return;
    }

   // vout << debug << "Re-parametrizing positions ..." << endl << high;

    // re-optimize path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->PPos = Beads[b]->SPos;
        Beads[b]->FPos = Beads[b]->SPos;
    }

    // update segments
    InitPathSegments();

    // per path segments
    std::vector< std::vector<CBeadPtr> >::iterator it = PathSegments.begin();
    std::vector< std::vector<CBeadPtr> >::iterator ie = PathSegments.end();

    while( it != ie ){
        std::vector<CBeadPtr>& bl = *it;

        if( bl.size() < 3 ) continue;

        OptimizePath(bl);

        for(size_t b=1; b < bl.size()-1; b++){  // skip terminals

            double alpha = (double)b / ((double)bl.size()-1.0);
            for(int i=0; i < NumOfCVs; i++){
                if( Beads[b]->BeadType == BTY_NORMAL ){
                    Beads[b]->FPos[i] = CVSplines[i]->GetCV(alpha);
                }
            }
        }

        it++;
    }
}

//------------------------------------------------------------------------------

void CSTMPath::InitPathSegments(void)
{
    PathSegments.clear();

    int nbeads = Beads.size();
    int i = 0;

    while( i < nbeads ){
        if( (i != 0) && (Beads[i]->BeadType != BTY_FREE) ){
            CSmallString error;
            error << "path segment must start with the free bead or the first bead of the path, bidx: " << i+1;
            RUNTIME_ERROR(error)
        }
        int seg_first = i;

        i++;

        // Find the end of the segment
        while( (i < nbeads) && (Beads[i]->BeadType != BTY_FREE) ) {
            i++;
        }
        int seg_last = i;

        // make a list
        std::vector<CBeadPtr> segment;
        for(int idx = seg_first; idx < seg_last; idx++){
            if( (idx < 0) || (idx >= nbeads) ) continue;
            segment.push_back(Beads[idx]);
        }
        PathSegments.push_back(segment);
    }

    // cout << "PSS: " << PathSegments.size() <<  endl; 
}

//------------------------------------------------------------------------------

void CSTMPath::CheckBoundaries(void)
{
    // now we are in scalled coordinates
    for(int b=0; b < NumOfBeads; b++){
        for(int i=0; i < NumOfCVs; i++){
            if( Beads[b]->FPos[i] < 0.0 ){
                Beads[b]->FPos[i] = 0.0;
            }
            if( Beads[b]->FPos[i] > 1.0 ){
                Beads[b]->FPos[i] = 1.0;
            }
        }
    }

    // get data about the final path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->PPos  = Beads[b]->FPos;
    }
    UpdatedPathLength = OptimizePath(Beads);
}

//------------------------------------------------------------------------------

void CSTMPath::IntegratePath(void)
{
    // get PMF projections along path
    double fes = 0.0;
    for(int b=0; b < NumOfBeads; b++){
        double a = 0.0;
        // get bead derivative along path
        for(int i=0; i < NumOfCVs; i++) {
            Beads[b]->dCVdAlpha[i] = CVSplines[i]->GetCVFirstDer(Beads[b]->Alpha);
            double sc = CVs[i]->GetMaxValue() - CVs[i]->GetMinValue();
            // dAdAlphas are per path segment, thus correction * CurrentPathLength / bead->SegLength
            a += Beads[b]->dCVdAlpha[i] * sc * Beads[b]->MF[i] * CurrentPathLength / Beads[b]->SegLength;
        }
        Beads[b]->dAdAlpha = a;
        if( b > 0 ){
            fes += 0.5*(Beads[b]->Alpha - Beads[b-1]->Alpha)*(Beads[b]->dAdAlpha + Beads[b-1]->dAdAlpha);
        }
        Beads[b]->A = fes;
    }

    if( ShifGlobalMinA2Zero ) {
        // get value of global minima
        double min = Beads[0]->A;
        for(int b=1; b < NumOfBeads; b++){
            if( min > Beads[b]->A ){
                min = Beads[b]->A;
            } 
        }
        for(int b=0; b < NumOfBeads; b++){
            Beads[b]->A = Beads[b]->A - min;
        }
    }
}

//------------------------------------------------------------------------------

void CSTMPath::SetServerTerminated(void)
{
    STMStatus = ESTMS_COMPLETED;
    RendezvousCond.BroadcastSignal();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CSTMPath::OptimizePath(std::vector<CBeadPtr>& beads)
{
    if( beads.size() < 2 ){
        RUNTIME_ERROR("beads.size() must be greater or equal 2");
    }

    // get initial path length from linear interpolation
    double tot_length = 0;
    for(size_t b=1; b < beads.size(); b++){
        double slen = 0;
        for(int i=0; i < NumOfCVs; i++){
            slen += (beads[b]->PPos[i]-beads[b-1]->PPos[i])*(beads[b]->PPos[i]-beads[b-1]->PPos[i]);
        }
        tot_length += sqrt(slen);
    }

    if( tot_length == 0 ){
        RUNTIME_ERROR("path has zero length");
    }

    // get initial alphas from linear interpolation
    beads[0]->Alpha = 0.0;
    double path_length = 0;
    for(size_t b=1; b < beads.size()-1; b++){
        double slen = 0;
        for(int i=0; i < NumOfCVs; i++){
            slen += (beads[b]->PPos[i]-beads[b-1]->PPos[i])*(beads[b]->PPos[i]-beads[b-1]->PPos[i]);
        }
        if( slen == 0 ){
            RUNTIME_ERROR("path segment has zero length");
        }
        path_length += sqrt(slen);
        beads[b]->Alpha = path_length/tot_length;
    }
    beads[beads.size()-1]->Alpha = 1.0;

    // interpolate CVS
    for(int i=0; i < NumOfCVs; i++){
        CVSplines[i]->Allocate(beads.size());
        for(size_t b=0; b < beads.size(); b++){
            CVSplines[i]->SetPoint(b,beads[b]->Alpha,beads[b]->PPos[i]);
        }
        CVSplines[i]->BuildSpline();
    }

   // vout << debug;
   // vout << "Path length = " << tot_length << endl;

    return(tot_length);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

