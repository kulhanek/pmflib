// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
* GD        - gradient descent (gradient)
* NGD       - normalized gradient descent (normalized gradient)
* NGD-AUTO  - gradient descent (switch between GD and NGD)
* ADAM      - Adaptive Moment Estimation
* AMSGrad   - AMSGrad
* AMSGradBC - AMSGrad + bias corrected estimates
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
    OutputPathSummary = "_stm.results";
    PathTrajectory = "_stm.traj";

    // intervals
    TrajInterval    = 0;
    OutInterval     = 100;
    SmoothInterval  = 0;
    ReparamInterval = 1;

    // stm
    InitPeriod = 10000;         // initialization period
    EquiPeriod =  1000;         // equilibration period
    AccuPeriod =  5000;         // accumulation period
    ProdPeriod = 50000;         // final production period

    MaxSTMSteps         = 100;
    OptMethod           = "amsgradbc";
    FinalMaxPLenChange  = 0.005;
    FinalMaxMovement    = 0.01;
    FinalAveMovement    = 0.01;
    FinalpMFSizeMax     = 4.00;
    FinalpMFSizeAve     = 1.00;

    MaxGNormForGD       = 5.0;
    MinGNormEps         = 1e-7;
    StepSize            = 0.1;
    AdamB1              = 0.9;
    AdamB2              = 0.999;
    ResetAdamAlg        = 0;
    MemoryLength        = 0;

    SmoothingFac        = 0.0;

    AsynchronousMode = false;    // update per bead or path

    STMStep         = 0;
    MaxMovement     = 0.0;          // current max path movement
    MaxMovementBead = 0;        // current max path movement is for given bead
    AveMovement     = 0;            // current average path movement
    pMFSizeMax      = 0;
    pMFSizeAve      = 0;

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
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary
             << "  (default)" << endl;
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory
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

    if(prmfile.GetStringByKey("summary",OutputPathSummary) == true) {
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary << endl;
    } else {
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary
             << "  (default)" << endl;
    }

    if(prmfile.GetStringByKey("trajectory",PathTrajectory) == true) {
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory << endl;
    } else {
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory
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
        vout << "Max number of STM steps (steps)                = " << setw(9) << MaxSTMSteps
             << left << "             (default)" << endl;
        vout << "Optimization method (optmethod)                = " << setw(9) << OptMethod
             << left << "             (default)" << endl;

        vout << "Max final path length change (maxfplch)        = " << setw(9) << FinalMaxPLenChange
             << left << "             (default)" << endl;
        vout << "Max final path movement (maxfmove)             = " << setw(9) << FinalMaxMovement
             << left << "             (default)" << endl;
        vout << "Average final path movement (avefmove)         = " << setw(9) << FinalAveMovement
             << left << "             (default)" << endl;

        vout << "Max perpendicular mean force (maxfpmf)         = " << setw(9) << FinalpMFSizeMax
             << left << "             (default)" << endl;
        vout << "Average perpendicular mean force (avefpmf)     = " << setw(9) << FinalpMFSizeAve
             << left << "             (default)" << endl;

        vout << "Initialization period (init)                   = " << setw(9) << InitPeriod
             << left << "             (default)" << endl;
        vout << "Accumulation period (accu)                     = " << setw(9) << AccuPeriod
             << left << "             (default)" << endl;
        vout << "Equilibration period (equi)                    = " << setw(9) << EquiPeriod
             << left << "             (default)" << endl;
        vout << "Final production period (prod)                 = " << setw(9) << ProdPeriod
             << left << "             (default)" << endl;
        vout << "Asynchronous mode (async)                      = " << setw(9) << right << PrmFileOnOff(AsynchronousMode)
             << left << "             (default)" << endl;
        vout << "Path smoothing factor (sfac)                   = " << setw(9) << SmoothingFac
             << left << "             (default)" << endl;
        vout << "Asynchronous mode (async)                      = " << setw(9) << right << PrmFileOnOff(AsynchronousMode)
             << left << "             (default)" << endl;

        // opt method - GD
        bool result  = ProcessGDOptMethodSetup(prmfile);
        return(result);
    }

    if(prmfile.GetIntegerByKey("steps",MaxSTMSteps) == true) {
        vout << "Max number of STM steps (steps)                = " << setw(9) << MaxSTMSteps << left << endl;
    } else {
        vout << "Max number of STM steps (steps)                = " << setw(9) << MaxSTMSteps
             << left << "             (default)" << endl;
    }

    if(prmfile.GetStringByKey("optmethod",OptMethod) == true) {
        vout << "Optimization method (optmethod)                = " << setw(9) << OptMethod << left << endl;
    } else {
        vout << "Optimization method (optmethod)                = " << setw(9) << OptMethod
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("maxfplch",FinalMaxPLenChange) == true) {
        vout << "Max final path length change (maxfplch)        = " << setw(9) << FinalMaxPLenChange << left << endl;
    } else {
        vout << "Max final path length change (maxfplch)        = " << setw(9) << FinalMaxPLenChange
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("maxfmove",FinalMaxMovement) == true) {
        vout << "Max final path movement (maxfmove)             = " << setw(9) << FinalMaxMovement << left << endl;
    } else {
        vout << "Max final path movement (maxfmove)             = " << setw(9) << FinalMaxMovement
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("avefmove",FinalAveMovement) == true) {
        vout << "Average final path movement (avefmove)         = " << setw(9) << FinalAveMovement << left << endl;
    } else {
        vout << "Average final path movement (avefmove)         = " << setw(9) << FinalAveMovement
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("avefpmf",FinalpMFSizeAve) == true) {
        vout << "Average perpendicular mean force (avefpmf)     = " << setw(9) << FinalpMFSizeAve << left << endl;
    } else {
        vout << "Average perpendicular mean force (avefpmf)     = " << setw(9) << FinalpMFSizeAve
             << left << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("maxfpmf",FinalpMFSizeMax) == true) {
        vout << "Max perpendicular mean force (maxfpmf)         = " << setw(9) << FinalpMFSizeMax << left << endl;
    } else {
        vout << "Max perpendicular mean force (maxfpmf)         = " << setw(9) << FinalpMFSizeMax
             << left << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("init",InitPeriod) == true) {
        vout << "Initialization period (init)                   = " << setw(9) << InitPeriod << endl;
    } else {
        vout << "Initialization period (init)                   = " << setw(9) << InitPeriod
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("accu",AccuPeriod) == true) {
        vout << "Accumulation period (accu)                     = " << setw(9) << AccuPeriod << endl;
    } else {
        vout << "Accumulation period (accu)                     = " << setw(9) << AccuPeriod
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("equi",EquiPeriod) == true) {
        vout << "Equilibration period (equi)                    = " << setw(9) << EquiPeriod << endl;
    } else {
        vout << "Equilibration period (equi)                    = " << setw(9) << EquiPeriod
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("prod",ProdPeriod) == true) {
        vout << "Final production period (prod)                 = " << setw(9) << ProdPeriod << endl;
    } else {
        vout << "Final production period (prod)                 = " << setw(9) << ProdPeriod
             << "             (default)" << endl;
    }

    if(prmfile.GetDoubleByKey("sfac",SmoothingFac) == true) {
        vout << "Path smoothing factor (sfac)                   = " << setw(9) << SmoothingFac << left << endl;
    } else {
        vout << "Path smoothing factor (sfac)                   = " << setw(9) << SmoothingFac
             << left << "             (default)" << endl;
    }

    if(prmfile.GetLogicalByKey("async",AsynchronousMode) == true) {
        vout << "Asynchronous mode (async)                      = " << setw(9) << right << PrmFileOnOff(AsynchronousMode) << left << endl;
    } else {
        vout << "Asynchronous mode (async)                      = " << setw(9) << right << PrmFileOnOff(AsynchronousMode)
             << left << "             (default)" << endl;
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
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval
             << "             (default)" << endl;
        vout << "Output path update (output)                    = " << setw(9) << OutInterval
             << "             (default)" << endl;
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval
             << "             (default)" << endl;
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval
             << "             (default)" << endl;
        return(true);
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

    if(prmfile.GetIntegerByKey("smooth",SmoothInterval) == true) {
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval << endl;
    } else {
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval
             << "             (default)" << endl;
    }

    if(prmfile.GetIntegerByKey("reparam",ReparamInterval) == true) {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval << endl;
    } else {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval
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
//    ! flexible   cv1 cv2 cv3 ... cvn
//    ! permanent  cv1 cv2 cv3 ... cvn

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

    // optimize path
    for(int i=0; i < num_of_user_beads; i++){
        InputBeads[i]->PPos = InputBeads[i]->Pos;
    }
    OptimizePath(InputBeads);

    // generate missing points or re-optimize path
    Beads[0]->Alpha = 0.0;
    Beads[0]->Permanent = InputBeads[0]->Permanent;
    Beads[0]->BeadID = 1;
    Beads[NumOfBeads-1]->Alpha = 1.0;
    Beads[NumOfBeads-1]->Permanent = InputBeads[num_of_user_beads-1]->Permanent;
    Beads[NumOfBeads-1]->BeadID = NumOfBeads;
    for(int i=0; i < NumOfCVs; i++){
        Beads[0]->Pos[i] = CVSplines[i]->GetCV(0.0);
        Beads[NumOfBeads-1]->Pos[i] = CVSplines[i]->GetCV(1.0);
        for(int b=1; b < NumOfBeads-1; b++){
            double alpha = (double)b / ((double)NumOfBeads-1.0);
            Beads[b]->Pos[i] = CVSplines[i]->GetCV(alpha);
            Beads[b]->Alpha = alpha;
            Beads[b]->BeadID = b + 1;
        }
    }

    // check boundaries
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->FPos = Beads[b]->Pos;
    }
    CheckBoundaries();

    // and again re-optimize path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->PPos = Beads[b]->FPos;
        Beads[b]->Pos  = Beads[b]->FPos;
    }
    OptimizePath(Beads);

    // and final correct positions
    Beads[0]->Alpha = 0.0;
    Beads[NumOfBeads-1]->Alpha = 1.0;
    for(int i=0; i < NumOfCVs; i++){
        Beads[0]->Pos[i] = CVSplines[i]->GetCV(0.0);
        Beads[NumOfBeads-1]->Pos[i] = CVSplines[i]->GetCV(1.0);
        for(int b=1; b < NumOfBeads-1; b++){
            double alpha = (double)b / ((double)NumOfBeads-1.0);
            Beads[b]->Pos[i] = CVSplines[i]->GetCV(alpha);
            Beads[b]->Alpha = alpha;
        }
    }

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
            double max = CSmallString(tokens[i]).ToDouble();
            CVs[i]->SetMaxMovement(max);
            vout << " " << setw(12) << max;
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

        if( (tokens[0] != "flexible") && (tokens[0] != "permanent") ){
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

        if( (tokens[0] != "flexible") && (tokens[0] != "permanent") ){
            CSmallString error;
            error << "unsupported key '" << tokens[0] << "'";
            RUNTIME_ERROR(error)
        }
        if( beadid > NumOfBeads ){
            RUNTIME_ERROR("more beads specification than nbeads");
        }

        // process permanent or flexible point definition
        beads[beadid]->InitBead(this,NumOfCVs);
        for(int i=0; i < NumOfCVs; i++){
            beads[beadid]->Pos[i] = CSmallString(tokens[i+1]).ToDouble();
        }
        beads[beadid]->Permanent = tokens[0] != "flexible";

        if( tokens[0] == "flexible" ) {
            vout << setw(4) << beadid+1 << " F     " << scientific << setprecision(5);
        } else {
            vout << setw(4) << beadid+1 << " P     " << scientific << setprecision(5);
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << beads[beadid]->Pos[i];
        }
        vout << endl;
        beadid++;
    }
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

    Trajectory.open(PathTrajectory);
    if( ! Trajectory ){
        ES_ERROR("unable to open path trajectory file");
        return(false);
    }

    // write header
    Trajectory << "# STMTRAJ " << NumOfCVs << " " << NumOfBeads << endl;
    PrintPathSummaryHeader(Trajectory);

    return(true);
}

//------------------------------------------------------------------------------

void CSTMPath::SaveSnapshot(void)
{
    if( TrajInterval <= 0 ) return;
    Trajectory << "# STMSNAP " << STMStep / TrajInterval << endl;
    PrintPathSummaryData(Trajectory);
    Trajectory << endl; // necessary for gnuplot
}

//------------------------------------------------------------------------------

void CSTMPath::CloseTrajectory(void)
{
    Trajectory.close();
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
                            UsedStepSize = StepSize;
                            for(int i=0; i < 5; i++){
                                UpdateAllPositions();
                                SmoothAllPositions();
                                ReparametrizeAllPositions();
                                CheckBoundaries();
                                if( ! isnan(UpdatedPathLength) ) break;
                                vout << ">> WARNING: Stability problem - reducing step size!" <<  endl;
                                UsedStepSize = UsedStepSize / 2.0;
                            }
                            UpdateAllPositionsFinalize();
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
            UsedStepSize = StepSize;
            for(int i=0; i < 5; i++){
                UpdateAllPositions();
                SmoothAllPositions();
                ReparametrizeAllPositions();
                CheckBoundaries();
                if( ! isnan(UpdatedPathLength) ) break;
                vout << ">> WARNING: Stability problem - reducing step size!" <<  endl;
                UsedStepSize = UsedStepSize / 2.0;
            }
            UpdateAllPositionsFinalize();
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
        p_bead->WaitForRendezvous();
        // shift to next mode is processed in Launcher
    } else {
        p_bead->SkipProductionData();
    }

// FIXME
// should it be here?
//    if( (STMStatus == ESTMS_MAX_STEPS_REACHED) || (STMStatus == ESTMS_COMPLETED) ) return;

    // terminate client
    TerminateClient(p_rele);
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
    vout << "#  ID   Type  ST  alpha  dA/dalpha       A            CID Updates";
    for(int i=0; i < NumOfCVs; i++){
        vout << "     CV" << left << setw(2) << i+1 << "    ";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << "   dA/dCV" << left << setw(2) << i+1 << "  ";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " dCV" << left << setw(2) << i+1 << "/dalpha";
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " -|F" << left << setw(2) << i+1 << "/dalpha";
    }
    vout << endl;

    // delimiters
    vout << "# ---- ------ -- ------ ------------ ------------ ------- -------";
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
    vout << left << "#      names                                                    " << right;
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
    vout << left << "#      types                                                    " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetType();
    }
    vout << endl;
    vout << left << "#      min                                                      " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMinValue();
    }
    vout << endl;
    vout << left << "#      max                                                      " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i]->GetMaxValue();
    }
    vout << endl;
    vout << left << "#      maxmov                                                   " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        if( CVs[i]->GetMaxMovement() > 0 ){
            vout << " " << setw(12) << CVs[i]->GetMaxMovement();
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    vout << endl;

    vout << "# ---- ------ -- ------ ------------ ------------ ------- -------";
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

    vout << "#    1      2  3      4            5            6       7       8";
    int id = 9;
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
    vout << "# ---- ------ -- ------ ------------ ------------ ------- -------";
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
        if( Beads[b]->Permanent ) {
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " P     ";
        } else {
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " F     ";
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
            vout << " " << setw(12) << Beads[b]->Pos[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b]->MF[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b]->dCV[i];
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
    vout << "#  ID   Type  ST Nalpha     CID Updates  ";
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
    vout << "# ---- ------ -- ------ ------- -------";
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
        if( CVs[i]->GetMaxMovement() > 0 ){
            vout << " " << setw(12) << CVs[i]->GetMaxMovement();
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    vout << endl;

    vout << "# ---- ------ -- ------ ------- -------";
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

    vout << "#    1      2  3      4       5       6";
    int id = 7;
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
    vout << "# ---- ------ -- ------ ------- -------";
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
        if( Beads[b]->Permanent ) {
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " P     ";
        } else {
            vout << "  " << setw(4) << Beads[b]->GetBeadID() << setw(7) << " F     ";
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
            vout << " " << setw(12) << Beads[b]->OPos[i];
        }
        if( Beads[b]->NumOfUpdates > 0 ){
            for(int i=0; i < NumOfCVs; i++){
                vout << " " << setw(12) << Beads[b]->Pos[i];
            }
            for(int i=0; i < NumOfCVs; i++){
                vout << " " << setw(12) << Beads[b]->Pos[i] - Beads[b]->OPos[i];
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
        vout << " " << setw(12) << CVs[i]->GetMaxMovement();
    }
    vout << endl;

    for(int b=0; b < NumOfBeads; b++){
        vout << right;
        if( Beads[b]->Permanent ) {
            vout << "permanent";
        } else {
            vout << "flexible ";
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b]->Pos[i];
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
    vout << endl;
    vout << "# Step|  Path length  Length change | Max movement  BID  Ave movement | Max pMF size  BID  Ave pMF size |Term" << endl;
    vout << "# ----|-------------- --------------|-------------- --- --------------|-------------- --- --------------|----" << endl;

    HeaderPrinted = true;
}

//------------------------------------------------------------------------------

void CSTMPath::PrintSTMStepInfo(void)
{
    int termcrit = 0;

    vout << left << setw(6)  << STMStep << " " << right;
    vout << setw(14) << setprecision(7) << scientific << CurrentPathLength << " ";

    if( fabs(UpdatedPathLength-CurrentPathLength) < FinalMaxPLenChange ){
        vout << "<green>" << setw(14) << setprecision(7) << scientific << UpdatedPathLength-CurrentPathLength << "</green> ";
        termcrit++;
    } else {
        vout << setw(14) << setprecision(7) << scientific << UpdatedPathLength-CurrentPathLength << " ";
    }

    MaxMovement = 0;
    AveMovement = 0;
    pMFSizeAve = 0;
    pMFSizeMax = 0;

    for(int b=0; b < NumOfBeads; b++){
        double mov = 0;
        double mfsize = 0.0;
        for(int i=0; i < NumOfCVs; i++){
            mov += (Beads[b]->FPos[i]-Beads[b]->OPos[i])*(Beads[b]->FPos[i]-Beads[b]->OPos[i]);
            mfsize += (Beads[b]->pMF[i])*(Beads[b]->pMF[i]);
            //mfsize += (Beads[b]->mkold[i])*(Beads[b]->mkold[i]);
        }

        mov = sqrt(mov);
        AveMovement += mov; // add mov square
        if( mov > MaxMovement ){
            MaxMovement = mov;
            MaxMovementBead = b+1;
        }

        mfsize = sqrt(mfsize);
        pMFSizeAve += mfsize;
        if( mfsize > pMFSizeMax ){
            pMFSizeMax = mfsize;
            MaxpMFBead = b+1;
        }
    }
    AveMovement = AveMovement / (double)NumOfBeads;
    pMFSizeAve = pMFSizeAve / (double)NumOfBeads;

    if( MaxMovement < FinalMaxMovement ){
        vout << "<green>" << setw(14) << setprecision(7) << scientific <<  MaxMovement << "</green> ";
        termcrit++;
    } else {
        vout << setw(14) << setprecision(7) << scientific << MaxMovement << " ";
    }

    vout << setw(3) << MaxMovementBead << " ";

    if( AveMovement < FinalAveMovement ){
        vout << "<green>" << setw(14) << setprecision(7) << scientific << AveMovement << "</green> ";
        termcrit++;
    } else {
        vout << setw(14) << setprecision(7) << scientific << AveMovement << " ";
    }

    if( pMFSizeMax < FinalpMFSizeMax  ){
        vout << "<green>" << setw(14) << setprecision(7) << scientific <<  pMFSizeMax << "</green> ";
        termcrit++;
    } else {
        vout << setw(14) << setprecision(7) << scientific << pMFSizeMax << " ";
    }

    vout << setw(3) << MaxpMFBead << " ";

    if( pMFSizeAve < FinalpMFSizeAve ){
        vout << "<green>" << setw(14) << setprecision(7) << scientific <<  pMFSizeAve << "</green> ";
        termcrit++;
    } else {
        vout << setw(14) << setprecision(7) << scientific << pMFSizeAve << " ";
    }

    vout << " " << setw(1) << termcrit << "/" << "5";
    vout << endl;

    if( (termcrit == 5) && (STMStatus != ESTMS_PATH_FOUND) ){
        STMStatus = ESTMS_PATH_FOUND;
        vout << endl;
        vout << ">> INFO: The path have converged." << endl;
        if( ProdPeriod > 0 ){
            vout << ">> INFO: Entering production accumulation (" << ProdPeriod <<" steps)." <<  endl;
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
        SaveSnapshot();
    }
}

//------------------------------------------------------------------------------

void CSTMPath::CompletePathData(void)
{
    for(int i=0; i < NumOfBeads; i++){
        Beads[i]->ResetPosUpdates();
    }

    // re-optimize path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b]->PPos = Beads[b]->Pos;
    }
    CurrentPathLength = OptimizePath(Beads);

    for(int i=0; i < NumOfBeads; i++){
        Beads[i]->CalcProjector();
    }
}

//------------------------------------------------------------------------------

void CSTMPath::UpdateAllPositions(void)
{
    STMStep++;

    if( OptMethod == "gd" ){
        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->UpdatePositionGD(UsedStepSize);
        }
    } else if ( OptMethod == "ngd" ){
        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->UpdatePositionNGD(UsedStepSize,MinGNormEps);
        }
    } else if ( OptMethod == "ngd-auto" ){
        for(int i=0; i < NumOfBeads; i++){
            Beads[i]->UpdatePositionNGDAuto(UsedStepSize,MaxGNormForGD,MinGNormEps);
        }
    } else if ( OptMethod == "adam" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionADAM(UsedStepSize,AdamB1,AdamB2,MinGNormEps);
        }
        if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ){
            vout << ">> INFO: Reset ADAM memory." << endl;
            ResetAdamAlg--;
        }
    } else if ( OptMethod == "amsgrad" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionAMSGrad(UsedStepSize,AdamB1,AdamB2,MinGNormEps);
        }
        if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ){
            vout << ">> INFO: Reset AMSGrad memory." << endl;
            ResetAdamAlg--;
        }
    } else if ( OptMethod == "amsgradbc" ){
        for(int i=0; i < NumOfBeads; i++){
            if( (ResetAdamAlg > 0) && (MemoryLength > 0) && (STMStep % MemoryLength == 0) ) Beads[i]->ResetADAM();
            Beads[i]->UpdatePositionAMSGradBC(UsedStepSize,AdamB1,AdamB2,MinGNormEps);
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
        if( (i == 0) || (i == NumOfBeads-1) || (Beads[i]->Permanent) ){
            for(int j=0; j < NumOfCVs; j++){
                Beads[i]->SPos[j] = Beads[i]->NPos[j];
            }
        } else {
            for(int j=0; j < NumOfCVs; j++){
                Beads[i]->SPos[j] = (1.0-SmoothingFac)*Beads[i]->NPos[j]
                                + 0.5*SmoothingFac*(Beads[i-1]->NPos[j]+Beads[i+1]->NPos[j]);
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
    }
    OptimizePath(Beads);

    // and correct positions
    for(int i=0; i < NumOfCVs; i++){
        if( Beads[0]->Permanent ){
            Beads[0]->FPos[i] = Beads[0]->SPos[i];
        } else {
            Beads[0]->FPos[i] = CVSplines[i]->GetCV(0.0);
        }
        if( Beads[NumOfBeads-1]->Permanent ){
            Beads[NumOfBeads-1]->FPos[i] = Beads[NumOfBeads-1]->SPos[i];
        } else {
            Beads[NumOfBeads-1]->FPos[i] = CVSplines[i]->GetCV(1.0);
        }

        for(int b=1; b < NumOfBeads-1; b++){
            double alpha = (double)b / ((double)NumOfBeads-1.0);
            if( Beads[b]->Permanent ){
                Beads[b]->FPos[i] = Beads[b]->SPos[i];
            } else {
                Beads[b]->FPos[i] = CVSplines[i]->GetCV(alpha);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CSTMPath::CheckBoundaries(void)
{
    for(int b=0; b < NumOfBeads; b++){
        for(int i=0; i < NumOfCVs; i++){
            if( Beads[b]->FPos[i] < CVs[i]->GetMinValue() ){
                Beads[b]->FPos[i] = CVs[i]->GetMinValue();
            }
            if( Beads[b]->FPos[i] > CVs[i]->GetMaxValue() ){
                Beads[b]->FPos[i] = CVs[i]->GetMaxValue();
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
            Beads[b]->dCV[i] = CVSplines[i]->GetCVFirstDer(Beads[b]->Alpha);
            a += Beads[b]->dCV[i]*Beads[b]->MF[i];
        }
        Beads[b]->dAdAlpha = a;
        if( b > 0 ){
            fes += 0.5*(Beads[b]->Alpha - Beads[b-1]->Alpha)*(Beads[b]->dAdAlpha + Beads[b-1]->dAdAlpha);
        }
        Beads[b]->A = fes;
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

   // vout << debug;
   // vout << "Initial path length = " << tot_length << endl;

    double prev_length = 0;

    for(int s=0; s < 1000; s++){
        // interpolate CVS
        for(int i=0; i < NumOfCVs; i++){
            CVSplines[i]->Allocate(beads.size());
            for(size_t b=0; b < beads.size(); b++){
                CVSplines[i]->SetPoint(b,beads[b]->Alpha,beads[b]->PPos[i]);
            }
            CVSplines[i]->BuildSpline();
        }

        prev_length = tot_length;

        // determine new path length
        tot_length = 0;
        for(size_t b=1; b < beads.size(); b++){
            tot_length += GetSegmentLength(beads[b-1]->Alpha,beads[b]->Alpha);
        }

        // vout << "Optimized path length = " << tot_length << endl;

        if( fabs(tot_length-prev_length) < 1e-7 ){
        //     vout << "Converged path length = " << tot_length << endl;
            return(tot_length);
        }

        // determine new alphas
        beads[0]->Alpha = 0.0;
        double path_length = 0;
        double prev_alpha = beads[0]->Alpha;
        for(size_t b=1; b < beads.size()-1; b++){
            path_length += GetSegmentLength(prev_alpha,beads[b]->Alpha);
            beads[b]->Alpha = path_length/tot_length;
            prev_alpha = beads[b]->Alpha;
        }
        beads[beads.size()-1]->Alpha = 1.0;
    }

    return(tot_length);
}

//------------------------------------------------------------------------------

double CSTMPath::GetSegmentLength(double alpha1,double alpha2)
{

   // cout << "cv-splines: " << CVSplines.size() << endl;
   // cout << "SPOS: " << SPos.GetLength() << endl;

    double len = 0;
    for(int i=0; i < NumOfCVs; i++){
        SPos[i] = CVSplines[i]->GetCV(alpha1);
    }
    double step = (alpha2-alpha1)/SegmentDiscretization;
    double alpha = alpha1 + step;
    while( alpha < alpha2 ){
        double slen2 = 0;
        for(int i=0; i < NumOfCVs; i++){
            double curr = CVSplines[i]->GetCV(alpha);
            slen2 +=  (curr-SPos[i])*(curr-SPos[i]);
            SPos[i] = curr;
        }
        len += sqrt(slen2);
        alpha += step;
    }

    double slen2 = 0;
    for(int i=0; i < NumOfCVs; i++){
        double last = CVSplines[i]->GetCV(alpha2);
        slen2 +=  (last-SPos[i])*(last-SPos[i]);
    }
    len += sqrt(slen2);

    return(len);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

