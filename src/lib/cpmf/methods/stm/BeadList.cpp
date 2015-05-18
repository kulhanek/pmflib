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
#include "BeadList.hpp"
#include <TemplIterator.hpp>
#include <ErrorSystem.hpp>
#include <iomanip>
#include <string>
#include <iterator>
#include <math.h>
#include <XMLIterator.hpp>
#include <PrmUtils.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CBeadList::CBeadList(void)
{
    // files
    InputPath = "{PATHS}";
    OutputPath = "_stm.path";
    OutputPathSummary = "_stm.results";
    PathTrajectory = "_stm.traj";

    // intervals
    TrajInterval = 0;
    OutInterval = 100;
    SmoothInterval = 1;
    ReparamInterval = 1;

    // stm
    InitPeriod = 10000;         // initialization period
    AccuPeriod =  5000;         // accumulation period
    EquiPeriod =  1000;         // equilibration period
    ProdPeriod = 50000;         // final production period

    MaxSTMSteps         = 250;
    StepSize            = 0.001;
    FinalMaxPLenChange  = 0.01;
    FinalMaxMovement    = 0.01;
    FinalAveMovement    = 0.01;

    SmoothingFac = 0.1;
    AsynchronousMode = false;    // update per bead or path

    STMStep = 0;
    MaxMovement = 0.0;          // current max path movement
    MaxMovementBead = 0;        // current max path movement is for given bead
    AveMovement = 0;            // current average path movement

    // control
    NumOfRendezvousBeads = 0;
    STMStatus = ESTMS_INITIALIZED;
    HeaderPrinted = false;
    Terminate = false;

    // how many points are used to calculate path segment length
    SegmentDiscretization = 10;

    ClearPath();
}

//------------------------------------------------------------------------------

void CBeadList::AllocatePath(void)
{
    if( NumOfCVs < 2 ){
        RUNTIME_ERROR("number of CVs must be larger than or equal to 2");
    }
    if( NumOfBeads < 3 ){
        RUNTIME_ERROR("number of beads must be larger than or equal to 3");
    }

    CVs.CreateVector(NumOfCVs);
    CVSplines.CreateVector(NumOfCVs);
    Beads.CreateVector(NumOfBeads);
    SPos.CreateVector(NumOfCVs);

    for(int b=0; b < NumOfBeads; b++){
        Beads[b].InitBead(this,NumOfCVs);
    }
}

//------------------------------------------------------------------------------

void CBeadList::ClearPath(void)
{
    NumOfBeads = 0;
    NumOfCVs = 0;
    CVs.FreeVector();
    Beads.FreeVector();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CBeadList::AttachVerboseStream(std::ostream& str,bool verbose)
{
    vout.Attach(str);
    if( verbose ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }
}

//------------------------------------------------------------------------------

void CBeadList::ProcessFilesControl(CPrmFile& file)
{
    vout << endl;
    vout << "=== [files] ====================================================================" << endl;
    if(file.OpenSection("files") == false) {
        vout << "Input path (input)                             = " << setw(20) << InputPath
             << "  (default)" << endl;
        vout << "Output path (output)                           = " << setw(20) << OutputPath
             << "  (default)" << endl;
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary
             << "  (default)" << endl;
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory
             << "  (default)" << endl;
        return;
    }

    if(file.GetStringByKey("input",InputPath) == true) {
        vout << "Input path (input)                             = " << setw(20) << InputPath << endl;
    } else {
        vout << "Input path (input)                             = " << setw(20) << InputPath
             << "  (default)" << endl;
    }

    if(file.GetStringByKey("output",OutputPath) == true) {
        vout << "Output path (output)                           = " << setw(20) << OutputPath << endl;
    } else {
        vout << "Output path (output)                           = " << setw(20) << OutputPath
             << "  (default)" << endl;
    }

    if(file.GetStringByKey("summary",OutputPathSummary) == true) {
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary << endl;
    } else {
        vout << "Output path summary (summary)                  = " << setw(20) << OutputPathSummary
             << "  (default)" << endl;
    }

    if(file.GetStringByKey("trajectory",PathTrajectory) == true) {
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory << endl;
    } else {
        vout << "Path trajectory (trajectory)                   = " << setw(20) << PathTrajectory
             << "  (default)" << endl;
    }
}

//------------------------------------------------------------------------------

void CBeadList::ProcessSTMControl(CPrmFile& file)
{
    vout << endl;
    vout << "=== [stm] ======================================================================" << endl;
    if(file.OpenSection("stm") == false) {
        vout << "Max number of STM steps (steps)                = " << setw(9) << MaxSTMSteps
             << left << "             (default)" << endl;
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
        vout << "Max final path length change (maxfplch)        = " << setw(9) << FinalMaxPLenChange
             << left << "             (default)" << endl;
        vout << "Max final path movement (maxfmove)             = " << setw(9) << FinalMaxMovement
             << left << "             (default)" << endl;
        vout << "Average final path movement (avefmove)         = " << setw(9) << FinalAveMovement
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
        return;
    }

    if(file.GetIntegerByKey("steps",MaxSTMSteps) == true) {
        vout << "Max number of STM steps (steps)                = " << setw(9) << MaxSTMSteps << left << endl;
    } else {
        vout << "Max number of STM steps (steps)                = " << setw(9) << MaxSTMSteps
             << left << "             (default)" << endl;
    }

    if(file.GetDoubleByKey("stepsize",StepSize) == true) {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize << left << endl;
    } else {
        vout << "Step size (stepsize)                           = " << setw(9) << StepSize
             << left << "             (default)" << endl;
    }

    if(file.GetDoubleByKey("maxfplch",FinalMaxPLenChange) == true) {
        vout << "Max final path length change (maxfplch)        = " << setw(9) << FinalMaxPLenChange << left << endl;
    } else {
        vout << "Max final path length change (maxfplch)        = " << setw(9) << FinalMaxPLenChange
             << left << "             (default)" << endl;
    }

    if(file.GetDoubleByKey("maxfmove",FinalMaxMovement) == true) {
        vout << "Max final path movement (maxfmove)             = " << setw(9) << FinalMaxMovement << left << endl;
    } else {
        vout << "Max final path movement (maxfmove)             = " << setw(9) << FinalMaxMovement
             << left << "             (default)" << endl;
    }

    if(file.GetDoubleByKey("avefmove",FinalAveMovement) == true) {
        vout << "Average final path movement (avefmove)         = " << setw(9) << FinalAveMovement << left << endl;
    } else {
        vout << "Average final path movement (avefmove)         = " << setw(9) << FinalAveMovement
             << left << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("init",InitPeriod) == true) {
        vout << "Initialization period (init)                   = " << setw(9) << InitPeriod << endl;
    } else {
        vout << "Initialization period (init)                   = " << setw(9) << InitPeriod
             << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("accu",AccuPeriod) == true) {
        vout << "Accumulation period (accu)                     = " << setw(9) << AccuPeriod << endl;
    } else {
        vout << "Accumulation period (accu)                     = " << setw(9) << AccuPeriod
             << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("equi",EquiPeriod) == true) {
        vout << "Equilibration period (equi)                    = " << setw(9) << EquiPeriod << endl;
    } else {
        vout << "Equilibration period (equi)                    = " << setw(9) << EquiPeriod
             << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("prod",ProdPeriod) == true) {
        vout << "Final production period (prod)                 = " << setw(9) << ProdPeriod << endl;
    } else {
        vout << "Final production period (prod)                 = " << setw(9) << ProdPeriod
             << "             (default)" << endl;
    }

    if(file.GetDoubleByKey("sfac",SmoothingFac) == true) {
        vout << "Path smoothing factor (sfac)                   = " << setw(9) << SmoothingFac << left << endl;
    } else {
        vout << "Path smoothing factor (sfac)                   = " << setw(9) << SmoothingFac
             << left << "             (default)" << endl;
    }

    if(file.GetLogicalByKey("async",AsynchronousMode) == true) {
        vout << "Asynchronous mode (async)                      = " << setw(9) << right << PrmFileOnOff(AsynchronousMode) << left << endl;
    } else {
        vout << "Asynchronous mode (async)                      = " << setw(9) << right << PrmFileOnOff(AsynchronousMode)
             << left << "             (default)" << endl;
    }
}

//------------------------------------------------------------------------------

void CBeadList::ProcessIntervalsControl(CPrmFile& file)
{
    vout << endl;
    vout << "=== [intervals] ================================================================" << endl;
    if(file.OpenSection("intervals") == false) {
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval
             << "             (default)" << endl;
        vout << "Output path update (output)                    = " << setw(9) << OutInterval
             << "             (default)" << endl;
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval
             << "             (default)" << endl;
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval
             << "             (default)" << endl;
        return;
    }

    if(file.GetIntegerByKey("trajectory",TrajInterval) == true) {
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval << endl;
    } else {
        vout << "Trajectory interval (trajectory)               = " << setw(9) << TrajInterval
             << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("output",OutInterval) == true) {
        vout << "Output path update (output)                    = " << setw(9) << OutInterval << endl;
    } else {
        vout << "Output path update (output)                    = " << setw(9) << OutInterval
             << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("reparam",SmoothInterval) == true) {
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval << endl;
    } else {
        vout << "Path smoothing interval (smooth)               = " << setw(9) << SmoothInterval
             << "             (default)" << endl;
    }

    if(file.GetIntegerByKey("smooth",ReparamInterval) == true) {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval << endl;
    } else {
        vout << "Path reparametrization interval (reparam)      = " << setw(9) << ReparamInterval
             << "             (default)" << endl;
    }
}

//------------------------------------------------------------------------------

void CBeadList::ProcessPathControl(CPrmFile& file)
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
    if(file.OpenSection("PATH") == false) {
        RUNTIME_ERROR("[PATH] section is mandatory for a path specification");
    }

    if(file.GetStringByKey("name",PathName) == true) {
        vout << "Path name                         = " << PathName << endl;
    } else {
        RUNTIME_ERROR("path name (name) is not specified");
    }

    // read number of beads and CVs
    if(file.GetIntegerByKey("nbeads",NumOfBeads) == true) {
        vout << "Number of beads (nbeads)          = " << NumOfBeads << endl;
    } else {
        RUNTIME_ERROR("number of beads (nbeads) is not specified");
    }
    if(file.GetIntegerByKey("ncvs",NumOfCVs) == true) {
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
    ReadPathControls(file);

    // read temporary path - only points specified by user
    int num_of_user_beads = ReadPathNumberOfUserBeads(file);

    vout << debug << "Number of user provided beads: " << num_of_user_beads << endl;
    vout << high;

    CSimpleVector<CBead>    beads;
    beads.CreateVector(num_of_user_beads);

    ReadPathUserBeads(file,beads);

    // optimize path
    for(int i=0; i < num_of_user_beads; i++){
        beads[i].PPos = beads[i].Pos;
    }
    OptimizePath(beads);

    // generate missing points or reoptimize path
    Beads[0].Alpha = 0.0;
    Beads[0].Permanent = beads[0].Permanent;
    Beads[0].BeadID = 1;
    Beads[NumOfBeads-1].Alpha = 1.0;
    Beads[NumOfBeads-1].Permanent = beads[num_of_user_beads-1].Permanent;
    Beads[NumOfBeads-1].BeadID = NumOfBeads;
    for(int i=0; i < NumOfCVs; i++){
        Beads[0].Pos[i] = CVSplines[i].GetCV(0.0);
        Beads[NumOfBeads-1].Pos[i] = CVSplines[i].GetCV(1.0);
        for(int b=1; b < NumOfBeads-1; b++){
            double alpha = (double)b / ((double)NumOfBeads-1.0);
            Beads[b].Pos[i] = CVSplines[i].GetCV(alpha);
            Beads[b].Alpha = alpha;
            Beads[b].BeadID = b + 1;
        }
    }

    // check boundaries
    for(int b=0; b < NumOfBeads; b++){
        Beads[b].RPos = Beads[b].Pos;
    }
    CheckBoundaries();

    // and again reoptimize path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b].PPos = Beads[b].RPos;
    }
    OptimizePath(Beads);

    // and final correct positions
    Beads[0].Alpha = 0.0;
    Beads[NumOfBeads-1].Alpha = 1.0;
    for(int i=0; i < NumOfCVs; i++){
        Beads[0].Pos[i] = CVSplines[i].GetCV(0.0);
        Beads[NumOfBeads-1].Pos[i] = CVSplines[i].GetCV(1.0);
        for(int b=1; b < NumOfBeads-1; b++){
            double alpha = (double)b / ((double)NumOfBeads-1.0);
            Beads[b].Pos[i] = CVSplines[i].GetCV(alpha);
            Beads[b].Alpha = alpha;
        }
    }
}

//------------------------------------------------------------------------------

void CBeadList::ReadPathControls(CPrmFile& file)
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
            CVs[i].ID = i;
            CVs[i].SetName(tokens[i]);
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
            CVs[i].SetType(tokens[i]);
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
            CVs[i].SetMinValue(min);
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
            CVs[i].SetMaxValue(max);
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
            CVs[i].SetMaxMovement(max);
            vout << " " << setw(12) << max;
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

int CBeadList::ReadPathNumberOfUserBeads(CPrmFile& file)
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

void CBeadList::ReadPathUserBeads(CPrmFile& file,CSimpleVector<CBead>& beads)
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
        beads[beadid].InitBead(this,NumOfCVs);
        for(int i=0; i < NumOfCVs; i++){
            beads[beadid].Pos[i] = CSmallString(tokens[i+1]).ToDouble();
        }
        beads[beadid].Permanent = tokens[0] != "flexible";

        if( tokens[0] == "flexible" ) {
            vout << setw(4) << beadid+1 << " F     " << scientific << setprecision(5);
        } else {
            vout << setw(4) << beadid+1 << " P     " << scientific << setprecision(5);
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << beads[beadid].Pos[i];
        }
        vout << endl;
        beadid++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CBeadList::LoadPath(CPrmFile& file)
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

bool CBeadList::SavePath(void)
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

void CBeadList::SavePath(const CSmallString& name)
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

bool CBeadList::SavePathSummary(void)
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

void CBeadList::SavePathSummary(const CSmallString& name)
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

void CBeadList::FlushPath(void)
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

bool CBeadList::OpenTrajectory(void)
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

    // write initial bead positions
    SaveSnapshot();

    return(true);
}

//------------------------------------------------------------------------------

void CBeadList::SaveSnapshot(void)
{
    if( TrajInterval <= 0 ) return;
    Trajectory << "# STMSNAP " << STMStep / TrajInterval << endl;
    PrintPathSummaryData(Trajectory);
    Trajectory << endl; // necessary for gnuplot
}

//------------------------------------------------------------------------------

void CBeadList::CloseTrajectory(void)
{
    Trajectory.close();
}

//------------------------------------------------------------------------------

void CBeadList::ForceTermination(void)
{
    ProcessingMutex.Lock();
        Terminate = true;
        vout << ">>> INFO: Received soft termination signal." << endl;
    ProcessingMutex.Unlock();
}

//------------------------------------------------------------------------------

void CBeadList::SetAsynchronousMode(bool set)
{
    AsynchronousMode = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CBeadList::CheckClient(CXMLElement* p_cele)
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

void CBeadList::RegisterBead(int bead_id,int client_id)
{
    try {
        ProcessingMutex.Lock();

        CBead* p_bead = GetBead(bead_id);
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
                    Beads[i].MoveToNextMode();
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

void CBeadList::BeginAsynchronousMode(void)
{
    // move to first mode
    for(int i=0; i < NumOfBeads; i++){
        Beads[i].MoveToNextMode();
    }
}

//------------------------------------------------------------------------------

void CBeadList::ExchangeData(CXMLElement* p_cele,CXMLElement* p_rele)
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

void CBeadList::ExchangeDataSynchronously(CXMLElement* p_cele,CXMLElement* p_rele)
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
    CBead* p_bead = GetBead(bead_id);
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

void CBeadList::ProcessProductionData(CBead* p_bead)
{
    // how many beads are waiting
    RendezvousMutex.Lock();
        NumOfRendezvousBeads++;
        if( NumOfRendezvousBeads != NumOfBeads ){
            p_bead->Mode = BMO_WAITFORRENDEZVOUS;
            RendezvousCond.WaitForSignal(RendezvousMutex);
        } else {
            // do all operations on the whole path
            try{
                if( STMStatus == ESTMS_PATH_FOUND ){
                    // this will happen if the path was found and final production runs are required
                    STMStatus = ESTMS_COMPLETED;
                    vout << ">> INFO: The server is terminated since all data were acquired." <<  endl;
                } else {
                    ProcessingMutex.Lock();
                    switch( p_bead->GetMode() ){
                        case BMO_ACCUMULATION:
                            p_bead->Mode = BMO_WAITFORRENDEZVOUS;
                            // regular data accusition
                            UpdateAllPositions();
                            SmoothAllPositions();
                            ReparametrizeAllPositions();
                            CheckBoundaries();
                            PrintSTMStepInfo();
                            if( STMStatus == ESTMS_PATH_FOUND ){
                                if( ProdPeriod <= 0 ){
                                    for(int b=0; b < NumOfBeads; b++){
                                        Beads[b].Mode = BMO_ACCUMULATION;
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
                            break;
                        case BMO_PRODUCTION:
                            for(int b=0; b < NumOfBeads; b++){
                                Beads[b].Mode = BMO_PRODUCTION;
                            }
                            // this can happen only when STM with production period is run
                            // e.g. init, equi, accu periods are zero
                            STMStatus = ESTMS_COMPLETED;
                            vout << ">> INFO: The server is terminated since all data were acquired." <<  endl;
                            break;
                    }
                    if( Terminate ){
                        if( STMStatus == ESTMS_COMPLETED ){
                            vout << ">> INFO: The server is terminated since it was requested by stm-admin." <<  endl;
                            STMStatus = ESTMS_COMPLETED;
                        }
                    }
                    ProcessingMutex.Unlock();
                }
            } catch(...) {
                RendezvousCond.BroadcastSignal();
                NumOfRendezvousBeads = 0;
                ProcessingMutex.Unlock();
                RendezvousMutex.Unlock();
                throw;
            }
            // unblock all other beads
            RendezvousCond.BroadcastSignal();
            NumOfRendezvousBeads = 0;
        }
    RendezvousMutex.Unlock();
}

//------------------------------------------------------------------------------

void CBeadList::ProcessPathAsynchronously(void)
{
    ProcessingMutex.Lock();

        if( STMStatus == ESTMS_OPTIMIZING ){
            UpdateAllPositions();
            SmoothAllPositions();
            ReparametrizeAllPositions();
            CheckBoundaries();
            PrintSTMStepInfo();

            if( STMStatus == ESTMS_PATH_FOUND ){
                if( ProdPeriod <= 0 ){
                    for(int b=0; b < NumOfBeads; b++){
                        Beads[b].Mode = BMO_ACCUMULATION;
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
                    Beads[i].MoveToNextMode();
                }
            }
        } else if( STMStatus == ESTMS_PATH_FOUND ){
            STMStatus = ESTMS_COMPLETED;
            vout << ">> INFO: The server is terminated since all data were acquired." <<  endl;
        } else {
            RUNTIME_ERROR("should not happen")
        }

        if( Terminate ){
            STMStatus = ESTMS_COMPLETED;
        }

    ProcessingMutex.Unlock();
}

//------------------------------------------------------------------------------

void CBeadList::ExchangeDataAsynchronously(CXMLElement* p_cele,CXMLElement* p_rele)
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
    CBead* p_bead = GetBead(bead_id);
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
        // first step or bead was reeased
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

void CBeadList::TerminateClient(CXMLElement* p_rele)
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

void CBeadList::LoadInfo(CXMLElement* p_ele)
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
        CVs[i].LoadInfo(p_iele);
        p_iele = p_iele->GetNextSiblingElement("COORD");
    }

    p_iele = p_mele->GetFirstChildElement("BEAD");
    for(int b=0; b < NumOfBeads; b++) {
        Beads[b].LoadInfo(p_iele);
        p_iele = p_iele->GetNextSiblingElement("BEAD");
    }
}

//------------------------------------------------------------------------------

void CBeadList::SaveInfo(CXMLElement* p_ele)
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
        CVs[i].SaveInfo(p_iele);
    }

    for(int b=0; b < NumOfBeads; b++) {
        CXMLElement* p_iele = p_mele->CreateChildElement("BEAD");
        Beads[b].SaveInfo(p_iele);
    }
}

//------------------------------------------------------------------------------

bool CBeadList::CheckCoords(CXMLElement* p_ele)
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
        if( CVs[id].CheckInfo(p_nele) == false ){
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

void CBeadList::PrintPathSummary(std::ostream& vout)
{
    PrintPathSummaryHeader(vout);
    PrintPathSummaryData(vout);
}

//------------------------------------------------------------------------------

void CBeadList::PrintPathSummaryHeader(std::ostream& vout)
{
    vout << "# === [PATH] ===================================================================" << endl;
    vout << "# Path name       = " << PathName << endl;
    vout << "# Number of CVs   = " << NumOfCVs << endl;
    vout << "# Number of beads = " << NumOfBeads << endl;

// header --------------------
    // legends
    vout << "#  ID   Type  ST  alpha  dA/dalpha       A         CID Updates";
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
    vout << "# ---- ------ -- ------ ------------ ------------ ---- -------";
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
    vout << left << "#      names                                                 " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    vout << endl;
    vout << left << "#      types                                                 " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetType();
    }
    vout << endl;
    vout << left << "#      min                                                   " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMinValue();
    }
    vout << endl;
    vout << left << "#      max                                                   " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMaxValue();
    }
    vout << endl;
    vout << left << "#      maxmov                                                " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        if( CVs[i].GetMaxMovement() > 0 ){
            vout << " " << setw(12) << CVs[i].GetMaxMovement();
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    vout << endl;

    vout << "# ---- ------ -- ------ ------------ ------------ ---- -------";
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

    vout << "#    1      2  3      4            5            6    7       8";
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
    vout << "# ---- ------ -- ------ ------------ ------------ ---- -------";
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

void CBeadList::PrintPathSummaryData(std::ostream& vout)
{
    CalculatePathData();

    for(int b=0; b < NumOfBeads; b++){
        vout << right;
        if( Beads[b].Permanent ) {
            vout << "  " << setw(4) << Beads[b].GetBeadID() << setw(7) << " P     ";
        } else {
            vout << "  " << setw(4) << Beads[b].GetBeadID() << setw(7) << " F     ";
        }
        switch(Beads[b].GetMode()){
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
        vout << " " << setw(6) << Beads[b].Alpha;

        vout << scientific << setprecision(5);
        vout << " " << setw(12) << Beads[b].dAdAlpha;
        vout << " " << setw(12) << Beads[b].A;
        if( Beads[b].GetClientID() > 0 ){
            vout << " " << setw(4) << Beads[b].GetClientID();
        } else {
            vout << " " << setw(4) << "----";
        }
        vout << setw(8) << Beads[b].NumOfUpdates;
        vout << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b].Pos[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b].PMF[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b].dCV[i];
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b].pPMF[i];
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CBeadList::PrintPathUpdate(std::ostream& vout)
{
    int num_of_updates = 0;
    for(int b=0; b < NumOfBeads; b++){
        num_of_updates += Beads[b].NumOfUpdates;
    }

    if( num_of_updates > 0 ){
        // optimize path and alphas for current position
        for(int b=0; b < NumOfBeads; b++){
            Beads[b].PPos = Beads[b].RPos;
        }
        OptimizePath(Beads);
    } else {
        for(int b=0; b < NumOfBeads; b++){
            Beads[b].Alpha = 0.0;
        }
    }

    vout << "# === [PATH UPDATE] ============================================================" << endl;
    vout << "# Path name       = " << PathName << endl;
    vout << "# Number of CVs   = " << NumOfCVs << endl;
    vout << "# Number of beads = " << NumOfBeads << endl;

// header --------------------
    // legends
    vout << "#  ID   Type  ST Nalpha CID Updates  ";
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
    vout << "# ---- ------ -- ------ --- -------";
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
    vout << left << "#      names                       " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    vout << endl;
    vout << left << "#      types                       " << right;
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetType();
    }
    vout << endl;
    vout << left << "#      min                         " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMinValue();
    }
    vout << endl;
    vout << left << "#      max                         " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMaxValue();
    }
    vout << endl;
    vout << left << "#      maxmov                      " << right << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        if( CVs[i].GetMaxMovement() > 0 ){
            vout << " " << setw(12) << CVs[i].GetMaxMovement();
        } else {
            vout << " " << setw(12) << "--";
        }
    }
    vout << endl;

    vout << "# ---- ------ -- ------ --- -------";
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

    vout << "#    1      2  3      4   5       6";
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
    vout << "# ---- ------ -- ------ --- -------";
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
        if( Beads[b].Permanent ) {
            vout << "  " << setw(4) << Beads[b].GetBeadID() << setw(7) << " P     ";
        } else {
            vout << "  " << setw(4) << Beads[b].GetBeadID() << setw(7) << " F     ";
        }
        switch(Beads[b].GetMode()){
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
        vout << " " << setw(6) << Beads[b].Alpha;

        if( Beads[b].GetClientID() > 0 ){
            vout << setw(4) << Beads[b].GetClientID();
        } else {
            vout << setw(4) << " --";
        }
        vout << setw(8) << Beads[b].NumOfUpdates;
        vout << scientific << setprecision(5);
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b].Pos[i];
        }
        if( Beads[b].NumOfUpdates > 0 ){
            for(int i=0; i < NumOfCVs; i++){
                vout << " " << setw(12) << Beads[b].RPos[i];
            }
            for(int i=0; i < NumOfCVs; i++){
                vout << " " << setw(12) << Beads[b].RPos[i] - Beads[b].Pos[i];
            }
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CBeadList::PrintPath(std::ostream& vout)
{
    vout << "[PATH]" << endl;
    vout << "name     " << PathName << endl;
    vout << "ncvs     " << NumOfCVs << endl;
    vout << "nbeads   " << NumOfBeads << endl;
    vout << "names    ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetName();
    }
    vout << endl;
    vout << "types    ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetType();
    }
    vout << endl;
    vout << "min      ";
    vout << scientific << setprecision(5);
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMinValue();
    }
    vout << endl;
    vout << "max      ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMaxValue();
    }
    vout << endl;
    vout << "maxmov   ";
    for(int i=0; i < NumOfCVs; i++){
        vout << " " << setw(12) << CVs[i].GetMaxMovement();
    }
    vout << endl;

    for(int b=0; b < NumOfBeads; b++){
        vout << right;
        if( Beads[b].Permanent ) {
            vout << "permanent";
        } else {
            vout << "flexible ";
        }
        for(int i=0; i < NumOfCVs; i++){
            vout << " " << setw(12) << Beads[b].Pos[i];
        }
        vout << endl;
    }
}

//------------------------------------------------------------------------------

void CBeadList::SplitString(string text,vector<string>& words)
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

int CBeadList::GetNumOfBeads(void)
{
    return(NumOfBeads);
}

//------------------------------------------------------------------------------

int CBeadList::GetNumOfBeadsInRendezvousState(void)
{
    int count = 0;
    ProcessingMutex.Lock();

    for(int i=0; i < NumOfBeads; i++){
        if( Beads[i].GetMode() == BMO_WAITFORRENDEZVOUS ) count++;
    }

    ProcessingMutex.Unlock();
    return(count);
}

//------------------------------------------------------------------------------

ESTMState CBeadList::GetSTMStatus(void)
{
    return(STMStatus);
}

//------------------------------------------------------------------------------

bool CBeadList::IsAsynchronous(void)
{
    return(AsynchronousMode);
}

//------------------------------------------------------------------------------

CBead* CBeadList::GetBead(int bead_id)
{
    if( (bead_id <= 0) || (bead_id > NumOfBeads) ){
        LOGIC_ERROR("bead_id out-of-legal range");
    }
    CBead* p_bead = &Beads[bead_id-1];
    return(p_bead);
}

//------------------------------------------------------------------------------

CBead* CBeadList::GetBeadByClientID(int client_id)
{
    for(int b=0; b < NumOfBeads; b++){
        CBead* p_bead = &Beads[b];
        if( p_bead->GetClientID() == client_id ){
            return(p_bead);
        }
    }
    return(NULL);
}

//------------------------------------------------------------------------------

CBead* CBeadList::GetNextFreeBead(void)
{
    CBead* p_bead = NULL;
    for(int id=0; id < NumOfBeads; id++){
        if( Beads[id].GetClientID() == -1 ){
            Beads[id].SetClientID(0);
            p_bead = &Beads[id];
            return(p_bead);
        }
    }
    return(NULL);
}

//------------------------------------------------------------------------------

int CBeadList::GetSTMStep(void)
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

void CBeadList::PrintSTMHeader(void)
{
    vout << endl;
    vout << "# Step   Path length  Length change   Max movement  BID  Ave movement  Term" << endl;
    vout << "# ---- -------------- -------------- -------------- --- -------------- ----" << endl;

    HeaderPrinted = true;
}

//------------------------------------------------------------------------------

void CBeadList::PrintSTMStepInfo(void)
{
    vout << setw(6)  << STMStep << " ";
    vout << setw(14) << setprecision(7) << scientific << CurrentPathLength << " ";
    vout << setw(14) << setprecision(7) << scientific << UpdatedPathLength-CurrentPathLength << " ";

    MaxMovement = 0;
    AveMovement = 0;
    for(int b=0; b < NumOfBeads; b++){
        double mov = 0;
        for(int i=0; i < NumOfCVs; i++){
            mov += (Beads[b].RPos[i]-Beads[b].Pos[i])*(Beads[b].RPos[i]-Beads[b].Pos[i]);
        }
        AveMovement += mov; // add mov square
        mov = sqrt(mov);
        if( mov > MaxMovement ){
            MaxMovement = mov;
            MaxMovementBead = b+1;
        }
    }
    AveMovement = sqrt(AveMovement / NumOfBeads);

    vout << setw(14) << setprecision(7) << scientific << MaxMovement << " ";
    vout << setw(3) << MaxMovementBead << " ";
    vout << setw(14) << setprecision(7) << scientific << AveMovement << " ";

    int termcrit = 0;
    if( AveMovement < FinalAveMovement ) termcrit++;
    if( MaxMovement < FinalMaxMovement ) termcrit++;
    if( fabs(UpdatedPathLength-CurrentPathLength) < FinalMaxPLenChange ) termcrit++;

    vout << " " << setw(1) << termcrit << "/" << "3";
    vout << endl;

    if( termcrit == 3 ){
        STMStatus = ESTMS_PATH_FOUND;
        vout << endl;
        vout << ">> INFO: The path have converged." << endl;
        if( ProdPeriod > 0 ){
            vout << ">> INFO: Entering production accumulation (" << ProdPeriod <<" steps)." <<  endl;
        }
    }

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

void CBeadList::UpdateAllPositions(void)
{
    STMStep++;

    // vout << debug << "Updating positions ..." << endl << high;

    // reoptimize path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b].PPos = Beads[b].Pos;
    }
    CurrentPathLength = OptimizePath(Beads);

    for(int i=0; i < NumOfBeads; i++){
        CBead* p_bead = &Beads[i];
        p_bead->CalcProjector();
        p_bead->UpdatePosition();
    }
}

//------------------------------------------------------------------------------

void CBeadList::SmoothAllPositions(void)
{   
    if( (SmoothInterval == 0) || (STMStep % SmoothInterval != 0) ){
        for(int b=0; b < NumOfBeads; b++) {
            Beads[b].SPos = Beads[b].NPos;
        }
        return;
    }

    // vout << debug << "Smoothing positions ..." << endl << high;

    // smooth path
    for(int i=0; i < NumOfBeads; i++){
        CBead* p_bead = &Beads[i];

        if( (i == 0) || (i == NumOfBeads-1) || (p_bead->Permanent) ){
            for(int i=0; i < NumOfCVs; i++){
                p_bead->SPos[i] = p_bead->NPos[i];
            }
        } else {
            CBead* p_pnb = &Beads[i-1];
            CBead* p_nnb = &Beads[i+1];
            for(int i=0; i < NumOfCVs; i++){
                p_bead->SPos[i] = (1.0-SmoothingFac)*p_bead->NPos[i]
                                + 0.5*SmoothingFac*(p_pnb->NPos[i]+p_nnb->NPos[i]);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CBeadList::ReparametrizeAllPositions(void)
{   
    if( (ReparamInterval == 0) || (STMStep % ReparamInterval != 0) ){
        for(int b=0; b < NumOfBeads; b++) {
            Beads[b].RPos = Beads[b].SPos;
        }
        return;
    }

   // vout << debug << "Reparametrizing positions ..." << endl << high;

    // reoptimize path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b].PPos = Beads[b].SPos;
    }
    OptimizePath(Beads);

    // and correct positions
    for(int i=0; i < NumOfCVs; i++){
        Beads[0].RPos[i] = CVSplines[i].GetCV(0.0);
        Beads[NumOfBeads-1].RPos[i] = CVSplines[i].GetCV(1.0);
        for(int b=1; b < NumOfBeads-1; b++){
            double alpha = (double)b / ((double)NumOfBeads-1.0);
            Beads[b].RPos[i] = CVSplines[i].GetCV(alpha);
        }
    }
}

//------------------------------------------------------------------------------

void CBeadList::CheckBoundaries(void)
{
    for(int b=0; b < NumOfBeads; b++){
        CBead* p_bead = &Beads[b];

        for(int i=0; i < NumOfCVs; i++){
            if( p_bead->RPos[i] < CVs[i].GetMinValue() ){
                p_bead->RPos[i] = CVs[i].GetMinValue();
            }
            if( p_bead->RPos[i] > CVs[i].GetMaxValue() ){
                p_bead->RPos[i] = CVs[i].GetMaxValue();
            }
        }
    }

    // get data about the final path
    for(int b=0; b < NumOfBeads; b++){
        Beads[b].PPos = Beads[b].RPos;
    }
    UpdatedPathLength = OptimizePath(Beads);
}

//------------------------------------------------------------------------------

void CBeadList::CalculatePathData(void)
{
    // optimize path and alphas for current position
    for(int b=0; b < NumOfBeads; b++){
        Beads[b].PPos = Beads[b].Pos;
    }
    OptimizePath(Beads);

    // get PMF projections along path
    double fes = 0.0;
    for(int b=0; b < NumOfBeads; b++){
        double a = 0.0;
        // get bead derivative along path
        for(int i=0; i < NumOfCVs; i++) {
            Beads[b].dCV[i] = CVSplines[i].GetCVFirstDer(Beads[b].Alpha);
            a += Beads[b].dCV[i]*Beads[b].PMF[i];
        }
        Beads[b].dAdAlpha = a;
        if( b > 0 ){
            fes += 0.5*(Beads[b].Alpha - Beads[b-1].Alpha)*(Beads[b].dAdAlpha + Beads[b-1].dAdAlpha);
        }
        Beads[b].A = fes;
    }
}

//------------------------------------------------------------------------------

void CBeadList::SetServerTerminated(void)
{
    STMStatus = ESTMS_COMPLETED;
    RendezvousCond.BroadcastSignal();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CBeadList::OptimizePath(CSimpleVector<CBead>& beads)
{
    if( beads.GetLength() < 2 ){
        RUNTIME_ERROR("beads.GetLength() must be greater or equal 2");
    }

    // get initial path length from linear interpolation
    double tot_length = 0;
    for(int b=1; b < beads.GetLength(); b++){
        double slen = 0;
        for(int i=0; i < NumOfCVs; i++){
            slen += (beads[b].PPos[i]-beads[b-1].PPos[i])*(beads[b].PPos[i]-beads[b-1].PPos[i]);
        }
        tot_length += sqrt(slen);
    }

    if( tot_length == 0 ){
        RUNTIME_ERROR("path has zero length");
    }

    // get initial alphas from linear interpolation
    beads[0].Alpha = 0.0;
    double path_length = 0;
    for(int b=1; b < beads.GetLength()-1; b++){
        double slen = 0;
        for(int i=0; i < NumOfCVs; i++){
            slen += (beads[b].PPos[i]-beads[b-1].PPos[i])*(beads[b].PPos[i]-beads[b-1].PPos[i]);
        }
        if( slen == 0 ){
            RUNTIME_ERROR("path segment has zero length");
        }
        path_length += sqrt(slen);
        beads[b].Alpha = path_length/tot_length;
    }
    beads[beads.GetLength()-1].Alpha = 1.0;

   // vout << debug;
   // vout << "Initial path length = " << tot_length << endl;

    double prev_length = 0;

    for(;;){
        // interpolate CVS
        for(int i=0; i < NumOfCVs; i++){
            CVSplines[i].Allocate(beads.GetLength());
            for(int b=0; b < beads.GetLength(); b++){
                CVSplines[i].AddPoint(b,beads[b].Alpha,beads[b].PPos[i]);
            }
            CVSplines[i].Finalize();
        }

        prev_length = tot_length;

        // determine new path length
        tot_length = 0;
        for(int b=1; b < beads.GetLength(); b++){
            tot_length += GetSegmentLength(beads[b-1].Alpha,beads[b].Alpha);
        }

       // vout << "Optimized path length = " << tot_length << endl;

        if( fabs(tot_length-prev_length) < 1e-7 ){
       //     vout << "Converged path length = " << tot_length << endl;
            return(tot_length);
        }

        // determine new alphas
        beads[0].Alpha = 0.0;
        double path_length = 0;
        double prev_alpha = beads[0].Alpha;
        for(int b=1; b < beads.GetLength()-1; b++){
            path_length += GetSegmentLength(prev_alpha,beads[b].Alpha);
            beads[b].Alpha = path_length/tot_length;
            prev_alpha = beads[b].Alpha;
        }
        beads[beads.GetLength()-1].Alpha = 1.0;
    }
}

//------------------------------------------------------------------------------

double CBeadList::GetSegmentLength(double alpha1,double alpha2)
{
    double len = 0;
    for(int i=0; i < NumOfCVs; i++){
        SPos[i] = CVSplines[i].GetCV(alpha1);
    }
    double step = (alpha2-alpha1)/SegmentDiscretization;
    double alpha = alpha1 + step;
    while( alpha < alpha2 ){
        double slen2 = 0;
        for(int i=0; i < NumOfCVs; i++){
            double curr = CVSplines[i].GetCV(alpha);
            slen2 +=  (curr-SPos[i])*(curr-SPos[i]);
            SPos[i] = curr;
        }
        len += sqrt(slen2);
        alpha += step;
    }

    double slen2 = 0;
    for(int i=0; i < NumOfCVs; i++){
        double last = CVSplines[i].GetCV(alpha2);
        slen2 +=  (last-SPos[i])*(last-SPos[i]);
    }
    len += sqrt(slen2);

    return(len);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

