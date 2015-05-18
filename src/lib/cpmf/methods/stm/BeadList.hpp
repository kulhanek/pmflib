#ifndef BeadListH
#define BeadListH
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

#include "Bead.hpp"
#include <SimpleList.hpp>
#include <IndexCounter.hpp>
#include <SimpleMutex.hpp>
#include <SimpleCond.hpp>
#include <ostream>
#include <ColVariable.hpp>
#include <PrmFile.hpp>
#include <VerboseStr.hpp>
#include <CVSpline.hpp>

//------------------------------------------------------------------------------

enum ESTMState {
    ESTMS_INITIALIZED = 0,
    ESTMS_OPTIMIZING = 1,
    ESTMS_PATH_FOUND = 2,
    ESTMS_MAX_STEPS_REACHED = 3,
    ESTMS_COMPLETED = 4
};

//------------------------------------------------------------------------------

/// STM enegine

class CBeadList {
public:
    CBeadList(void);

// file support ----------------------------------------------------------------
    /// attach stream
    void AttachVerboseStream(std::ostream& str,bool verbose);

    /// load input path
    bool LoadPath(CPrmFile& file);

    /// save path
    bool SavePath(void);

    /// save path
    void SavePath(const CSmallString& name);

    /// save path summary
    bool SavePathSummary(void);

    /// save path summary
    void SavePathSummary(const CSmallString& name);

    /// flush path
    void FlushPath(void);

    /// open trajectory
    bool OpenTrajectory(void);

    /// save path
    void SaveSnapshot(void);

    /// close trajectory
    void CloseTrajectory(void);

    /// force termination
    void ForceTermination(void);

    /// set asynchrnous mode
    void SetAsynchronousMode(bool set);

// control file parsing --------------------------------------------------------
    /// load files setup
    void ProcessFilesControl(CPrmFile& file);

    /// load stm setup
    void ProcessSTMControl(CPrmFile& file);

    /// load intervals setup
    void ProcessIntervalsControl(CPrmFile& file);

    /// load path specification
    void ProcessPathControl(CPrmFile& file);

// network executive methods ---------------------------------------------------
    /// check if client is eligible to connect to STM server
    bool CheckClient(CXMLElement* p_cele);

    /// register client
    void RegisterBead(int bead_id,int client_id);

    /// move all beads to the first mode
    void BeginAsynchronousMode(void);

    /// exchange data with client
    void ExchangeData(CXMLElement* p_cele,CXMLElement* p_rele);

    /// process accumulated data
    void ProcessPathAsynchronously(void);

// i/o methods -----------------------------------------------------------------
    /// load entire path form XML stream
    void LoadInfo(CXMLElement* p_ele);

    /// save path to XML stream
    void SaveInfo(CXMLElement* p_ele);

    /// check CVs
    bool CheckCoords(CXMLElement* p_ele);

// information methods ---------------------------------------------------------
    /// get number of beads
    int GetNumOfBeads(void);

    /// return number of beads in rendezvous state
    int GetNumOfBeadsInRendezvousState(void);

    /// get STM status
    ESTMState GetSTMStatus(void);

    /// shoul server be terminated?
    bool IsAsynchronous(void);

    /// print summary
    void PrintPathSummary(std::ostream& vout);

    /// print path
    void PrintPath(std::ostream& vout);

    /// print path
    void PrintPathUpdate(std::ostream& vout);

    /// get bead by its ID (1,2,3,...,NumOfBeads)
    CBead* GetBead(int bead_id);

    /// get bead by client ID
    CBead* GetBeadByClientID(int client_id);

    /// get next free bead
    CBead* GetNextFreeBead(void);

    /// get current STM step
    int GetSTMStep(void);

// executive methods -----------------------------------------------------------
    /// update all positions
    void UpdateAllPositions(void);

    /// smooth all positions
    void SmoothAllPositions(void);

    /// reparametrize all positions
    void ReparametrizeAllPositions(void);

    /// check bead position boundaries
    void CheckBoundaries(void);

    /// calculate path data
    void CalculatePathData(void);

    /// server is terminated - unblock waiting beads
    void SetServerTerminated(void);

// section of private data -----------------------------------------------------
private:
    CSmallString                PathName;
    int                         NumOfCVs;
    CSimpleVector<CColVariable> CVs;        // CV definitions

    int                         NumOfBeads;
    CSimpleVector<CBead>        Beads;      // bead data
    CSimpleVector<CCVSpline>    CVSplines;  // interpolated CV

    // STM setup ---------------------------------
    int                 MaxSTMSteps;        // maximum number of STM steps
    double              StepSize;           // step size for bead update
    double              FinalMaxPLenChange; // termination criteria - max path movement
    double              FinalMaxMovement;   // termination criteria - max path movement
    double              FinalAveMovement;   // termination criteria - average path movement

    int                 InitPeriod;         // initialization period
    int                 AccuPeriod;         // accumulation period
    int                 EquiPeriod;         // equilibration period
    int                 ProdPeriod;         // final production period

    bool                AsynchronousMode;   // use asynchronous mode
    double              SmoothingFac;       // smoothing factor

    // intervals ---------------------------------
    int                 TrajInterval;       // how often to print snapshot to trajectory
    int                 OutInterval;        // how often to write current path
    int                 SmoothInterval;     // how often to smooth path
    int                 ReparamInterval;    // how often to reparametrize path

    // files -------------------------------------
    CSmallString        InputPath;
    CSmallString        OutputPath;
    CSmallString        OutputPathSummary;
    CSmallString        PathTrajectory;
    std::ofstream       Trajectory;

    CVerboseStr         vout;               // info channel
    int                 STMStep;            // current STM step
    double              MaxMovement;        // current max path movement
    int                 MaxMovementBead;    // current max path movement is for given bead
    double              AveMovement;        // current average path movement
    double              CurrentPathLength;
    double              UpdatedPathLength;  

    CSimpleMutex        ProcessingMutex;    // mutex for path processing accesses

    // STM engine state
    ESTMState           STMStatus;
    bool                HeaderPrinted;
    bool                Terminate;

    // rendezvous point for synchronous update
    CSimpleMutex        RendezvousMutex;
    CSimpleCond         RendezvousCond;
    int                 NumOfRendezvousBeads;

    // technicals
    CSimpleVector<double>   SPos;
    int                     SegmentDiscretization;  // how many subpoints are used to calculate segment length

    // allocate path
    void AllocatePath(void);

    // clear everything
    void ClearPath(void);

    // split text
    void SplitString(std::string text,std::vector<std::string>& words);

    // print header
    void PrintSTMHeader(void);

    // print stm step
    void PrintSTMStepInfo(void);

    // print path summary header
    void PrintPathSummaryHeader(std::ostream& vout);

    // print path summary data
    void PrintPathSummaryData(std::ostream& vout);

    // synchronous mode
    void ExchangeDataSynchronously(CXMLElement* p_cele,CXMLElement* p_rele);
    void ProcessProductionData(CBead* p_bead);

    // asynchronous mode
    void ExchangeDataAsynchronously(CXMLElement* p_cele,CXMLElement* p_rele);

    // terminate client
    void TerminateClient(CXMLElement* p_rele);

    // path optimization, it return path length and initialize CVSplines
    // positions that are optimized are in PPos
    double OptimizePath(CSimpleVector<CBead>& beads);

    // get path length
    double GetSegmentLength(double alpha1,double alpha2);

    // read path helpers
    void ReadPathControls(CPrmFile& file);
    int  ReadPathNumberOfUserBeads(CPrmFile& file);
    void ReadPathUserBeads(CPrmFile& file,CSimpleVector<CBead>& beads);

    friend class CBead;
};

//------------------------------------------------------------------------------

#endif
