#ifndef BeadH
#define BeadH
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
// =============================================================================

#include <XMLElement.hpp>
#include <SimpleVector.hpp>
#include <FortranMatrix.hpp>
#include <SimpleMutex.hpp>

//------------------------------------------------------------------------------

#define BMO_UNKNOWN             0
#define BMO_INITIALIZATION      1
#define BMO_ACCUMULATION        2
#define BMO_EQUILIBRATION       3
#define BMO_PRODUCTION          4
#define BMO_WAITFORRENDEZVOUS   5
#define BMO_TERMINATE           6

#define BMS_PREPARED            0
#define BMS_RUNNING             1
#define BMS_FINISHED            2

//------------------------------------------------------------------------------

class CBeadList;

//------------------------------------------------------------------------------

class CBead {
public:
    CBead(void);

// informative methods ---------------------------------------------------------
    /// get client ID
    int GetClientID(void);

    /// get bead ID
    int GetBeadID(void);

    // get bead position
    double GetPos(int cv);

    /// get bead mode
    int GetMode(void);

    /// get bead mode
    char GetModeString(void);

    /// get mode status
    int GetModeStatus(void);

    /// get length of current mode
    int GetModeLength(void);

// executive methods -----------------------------------------------------------
    /// set base list
    void InitBead(CBeadList* p_list,int ncvs);

    /// set bead data
    void SetBeadData(int beadid,const CSimpleVector<double>& pos,bool flexible);

    /// set client ID
    void SetClientID(int client_id);

    /// relase possibly crashed bead
    void ReleaseBead(void);

    /// next program
    void MoveToNextMode(void);

    /// calc projector - path must be optimized!
    void CalcProjector(void);

    /// update bead position
    void UpdatePosition(void);

// input/output methods --------------------------------------------------------
    /// load data
    void LoadInfo(CXMLElement* p_ele);

    /// save data
    void SaveInfo(CXMLElement* p_ele);

    /// get production data
    void GetProductionData(CXMLElement* p_ele);

    /// skip production data
    void SkipProductionData(void);

    /// wait for randezvous
    void WaitForRendezvous(void);

    /// get production data
    void SetNextStepData(CXMLElement* p_ele);

// section of private data -----------------------------------------------------
private:
    CBeadList*              BeadList;

    // bead data
    int                     BeadID;         // bead id
    int                     ClientID;       // client id
    int                     Mode;           // current bead mode
    int                     ModeStatus;     // what is status of current mode

    // bead data
    int                     NumOfCVs;       // number of CVs
    bool                    Permanent;      // is bead permanent
    CSimpleVector<double>   Pos;            // bead position
    CSimpleVector<double>   NPos;           // new bead position
    CSimpleVector<double>   SPos;           // smoothed position
    CSimpleVector<double>   RPos;           // reparametrized position
    CSimpleVector<double>   PPos;           // position for path optimization
    CSimpleVector<double>   PMF;            // force acting on the bead
    CSimpleVector<double>   dCV;            // path derivatives
    CSimpleVector<double>   pPMF;           // force/velocity acting perpendiculary to the path
    CFortranMatrix          MTZ;            // metric tensor
    CFortranMatrix          P;              // projector
    double                  Alpha;          // path position
    double                  dAdAlpha;       // free energy derivative
    double                  A;              // free energy
    int                     NumOfUpdates;   // how many updates was performed

    friend class CBeadList;
};

//------------------------------------------------------------------------------

#endif
