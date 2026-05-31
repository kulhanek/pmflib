#ifndef BeadH
#define BeadH
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
// =============================================================================

#include <PMFMainHeader.hpp>
#include <XMLElement.hpp>
#include <SimpleVector.hpp>
#include <FortranMatrix.hpp>
#include <SimpleMutex.hpp>
#include <memory>

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

class CSTMPath;

//------------------------------------------------------------------------------

class PMF_PACKAGE CBead {
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
    void InitBead(CSTMPath* p_list,int ncvs);

    /// set bead data
    void SetBeadData(int beadid,const CSimpleVector<double>& pos,bool flexible);

    /// set client ID
    void SetClientID(int client_id);

    /// release possibly crashed bead
    void ReleaseBead(void);

    /// next program
    void MoveToNextMode(void);

    /// reset position updates
    void ResetPosUpdates(void);

    /// calculate projector - path must be optimized!
    void CalcProjector(void);

    /// update bead position - gradient descent
    void UpdatePositionGD(double step);

    /// update bead position - normalized gradient descent
    void UpdatePositionNGD(double step,double mingnormeps);

    /// update bead position - normalized gradient descent vs gradient descent
    void UpdatePositionNGDAuto(double step,double maxgnorm,double mingnormeps);

    /// ResetADAM
    void ResetADAM(void);

    /// Adam (Adaptive Moment Estimation)
    void UpdatePositionADAM(double step,double beta1,double beta2,double mingnormeps);

    /// ADABelif
    void UpdatePositionADABelif(double step,double beta1,double beta2,double mingnormeps);

    /// AMSGrad
    void UpdatePositionAMSGrad(double step,double beta1,double beta2,double mingnormeps);

    /// AMSGrad + Bias Corrected Variant
    void UpdatePositionAMSGradBC(double step,double beta1,double beta2,double mingnormeps);

    /// currently only for ADAM
    void UpdatePositionFinalize(void);

// input/output methods --------------------------------------------------------
    /// load data
    void LoadInfo(CXMLElement* p_ele);

    /// save data
    void SaveInfo(CXMLElement* p_ele);

    /// get production data
    void GetProductionData(CXMLElement* p_ele);

    /// skip production data
    void SkipProductionData(void);

    /// wait for rendezvous
    void WaitForRendezvous(void);

    /// get production data
    void SetNextStepData(CXMLElement* p_ele);

// section of private data -----------------------------------------------------
private:
    CSTMPath*               BeadList;

// bead data
    int                     BeadID;         // bead id
    int                     ClientID;       // client id
    int                     Mode;           // current bead mode
    int                     ModeStatus;     // what is status of current mode

// bead data 
    int                     NumOfCVs;       // number of CVs
    bool                    Permanent;      // is bead permanent?
    int                     NumOfUpdates;   // how many updates was performed

// bead data - unscaled
    CSimpleVector<double>   BPos;           // bead position
    CSimpleVector<double>   MF;             // force acting on the bead
    CFortranMatrix          MTZ;            // metric tensor

// bead data - free energy
    double                  Alpha;          // path position
    double                  dAdAlpha;       // free energy derivative
    double                  A;              // free energy

// bead data - scaled
    CFortranMatrix          P;              // projector
    CSimpleVector<double>   Pos;            // bead position
    CSimpleVector<double>   pMF;            // force acting perpendicularly to the path
    CSimpleVector<double>   dCVdAlpha;

// helper positions - scaled
    CSimpleVector<double>   OPos;           // old bead position
    CSimpleVector<double>   NPos;           // new bead position
    CSimpleVector<double>   SPos;           // smoothed position
    CSimpleVector<double>   FPos;           // re-parametrized position
    CSimpleVector<double>   PPos;           // position for path optimization

    // Adam (Adaptive Moment Estimation) variants
    double                  beta1told;
    double                  beta2told;
    CSimpleVector<double>   mtold;
    CSimpleVector<double>   vtold;
    CSimpleVector<double>   vthatold;

    double                  beta1tnew;
    double                  beta2tnew;
    CSimpleVector<double>   mtnew;
    CSimpleVector<double>   vtnew;
    CSimpleVector<double>   vthatnew;

    friend class CSTMPath;
};

//------------------------------------------------------------------------------

typedef std::shared_ptr<CBead>   CBeadPtr;

//------------------------------------------------------------------------------

#endif
