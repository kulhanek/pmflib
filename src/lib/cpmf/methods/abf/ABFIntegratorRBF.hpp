#ifndef ABFIntegratorRBFH
#define ABFIntegratorRBFH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <ABFAccumulator.hpp>

//------------------------------------------------------------------------------

class CABFAccumulator;
class CEnergySurface;

// how to find solution for least square problem

enum ERBFLLSMethod {
    ERBFLLS_SVD     = 1,    // SVD - divide and conquer driver
    ERBFLLS_QR      = 2,    // we have overdetermined system, thus only QR is applicable
};

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing radial basis functions
*/

class PMF_PACKAGE CABFIntegratorRBF {
public:
// constructor and destructor -------------------------------------------------
    CABFIntegratorRBF(void);
    virtual ~CABFIntegratorRBF(void);

// setup methods --------------------------------------------------------------
    /// set input ABF accumulator, only ABF forces are integrated
    void SetInputABFAccumulator(CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

    /// multiply of bin sizes
    void SetWFac(const CSmallString& spec);

    /// set rcond for SVD
    void SetRCond(double rcond);

    /// set reduction factor for RBF
    void SetRFac(const CSmallString& spec);

    /// should we apply periodicity?
    void SetPeriodicity(bool set);

    /// set overhang, i.e. how many RBFs will be deposited outside of ABF accumulator
    void SetOverhang(int nrbfs);

    /// set algorithm for LLS
    void SetLLSMehod(ERBFLLSMethod set);

    /// should we include glued area to energy calculation?
    void IncludeGluedAreas(bool set);

    /// skip energy calculation, it also disables errors
    void SetNoEnergy(bool set);

    /// set position of global minimum
    void SetGlobalMin(const CSmallString& spec);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout);

    /// get root mean square residuals
    double GetRMSR(size_t cv);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

// this destroys the state of the integrator
    /// filter by MF outliers
    void FilterByMFZScore(double zscore,CVerboseStr& vout);

    /// print exec info
    static void PrintExecInfo(CVerboseStr& vout);

// section of private data ----------------------------------------------------
private:
    CABFAccumulator*        Accumulator;
    CEnergySurface*         FES;

    CSimpleVector<double>   WFac;
    CSimpleVector<double>   RFac;   // reduction factor for number of RBFs
    int                     Overhang;
    ERBFLLSMethod           Method;

    // RBF data
    size_t                  NCVs;
    size_t                  NumOfBins;
    size_t                  NumOfRBFs;
    CSimpleVector<double>   Weights;
    CSimpleVector<size_t>   NumOfRBFBins;
    CSimpleVector<double>   Sigmas;

    bool                    IncludeGluedBins;

    bool                    NoEnergy;
    bool                    GlobalMinSet;
    CSimpleVector<double>   GPos;           // global position, either detected or use

    size_t                  NumOfValues;
    CSimpleVector<size_t>   ValueMap;

    // SVD setup
    double                  RCond;

    bool IntegrateByLS(CVerboseStr& vout);

    void    GetRBFPosition(size_t index,CSimpleVector<double>& position);
    double  GetValue(const CSimpleVector<double>& position);
    double  GetMeanForce(const CSimpleVector<double>& position,size_t icoord);
    void    CalculateEnergy(CVerboseStr& vout);
};

//------------------------------------------------------------------------------

#endif
