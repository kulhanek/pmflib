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
    ERBFLLS_SVD   = 1,
    ERBFLLS_QR    = 2,    // we have overdetermined system, thus only QR is applicable
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
    void SetInputABFAccumulator(const CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

    /// multiply of bin sizes
    void SetWFac1(double wfac);

    /// multiply of bin sizes, if zero use wfac1
    void SetWFac2(double wfac);

    /// set rcond for SVD
    void SetRCond(double rcond);

    /// set reduction factor for RBF
    void SetRFac1(double rfac);

    /// set reduction factor for RBF, if zero use rfac1
    void SetRFac2(double rfac);

    /// should we apply periodicity?
    void SetPeriodicity(bool set);

    /// set overhang, i.e. how many RBFs will be deposited outside of ABF accumulator
    void SetOverhang(int nrbfs);

    /// set algorithm for LLS
    void SetLLSMehod(ERBFLLSMethod set);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout,bool errors);

    /// get root mean square residuals
    double GetRMSR(void);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

// section of private data ----------------------------------------------------
private:
    const CABFAccumulator*  Accumulator;
    CEnergySurface*         FES;

    double                  WFac1;
    double                  RFac1;   // reduction factor for number of RBFs
    double                  WFac2;
    double                  RFac2;   // reduction factor for number of RBFs
    int                     Overhang;
    ERBFLLSMethod           Method;
    EABFAccuValue           IntegratedRealm;

    // RBF data
    int                     NumOfRBFs;
    CSimpleVector<double>   Weights;
    CSimpleVector<int>      NumOfRBFBins;
    int                     NumOfCVs;
    CSimpleVector<double>   Sigmas;

    // SVD setup
    double                  RCond;

    bool IntegrateByLS(CVerboseStr& vout);

    void    GetRBFPosition(unsigned int index,CSimpleVector<double>& position);
    double  GetValue(const CSimpleVector<double>& position);
};

//------------------------------------------------------------------------------

#endif
