#ifndef AUSEngineCST_AH
#define AUSEngineCST_AH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <GPREngineAUS.hpp>
#include <IntegratorGPR.hpp>
#include <SmootherGPR.hpp>

//------------------------------------------------------------------------------

/** \brief base for AUS GPR engines
*/

class PMF_PACKAGE CAUSEngineCST_A: public CGPREngineAUS {
public:
// constructor and destructor -------------------------------------------------
    CAUSEngineCST_A(void);
    virtual ~CAUSEngineCST_A(void);

// setup input data ------------------------------------------------------------
    /// set accumulator
    virtual void SetAccumulator(CPMFAccumulatorPtr accu);

// setup -----------------------------------------------------------------------
    /// set include error
    virtual void SetIncludeError(bool set);

    /// skip energy calculation, it also disables errors
    virtual void SetNoEnergy(bool set);

    /// should we include glued area to energy calculation?
    virtual void IncludeGluedAreas(bool set);

    /// calc hyprms grd
    virtual void PrepForHyprmsGrd(bool set);

    /// calc logpl
    virtual void SetCalcLogPL(bool set);

// base methods ----------------------------------------------------------------
    /// run GPR
    virtual bool RunGPR(CVerboseStr& vout,bool nostat=false);

    /// return number of tasks
    virtual int GetNumOfTasks(void);

    /// write file with derivatives
    virtual bool WriteMFInfo(const CSmallString& name,int task);

    /// get log of Marginal Likelihood
    virtual double GetLogML(void);

    /// get derivative of logML wrt hyperparameters
    /// order sigmaf2, covar, wfac, ncorr, sigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    virtual void GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    virtual double GetLogPL(void);

    /// get derivative of logPL wrt hyperparameters
    /// order sigmaf2, covar, wfac, ncorr, sigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    virtual void GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

// section of protected data ---------------------------------------------------
protected:
    CIntegratorGPR  A;  // <lam>
    CIntegratorGPR  B;  // Cov(lam,Etot)
    CSmootherGPR    C;  // d<Etot>_FW
    CSmootherGPR    D;  // d<Etot>
    CSmootherGPR    E;  // -TdS{CST}corr

    CEnergySurfacePtr  A_ES;
    CEnergySurfacePtr  B_ES;
    CEnergySurfacePtr  C_ES;
    CEnergySurfacePtr  D_ES;
    CEnergySurfacePtr  E_ES;

    bool NoEnergy;

    void CalculateEnergy(CVerboseStr& vout);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CAUSEngineCST_A>    CAUSEngineCST_APtr;

//------------------------------------------------------------------------------

#endif
