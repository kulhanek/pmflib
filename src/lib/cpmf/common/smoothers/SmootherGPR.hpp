#ifndef SmootherGPRH
#define SmootherGPRH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <FortranMatrix.hpp>
#include <stddef.h>
#include <EnergySurface.hpp>
#include <EnergyProxy.hpp>
#include <IntegratorGPR.hpp>
#include <vector>

//------------------------------------------------------------------------------

/** \brief smooth energy by GPR
*/

class PMF_PACKAGE CSmootherGPR : public CGPREngine {
public:
// constructor and destructor -------------------------------------------------
    CSmootherGPR(void);
    virtual ~CSmootherGPR(void);

// setup methods --------------------------------------------------------------
    /// set accumulator
    void SetAccumulator(CPMFAccumulatorPtr accu);

    /// set input energy proxy
    void SetInputEnergyProxy(CEnergyProxyPtr p_proxy);

    /// set output energy surface
    void SetOutputES(CEnergySurfacePtr p_surf);

// setup
    /// skip energy calculation, it also disables errors
    virtual void SetNoEnergy(bool set);

    /// set include error
    virtual void SetIncludeError(bool set);

// execution method -----------------------------------------------------------
    /// run GPR
    virtual bool RunGPR(CVerboseStr& vout,bool nostat=false);

    /// interpolate data
    bool Interpolate(CVerboseStr& vout,bool nostat=false);

    /// prepare for subsequent call WriteMFInfo
    void PrepForMFInfo(void);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

// GPR model optimization ----------------------------------------------------
    /// calc hyprms grd
    virtual void PrepForHyprmsGrd(bool set);

    /// calc logpl
    virtual void SetCalcLogPL(bool set);

    /// get log of Marginal Likelihood
    virtual double GetLogML(void);

    /// get derivative of logML wrt hyperparameters
    /// order sigmaf2, wfac, nsigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    virtual void GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    virtual double GetLogPL(void);

    /// get derivative of logPL wrt hyperparameters
    /// order sigmaf2, wfac, nsigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    virtual void GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

// section of private data ----------------------------------------------------
private:
    CEnergyProxyPtr         EneProxy;
    CEnergySurfacePtr       EneSurface;

    // GPR data, sizes and index maps
    size_t                  GPRSize;
    std::vector<size_t>     SampledMap;
    size_t                  NumOfValues;
    std::vector<size_t>     ValueMap;

    // GPR model
    EGPRKernel              Kernel;
    CFortranMatrix          KS;             // kernel matrix
    double                  logdetK;
    double                  Mean;
    CSimpleVector<double>   Y;              // energy
    CSimpleVector<double>   GPRModel;       // weights
    CFortranMatrix          Cov;            // covariances

    // derivatives
    CFortranMatrix          Kder;           // derivative of kernels w.r.t. a hyperparameter

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);
    void CalculateCovs(CVerboseStr& vout);
    void CalculateErrorsFromCov(CVerboseStr& vout);

    double GetValue(const CSimpleVector<double>& position);

    // kernel matrix + noise
    void CreateKS(void);
    void CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& ky);

    // derivatives
    void CalcKderWRTSigmaF2(void);
    void CalcKderWRTWFac(size_t cv);
    void CalcKderWRTSigmaN2(void);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CSmootherGPR>    CSmootherGPRPtr;

//------------------------------------------------------------------------------

#endif
