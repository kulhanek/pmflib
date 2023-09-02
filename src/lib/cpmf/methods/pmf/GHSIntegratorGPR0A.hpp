#ifndef GHSIntegratorGPR0AH
#define GHSIntegratorGPR0AH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <EnergyDerProxy.hpp>
#include <EnergySurface.hpp>
#include <GPRKernel.hpp>
#include <GPRHyprms.hpp>

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CGHSIntegratorGPR0A : public CGPRHyprms {
public:
// constructor and destructor -------------------------------------------------
    CGHSIntegratorGPR0A(void);
    virtual ~CGHSIntegratorGPR0A(void);

// setup methods --------------------------------------------------------------
    /// set accumulator
    void SetAccumulator(CPMFAccumulatorPtr accu);

    /// set input derivative proxies
    void SetGDerProxy(CEnergyDerProxyPtr p_proxy);
    void SetHEneProxy(CEnergyProxyPtr p_proxy);
    void SetSDerProxy(CEnergyDerProxyPtr p_proxy);

    /// set output energy surfaces
    void SetOutputFES(CEnergySurfacePtr p_surf);
    void SetOutputHES(CEnergySurfacePtr p_surf);
    void SetOutputSES(CEnergySurfacePtr p_surf);

// setup
    /// set include error
    void SetIncludeError(bool iset);

    /// skip energy calculation, it also disables errors
    void SetNoEnergy(bool iset);

    /// calc hyprms grd
    void PrepForHyprmsGrd(bool iset);

    /// calc logpl
    void SetCalcLogPL(bool iset);

    /// prepare for subsequent call WriteMFInfo
    void PrepForMFInfo(void);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout,bool nostat=false);

    /// get log of Marginal Likelihood
    double GetLogML(void);

    /// get derivative of logML wrt hyperparameters
    /// order sigmaf2, covar, wfac, nsigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    void GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    double GetLogPL(void);

    /// get derivative of logPL wrt hyperparameters
    /// order sigmaf2, covar, wfac, nsigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    void GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name,int task);

// section of private data ----------------------------------------------------
private:
    CEnergyDerProxyPtr      GDerProxy;
    CEnergyProxyPtr         HEneProxy;
    CEnergyDerProxyPtr      SDerProxy;
    CEnergySurfacePtr       GSurface;
    CEnergySurfacePtr       HSurface;
    CEnergySurfacePtr       SSurface;

    // GPR data, sizes and index maps
    size_t                  GPRSize;
    size_t                  NumOfUsedBins;
    std::vector<size_t>     SampledMap;

    // setup
    bool                    NoEnergy;

    // GPR model
    double                  HMean;
    EGPRKernel              Kernel;
    CFortranMatrix          KS;             // kernel matrix with noise
    CFortranMatrix          TK;             // task covariances
    bool                    KSInverted;
    double                  logdetK;
    CSimpleVector<double>   Y;              // mean forces
    CSimpleVector<double>   GPRModel;       // weights
    CFortranMatrix          Cov;            // covariances

    // derivatives
    CFortranMatrix          Kder;           // derivative of kernels w.r.t. a hyperparameter
    CFortranMatrix          TKder;          // task covariances derivatives

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);

// output data
    double GetValue(const CSimpleVector<double>& position,int task);

// source data
    double GetTrainingValue(const CSimpleVector<double>& position,size_t icoord,int task);
    double GetTrainingValueVar(const CSimpleVector<double>& position,size_t icoord,int task);
    double GetRMSR(size_t cv,int task);

// GPR
    void CreateTK(void);
    void CreateKS(void);

    void CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& ky,int task);
    void CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& ky2,int task);

// derivatives
    void CalcKderWRTSigmaF2(size_t idx);
    void CalcKderWRTWFac(size_t cv);
    void CalcKderWRTSigmaN2(size_t idx);
    void CreateTKDerSigmaF2(size_t idx);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CGHSIntegratorGPR0A>    CGHSIntegratorGPR0APtr;

//------------------------------------------------------------------------------

#endif
