#ifndef IntegratorGPRH
#define IntegratorGPRH
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
#include <EnergyDerProxy.hpp>
#include <EnergySurface.hpp>
#include <GPRHyprms.hpp>

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CIntegratorGPR : public CGPRHyprms {
public:
// constructor and destructor -------------------------------------------------
    CIntegratorGPR(void);
    virtual ~CIntegratorGPR(void);

// setup methods --------------------------------------------------------------
    /// set accumulator
    void SetAccumulator(CPMFAccumulatorPtr accu);

    /// set input energy der proxy
    void AddInputEnergyDerProxy(CEnergyDerProxyPtr p_proxy);

    /// clear energy der proxies
    void ClearInputEnergyDerProxies(void);

    /// set output free energy surface
    void SetOutputES(CEnergySurfacePtr p_surf);

// setup
    /// set include error
    void SetIncludeError(bool set);

    /// skip energy calculation, it also disables errors
    void SetNoEnergy(bool set);

    /// should we include glued area to energy calculation?
    void IncludeGluedAreas(bool set);

    /// calc hyprms grd
    void PrepForHyprmsGrd(bool set);

    /// calc logpl
    void SetCalcLogPL(bool set);

    /// use fast error algorithm
    void SetFastError(bool set);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout,bool nostat=false);

    // get mean force
    double GetMeanForce(const CSimpleVector<double>& position,size_t icoord);

    // get mean force variance
    double GetMeanForceVar(const CSimpleVector<double>& position,size_t icoord);

    /// get root mean square residuals
    double GetRMSR(size_t cv);

    /// get log of Marginal Likelihood
    double GetLogML(void);

    /// get derivative of logML wrt hyperparameters
    /// order sigmaf2, wfac, ncorr, nsigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    void GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    double GetLogPL(void);

    /// get derivative of logPL wrt hyperparameters
    /// order sigmaf2, wfac, ncorr, nsigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    void GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// prepare for subsequent call WriteMFInfo
    void PrepForMFInfo(void);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

    /// remove mean force outliers from ABF data
    void FilterByMFZScore(double maxzscore,CVerboseStr& vout);

    /// perform statistical reweighting
    CEnergySurfacePtr ReduceFES(const std::vector<bool>& keepcvs);

// section of private data ----------------------------------------------------
private:
    std::vector<CEnergyDerProxyPtr> DerProxyItems;
    CEnergySurfacePtr               EneSurface;

    // GPR data, sizes and index maps
    size_t                  GPRSize;
    size_t                  NumOfUsedBins;
    std::vector<size_t>     SampledMap;
    std::vector<size_t>     DerProxyMap;
    size_t                  NumOfValues;
    std::vector<size_t>     ValueMap;

    // setup
    bool                    NoEnergy;
    bool                    IncludeError;
    bool                    FastErrors;     // use faster but more memory intensive algorithm
    bool                    IncludeGluedBins;

    // GPR model
    CFortranMatrix          KS;             // kernel matrix with noice
    bool                    KSInverted;
    double                  logdetK;
    CSimpleVector<double>   Y;              // mean forces
    CSimpleVector<double>   GPRModel;       // weights
    CFortranMatrix          Cov;            // covariances

    // derivatives
    CFortranMatrix          Kder;           // derivative of kernels w.r.t. a hyperparameter

    // statistical reweighting
    std::vector<bool>       IntegratedCVs;

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);
    void CalculateErrors(CVerboseStr& vout);
    void CalculateErrorsFromCov(CVerboseStr& vout);
    void CalculateCovs(CVerboseStr& vout);

    double GetValue(const CSimpleVector<double>& position);
    double GetVar(CSimpleVector<double>& lpos);

// optimized version var+cov
    void GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& llvar,double& lrcov);

// kernel matrix + noise
    void CreateKS(void);
    void CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& ky);
    void CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& ky2);

// derivatives
    void CalcKderWRTSigmaF2(void);
    void CalcKderWRTNCorr(void);
    void CalcKderWRTWFac(size_t cv);
    void CalcKderWRTSigmaN2(size_t cv);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CIntegratorGPR>    CIntegratorGPRPtr;

//------------------------------------------------------------------------------

#endif
