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

//------------------------------------------------------------------------------

// linear algebra pathway

enum EGPRLAMethod {
    EGPRLA_LU      = 1,    // DGETRI/DGETRF (via LU factorization)
    EGPRLA_SVD     = 2,    // via SVD factorization, divide and conquer driver
    EGPRLA_SVD2    = 3,    // via SVD factorization, simple driver
    EGPRLA_LL      = 4,    // DPOTRF/DPOTRI (via LL factorization)
};

//------------------------------------------------------------------------------

// supported kernel functions
// ARD = automatic relevance determination

enum EGPRKernel {
    EGPRK_ARDSE   = 1,    // ARD squared exponential (radial basis function)
    EGPRK_ARDMC52 = 2,    // ARD Matern class 5/2 function
    EGPRK_ARDMC32 = 3,    // ARD Matern class 3/2 function
    EGPRK_ARDMC12 = 4,    // ARD Matern class 1/2 function
};

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CIntegratorGPR {
public:
// constructor and destructor -------------------------------------------------
    CIntegratorGPR(void);
    virtual ~CIntegratorGPR(void);

// setup methods --------------------------------------------------------------
    /// set input energy der proxy
    void AddInputEnergyDerProxy(CEnergyDerProxyPtr p_proxy);

    /// clear energy der proxies
    void ClearInputEnergyDerProxies(void);

    /// set output free energy surface
    void SetOutputES(CEnergySurfacePtr p_surf);

// hyperparameters
    /// set sigmaf2
    void SetSigmaF2(double sigf2);

    /// set ncorr
    void SetNCorr(double value);

    /// set sigman2
    void SetSigmaN2(const CSmallString& spec);

    /// set sigman2
    void SetSigmaN2(CSimpleVector<double>& sigman2);

    /// set sigman2
    void SetSigmaN2(size_t cvind, double value);

    /// multiply of bin sizes
    void SetWFac(const CSmallString& spec);

    /// multiply of bin sizes
    void SetWFac(CSimpleVector<double>& wfac);

    /// multiply of bin sizes
    void SetWFac(size_t cvind, double value);

    /// load hyperparameters from file
    void LoadGPRHyprms(const CSmallString& name);

// setup
    /// set include error
    void SetIncludeError(bool set);

    /// skip energy calculation, it also disables errors
    void SetNoEnergy(bool set);

    /// switch to numerical evaluation of kernel fce derivatives
    void SetUseNumDiff(bool set);

    /// set algorithm for LA
    void SetLAMethod(EGPRLAMethod set);

    /// set algorithm for LA
    void SetLAMethod(const CSmallString& method);

    /// set kernel
    void SetKernel(const CSmallString& kernel);

    /// set rcond for SVD
    void SetRCond(double rcond);

    /// should we include glued area to energy calculation?
    void IncludeGluedAreas(bool set);

    /// set position of global minimum - spec in real units
    void SetGlobalMin(const CSmallString& spec);

    /// set position of global minimum - in internal units
    void SetGlobalMin(const CSimpleVector<double>& pos);

    /// get position of global minima - in internal units
    CSimpleVector<double> GetGlobalMin(void);

    /// get position of global minima - bin
    int GetGlobalMinBin(void);

    /// use inversion alg
    void SetUseInv(bool set);

    /// calc hyprms grd
    void PrepForHyprmsGrd(bool set);

    /// calc logpl
    void SetCalcLogPL(bool set);

    /// include zero-point at user provided global minimum
    void SetUseZeroPoint(bool set);

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
    /// order sigmaf2, ncorr, wfac, only requested ders are calculated
    /// derivatives are ADDED to der
    void GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    double GetLogPL(void);

    /// get derivative of logPL wrt hyperparameters
    /// order sigmaf2, ncorr, wfac, only requested ders are calculated
    /// derivatives are ADDED to der
    void GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// prepare for subsequent call WriteMFInfo
    void PrepForMFInfo(void);

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

    /// remove mean force outliers from ABF data
    void FilterByMFZScore(double maxzscore,CVerboseStr& vout);

    /// print exec info
    void PrintExecInfo(CVerboseStr& vout);

    /// get kernel name
    const CSmallString GetKernelName(void);

    /// perform statistical reweighting
    CEnergySurfacePtr ReduceFES(const std::vector<bool>& keepcvs);

// section of private data ----------------------------------------------------
private:
    std::vector<CEnergyDerProxyPtr> DerProxyItems;
    CEnergySurfacePtr               EneSurface;
    int                             NumOfThreads;

    // GPR data, sizes and index maps
    size_t                  NumOfCVs;
    size_t                  NumOfBins;
    size_t                  GPRSize;
    size_t                  NumOfUsedBins;
    std::vector<size_t>     SampledMap;
    std::vector<size_t>     DerProxyMap;
    size_t                  NumOfValues;
    std::vector<size_t>     ValueMap;

    // setup
    bool                    UseNumDiff;
    bool                    NoEnergy;
    bool                    IncludeError;
    bool                    IncludeGluedBins;
    EGPRLAMethod            Method;

    // hyperparameters
    double                  SigmaF2;
    double                  NCorr;
    CSimpleVector<double>   SigmaN2;
    CSimpleVector<double>   WFac;
    CSimpleVector<double>   CVLengths2;

    // GPR model
    EGPRKernel              Kernel;
    CFortranMatrix          KS;             // kernel matrix with noice
    bool                    KSInverted;
    double                  logdetK;
    CSimpleVector<double>   Y;              // mean forces
    CSimpleVector<double>   GPRModel;       // weights
    CFortranMatrix          Cov;            // covariances

    // derivatives
    CFortranMatrix          Kder;           // derivative of kernels w.r.t. a hyperparameter

    bool                    GlobalMinSet;   // true if gpos set by SetGlobalMin()
    CSimpleVector<double>   GPos;           // global position, either detected or use
    bool                    GPosSet;        // true is gpos set by any means, either SetGlobalMin() or from FES
    int                     GPosBin;

    // SVD setup
    double                  RCond;

    // internal flow
    bool                    UseZeroPoint;
    bool                    UseInv;         // calc all via inversion
    bool                    NeedInv;        // need inverted matrix - hyprms, error analysis
    bool                    FastErrors;     // use faster but more memory intensive algorithm

    // statistical reweighting
    std::vector<bool>       IntegratedCVs;

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);
    void CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout); // gpos - position of global minimum
    void CalculateCovs(CVerboseStr& vout);
    void CalculateErrorsFromCov(CVerboseStr& vout);

    double GetValue(const CSimpleVector<double>& position);
    double GetVar(CSimpleVector<double>& lpos);
    // optimized version var+cov
    void   GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& llvar,double& lrcov);

    // kernel matrix + noise
    void   CreateKS(void);
    void   CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& ky);
    void   CreateKff2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& ky2);
    double GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp);
    void   GetKernelDerAnaI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& der);
    void   GetKernelDerNumI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& der);
    void   GetKernelDerAnaJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& der);
    void   GetKernelDerNumJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& der);
    void   GetKernelDer2Ana(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);
    void   GetKernelDer2Num(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);
    void   GetKernelDer2AnaWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);
    void   GetKernelDer2NumWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);


    // derivatives
    void CalcKderWRTSigmaF2(void);
    void CalcKderWRTNCorr(void);
    void CalcKderWRTWFac(size_t cv);
    void CalcKderWRTSigmaN2(size_t cv);

    // parallel processing
    void RunBlasLapackSeq(void);
    void RunBlasLapackPar(void);
};

//------------------------------------------------------------------------------

#endif
