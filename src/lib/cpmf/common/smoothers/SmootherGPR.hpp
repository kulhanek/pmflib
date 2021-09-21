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

class PMF_PACKAGE CSmootherGPR {
public:
// constructor and destructor -------------------------------------------------
    CSmootherGPR(void);
    virtual ~CSmootherGPR(void);

// setup methods --------------------------------------------------------------
    /// set input energy proxy
    void AddInputEnergyProxy(CEnergyProxyPtr p_proxy);

    /// clear energy proxies
    void ClearInputEnergyProxies(void);

    /// set output energy surface
    void SetOutputES(CEnergySurfacePtr p_surf);

// hyperparameters
    /// set sigmaf2
    void SetSigmaF2(double sigf2);

    /// set ncorr
    void SetNCorr(const CSmallString& spec);

    /// set ncorr
    void SetNCorr(double value);

    /// multiply of bin sizes
    void SetWFac(const CSmallString& spec);

    /// multiply of bin sizes
    void SetWFac(CSimpleVector<double>& wfac);

    /// multiply of bin sizes
    void SetWFac(size_t cvind, double value);

// setup
    /// set include error
    void SetIncludeError(bool set);

    /// set algorithm for LA
    void SetLAMethod(EGPRLAMethod set);

    /// set algorithm for LA
    void SetLAMethod(const CSmallString& method);

    /// set kernel
    void SetKernel(const CSmallString& kernel);

    /// set rcond for SVD
    void SetRCond(double rcond);

    /// set position of global minimum
    void SetGlobalMin(const CSmallString& spec);

    /// set position of global minimum
    void SetGlobalMin(const CSimpleVector<double>& pos);

    /// get position of global minima
    CSimpleVector<double> GetGlobalMin(void);

    /// get position of global minima - bin
    int GetGlobalMinBin(void);

    /// get position of global minima
    double GetGlobalMinValue(void) const;

    /// use inversion alg
    void SetUseInv(bool set);

// execution method -----------------------------------------------------------
    /// interpolate data
    bool Interpolate(CVerboseStr& vout,bool nostat=false);

    /// print exec info
    void PrintExecInfo(CVerboseStr& vout);

    /// get kernel name
    const CSmallString GetKernelName(void);

// GPR model optimization ----------------------------------------------------
    /// calc hyprms grd
    void PrepForHyprmsGrd(bool set);

    /// calc logpl
    void SetCalcLogPL(bool set);

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

// section of private data ----------------------------------------------------
private:
    std::vector<CEnergyProxyPtr>    EneProxyItems;
    CEnergySurfacePtr               EneSurface;
    int                             NumOfThreads;

    // GPR data, sizes and index maps
    size_t                          NumOfCVs;
    size_t                          NumOfBins;
    size_t                          GPRSize;
    std::vector<size_t>             SampledMap;
    std::vector<size_t>             EneProxyMap;
    size_t                          NumOfValues;
    std::vector<size_t>             ValueMap;

    // setup
    bool                    IncludeError;
    EGPRLAMethod            Method;

    // hyperparameters
    double                  SigmaF2;
    double                  NCorr;
    CSimpleVector<double>   WFac;
    CSimpleVector<double>   CVLengths2;

    // GPR model
    EGPRKernel              Kernel;
    CFortranMatrix          KS;             // kernel matrix
    double                  logdetK;
    double                  Mean;
    CSimpleVector<double>   Y;              // enthalpy
    CSimpleVector<double>   GPRModel;       // weights
    CFortranMatrix          Cov;            // covariances

    bool                    GlobalMinSet;   // true if gpos set by SetGlobalMin()
    CSimpleVector<double>   GPos;           // global position, either detected or use
    bool                    GPosSet;        // true is gpos set by any means, either SetGlobalMin() or from FES
    double                  GlbMinValue;
    int                     GPosBin;

    // derivatives
    CFortranMatrix          Kder;           // derivative of kernels w.r.t. a hyperparameter

    // internal
    bool                    UseInv;         // calc all via inversion
    bool                    NeedInv;        // need inverted matrix - hyprms, error analysis

    // SVD setup
    double                  RCond;

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);
    void CalculateCovs(CVerboseStr& vout);
    void CalculateErrorsFromCov(CVerboseStr& vout);

    double GetValue(const CSimpleVector<double>& position);

    // kernel matrix + noise
    void   CreateKS(void);
    void   CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& ky);
    double GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp);

    // parallel processing
    void RunBlasLapackSeq(void);
    void RunBlasLapackPar(void);

    // derivatives
    void CalcKderWRTSigmaF2(void);
    void CalcKderWRTNCorr(void);
    void CalcKderWRTWFac(size_t cv);
    double GetKernelValueWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv);
};

//------------------------------------------------------------------------------

#endif
