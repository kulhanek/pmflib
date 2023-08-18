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
#include <IntegratorGPR.hpp>

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CGHSIntegratorGPR0A {
public:
// constructor and destructor -------------------------------------------------
    CGHSIntegratorGPR0A(void);
    virtual ~CGHSIntegratorGPR0A(void);

// setup methods --------------------------------------------------------------
    /// set input derivative proxies
    void SetGDerProxy(CEnergyDerProxyPtr p_proxy);
    void SetHDerProxy(CEnergyDerProxyPtr p_proxy);
    void SetSDerProxy(CEnergyDerProxyPtr p_proxy);

    /// set output energy surfaces
    void SetOutputFES(CEnergySurfacePtr p_surf);
    void SetOutputHES(CEnergySurfacePtr p_surf);
    void SetOutputSES(CEnergySurfacePtr p_surf);

// hyperparameters
    /// set sigmaf2
    void SetSigmaF2(const CSmallString& spec);

    /// set sigmaf2
    void SetSigmaF2(CSimpleVector<double>& sigmaf2);

    /// set sigmaf2
    void SetSigmaF2(size_t cvind, double value);

    /// set covariances
    void SetCoVar(const CSmallString& spec);

    /// set covariances
    void SetCoVar(CSimpleVector<double>& covar);

    /// set covariances
    void SetCoVar(size_t cvind, double value);

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
    /// enable constraints
    void EnableConstraints(bool set);

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

    /// use fast error algorithm
    void SetFastError(bool set);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout,bool nostat=false);

    /// get log of Marginal Likelihood
    double GetLogML(void);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    double GetLogPL(void);

    /// print exec info
    void PrintExecInfo(CVerboseStr& vout);

    /// get kernel name
    const CSmallString GetKernelName(void);

// section of private data ----------------------------------------------------
private:
    CEnergyDerProxyPtr      GDerProxy;
    CEnergyDerProxyPtr      HDerProxy;
    CEnergyDerProxyPtr      SDerProxy;
    CEnergySurfacePtr       GSurface;
    CEnergySurfacePtr       HSurface;
    CEnergySurfacePtr       SSurface;
    int                     NumOfThreads;

    // GPR data, sizes and index maps
    size_t                  NumOfCVs;
    size_t                  NumOfBins;
    size_t                  GPRSize;
    size_t                  NumOfUsedBins;
    std::vector<size_t>     SampledMap;

    // setup
    bool                    ConstrainedTK;
    bool                    UseNumDiff;
    bool                    NoEnergy;
    bool                    IncludeError;
    EGPRLAMethod            Method;

    // hyperparameters
    CSimpleVector<double>   SigmaF2;
    CSimpleVector<double>   CoVar;
    CSimpleVector<double>   SigmaN2;
    CSimpleVector<double>   WFac;
    CSimpleVector<double>   CVLengths2;

    // GPR model
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

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);

    void GetValues(const CSimpleVector<double>& position,double& dg, double& dh, double& mtds);

    void   CreateTK(void);
    void   CreateKS(void);
    void   CreateKff(const CSimpleVector<double>& ip,CSimpleVector<double>& ky,int task);
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

    // parallel processing
    void RunBlasLapackSeq(void);
    void RunBlasLapackPar(void);
};

//------------------------------------------------------------------------------

#endif
