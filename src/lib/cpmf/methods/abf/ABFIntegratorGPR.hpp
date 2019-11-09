#ifndef ABFIntegratorGPRH
#define ABFIntegratorGPRH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

//------------------------------------------------------------------------------

class CABFAccumulator;
class CEnergySurface;

// how to invert the matrix

enum EGPRINVMethod {
    EGPRINV_LU      = 1,    // DGETRI/DGETRF (via LU factorization)
    EGPRINV_SVD     = 2,    // via SVD factorization, divide and conquer driver
    EGPRINV_SVD2    = 3,    // via SVD factorization, simple driver
    EGPRINV_LL      = 4,    // DPOTRF/DPOTRI (via LL factorization)
};

//------------------------------------------------------------------------------

// supported kernel functions
// ARD = automatic relevance determination

enum EGPRKernel {
    EGPRK_ARDSE   = 1,    // ARD squared exponential (radial basis function)
    EGPRK_ARDMC52 = 3,    // ARD Matern class 5/2 function
};

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CABFIntegratorGPR {
public:
// constructor and destructor -------------------------------------------------
    CABFIntegratorGPR(void);
    virtual ~CABFIntegratorGPR(void);

// setup methods --------------------------------------------------------------
    /// set input ABF accumulator, only ABF forces are integrated
    void SetInputABFAccumulator(CABFAccumulator* p_accu);

    /// set output free energy surface
    void SetOutputFESurface(CEnergySurface* p_surf);

// hyperparameters
    /// set sigmaf2
    void SetSigmaF2(double sigf2);

    /// set ncorr
    void SetNCorr(const CSmallString& spec);

    /// set ncorr
    void SetNCorr(CSimpleVector<double>& wfac);

    /// set ncorr
    void SetNCorr(size_t cvind, double value);

    /// multiply of bin sizes
    void SetWFac(const CSmallString& spec);

    /// multiply of bin sizes
    void SetWFac(CSimpleVector<double>& wfac);

    /// multiply of bin sizes
    void SetWFac(size_t cvind, double value);

// setup
    /// set include error
    void SetIncludeError(bool set);

    /// skip energy calculation, it also disables errors
    void SetNoEnergy(bool set);

    /// switch to numerical evaluation of kernel fce derivatives
    void SetUseNumDiff(bool set);

    /// use split ncorr
    void SetSplitNCorr(bool set);

    /// set algorithm for LLS
    void SetINVMethod(EGPRINVMethod set);

    /// set algorithm for LLS
    void SetINVMethod(const CSmallString& method);

    /// set kernel
    void SetKernel(const CSmallString& kernel);

    /// set rcond for SVD
    void SetRCond(double rcond);

    /// should we include glued area to energy calculation?
    void IncludeGluedAreas(bool set);

    /// set position of global minimum
    void SetGlobalMin(const CSmallString& spec);

// execution method -----------------------------------------------------------
    /// integrate data
    bool Integrate(CVerboseStr& vout,bool nostat=false);

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

    /// write file with derivatives
    bool WriteMFInfo(const CSmallString& name);

    /// remove mean force outliers from ABF data
    void FilterByMFZScore(double maxzscore,CVerboseStr& vout);

    /// print exec info
    void PrintExecInfo(CVerboseStr& vout);

    /// get kernel name
    const CSmallString GetKernelName(void);

// section of private data ----------------------------------------------------
private:
    CABFAccumulator*        Accumulator;
    CEnergySurface*         FES;

    int                     NumOfThreads;

    // GPR data, sizes and index maps
    size_t                  NCVs;
    size_t                  NumOfBins;
    size_t                  GPRSize;
    size_t                  NumOfUsedBins;
    CSimpleVector<size_t>   SampledMap;
    size_t                  NumOfValues;
    CSimpleVector<size_t>   ValueMap;

    // setup
    bool                    UseNumDiff;
    bool                    NoEnergy;
    bool                    IncludeError;
    bool                    IncludeGluedBins;
    EGPRINVMethod           Method;

    // hyperparameters
    bool                    SplitNCorr;     // ncorr for each cv / all cvs
    double                  SigmaF2;
    CSimpleVector<double>   NCorr;
    CSimpleVector<double>   WFac;
    CSimpleVector<double>   CVLengths2;

    // GPR model
    EGPRKernel              Kernel;
    CFortranMatrix          KS;              // kernel matrix
    double                  logdetK;
    CSimpleVector<double>   Y;              // mean forces
    CSimpleVector<double>   GPRModel;       // weights

    // derivatives
    CFortranMatrix          Kder;           // derivative of kernels w.r.t. a hyperparameter

    bool                    GlobalMinSet;
    CSimpleVector<double>   GPos;           // global position, either detected or use

    // SVD setup
    double                  RCond;

    bool TrainGP(CVerboseStr& vout);
    void CalculateEnergy(CVerboseStr& vout);
    void CalculateErrors(CSimpleVector<double>& gpos,CVerboseStr& vout); // gpos - position of global minimum

    double GetValue(const CSimpleVector<double>& position);
    double GetMeanForce(const CSimpleVector<double>& position,size_t icoord);
    double GetVar(CSimpleVector<double>& lpos);
    // optimized version var+cov
    void   GetCovVar(CSimpleVector<double>& lpos,CSimpleVector<double>& rpos,double& llvar,double& lrcov);

    // kernel matrix + noise
    void   CreateKS(void);
    void   CreateKy(const CSimpleVector<double>& ip,CSimpleVector<double>& ky);
    void   CreateKy2(const CSimpleVector<double>& ip,size_t icoord,CSimpleVector<double>& ky2);
    double GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp);
    void   GetKernelDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& der);
    void   GetKernelDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& der);
    void   GetKernelDer2Ana(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);
    void   GetKernelDer2Num(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);
    void   GetKernelDer2AnaWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);
    void   GetKernelDer2NumWFac(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);


    // derivatives
    void CalcKderWRTSigmaF2(void);
    void CalcKderWRTNCorr(size_t cv);
    void CalcKderWRTWFac(size_t cv);

    // parallel processing
    void RunBlasLapackSeq(void);
    void RunBlasLapackPar(void);
};

//------------------------------------------------------------------------------

#endif
