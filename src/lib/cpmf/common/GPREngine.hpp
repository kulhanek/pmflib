#ifndef GPRHyprmsH
#define GPRHyprmsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <vector>
#include <SimpleVector.hpp>
#include <boost/shared_ptr.hpp>
#include <GPRKernel.hpp>
#include <VerboseStr.hpp>

//------------------------------------------------------------------------------

// linear algebra pathway

enum EGPRLAMethod {
    EGPRLA_LU      = 1,    // DGETRI/DGETRF (via LU factorization)
    EGPRLA_SVD     = 2,    // via SVD factorization, divide and conquer driver
    EGPRLA_SVD2    = 3,    // via SVD factorization, simple driver
    EGPRLA_LL      = 4,    // DPOTRF/DPOTRI (via LL factorization)
};

//------------------------------------------------------------------------------

/** \brief base for GPR
*/

class PMF_PACKAGE CGPREngine: public CGPRKernel {
public:
// constructor and destructor -------------------------------------------------
    CGPREngine(void);
    virtual ~CGPREngine(void);

// hyperparameters
    /// set sigmaf2
    void SetSigmaF2(const CSmallString& spec);

    /// set sigmaf2
    void SetSigmaF2(CSimpleVector<double>& sigmaf2);

    /// set sigmaf2
    void SetSigmaF2(size_t cvind, double value);

//-----
    /// set covariances
    void SetCoVar(const CSmallString& spec);

    /// set covariances
    void SetCoVar(CSimpleVector<double>& covar);

    /// set covariances
    void SetCoVar(size_t cvind, double value);

//-----
    // WFac is part of CGPRKernel

//-----
    /// set sigman2
    void SetSigmaN2(const CSmallString& spec);

    /// set sigman2
    void SetSigmaN2(CSimpleVector<double>& sigman2);

    /// set sigman2
    void SetSigmaN2(size_t cvind, double value);

//-----
    /// load hyper-parameters from file
    void LoadGPRHyprms(const CSmallString& name);

// parallel processing ---------------------------------------------------------
    /// switch to sequential mode
    void RunBlasLapackSeq(void);

    /// switch to parallel mode
    void RunBlasLapackPar(void);

    /// print exec info
    void PrintExecInfo(CVerboseStr& vout);

// linear algebra setup --------------------------------------------------------
    /// set algorithm for LA
    virtual void SetLAMethod(EGPRLAMethod set);

    /// set algorithm for LA
    virtual void SetLAMethod(const CSmallString& method);

    /// set rcond for SVD
    virtual void SetRCond(double rcond);

    /// use inversion alg
    virtual void SetUseInv(bool iset);

// get parameters --------------------------------------------------------------

    virtual int GetNumOfSigmaF2(void);
    virtual int GetNumOfCoVar(void);
    virtual int GetNumOfWFac(void);
    virtual int GetNumOfSigmaN2(void);
    virtual int GetNumOfHyprms(void);

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

    /// use the first derivatives of the kernel
    virtual void UseFirstKernelDerivatives(bool set);

// base methods ----------------------------------------------------------------
    /// run GPR
    virtual bool RunGPR(CVerboseStr& vout,bool nostat=false);

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
// hyper-parameters
    size_t                  NumOfSigmaF2;
    CSimpleVector<double>   SigmaF2;

    size_t                  NumOfCoVar;
    CSimpleVector<double>   CoVar;

    // CSimpleVector<double>   WFac;    // part of CGPRkernel

    size_t                  NumOfSigmaN2;
    CSimpleVector<double>   SigmaN2;

// linear algebra setup
    EGPRLAMethod            Method;
    double                  RCond;      // SVD setup
    bool                    UseInv;     // calc all via inversion
    bool                    NeedInv;    // need inverted matrix - hyprms, error analysis

// setup
    bool                    NoEnergy;
    bool                    IncludeError;

// section of private data ----------------------------------------------------
private:
    int                     NumOfThreads;
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CGPREngine>    CGPREnginePtr;

//------------------------------------------------------------------------------

#endif
