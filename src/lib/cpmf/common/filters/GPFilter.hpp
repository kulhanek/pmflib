#ifndef GPFilterH
#define GPFilterH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <Filter.hpp>
#include <FortranMatrix.hpp>
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>

//------------------------------------------------------------------------------

// supported kernel functions
// ARD = automatic relevance determination

enum EGPFilterKernel {
    EGPFK_ARDSE   = 1,    // ARD squared exponential (radial basis function)
    EGPFK_ARDMC52 = 2,    // ARD Matern class 5/2 function
    EGPFK_ARDMC32 = 3,    // ARD Matern class 3/2 function
};

//------------------------------------------------------------------------------

/// Gaussian Process filter

class PMF_PACKAGE CGPFilter : public CFilter {
public:
// constructor and destructor -------------------------------------------------
    CGPFilter(void);
    virtual ~CGPFilter(void);

// executive methods ----------------------------------------------------------
    /// set kernel
    void SetKernel(const CSmallString& kernel);

    /// get kernel name
    const CSmallString GetKernelName(void);

    /// setup Gaussian process
    void SetFilter(double timestep,int framelen, double wfac, double noise,CVerboseStr& vout);

    /// train process
    void TrainProcess(const CVectorDataPtr& in, size_t offset);

    /// predict data
    void PredictData(CVectorDataPtr& out, size_t offset);

    /// get logML
    double GetLogML(void);

    /// get logPL
    double GetLogPL(void);

    /// filter data
    virtual void RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out);

// section of private data ----------------------------------------------------
protected:
    int                     NumOfThreads;
    size_t                  GPRSize;
    double                  SigmaF2;
    double                  WFac;
    double                  Noise;
    double                  RCond;
    EGPFilterKernel         Kernel;

    CFortranMatrix          KS;              // kernel matrix
    double                  logdetK;
    double                  Offset;
    CSimpleVector<double>   Y;              // mean forces
    CSimpleVector<double>   GPRModel;       // weights
    CSimpleVector<double>   KFF;

    // parallel processing
    void PrintExecInfo(CVerboseStr& vout);
    void RunBlasLapackSeq(void);
    void RunBlasLapackPar(void);

    // kernel values
    double GetKernelValue(double indi,double indj);

    // update KFF
    void GetKFF(CSimpleVector<double>& kff, double indi);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CGPFilter>    CGPFilterPtr;

//------------------------------------------------------------------------------

#endif
