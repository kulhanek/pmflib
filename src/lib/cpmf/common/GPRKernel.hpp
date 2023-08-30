#ifndef GPRKernelH
#define GPRKernelH
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
#include <FortranMatrix.hpp>
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

// supported kernel functions
// ARD = automatic relevance determination

enum EGPRKernel {
    EGPRK_NONE      = 0,
    EGPRK_ARDSE     = 1,    // ARD squared exponential (radial basis function)
    EGPRK_ARDMC52   = 2,    // ARD Matern class 5/2 function
    EGPRK_ARDMC32   = 3,    // ARD Matern class 3/2 function
    EGPRK_ARDMC12   = 4,    // ARD Matern class 1/2 function
};

//------------------------------------------------------------------------------

/** \brief GPRKernel
*/

class PMF_PACKAGE CGPRKernel {
public:
// constructor and destructor -------------------------------------------------
    CGPRKernel(void);
    virtual ~CGPRKernel(void);

    /// set accumulator
    void SetAccumulator(CPMFAccumulatorPtr accu);

    /// switch to numerical evaluation of kernel fce derivatives
    void SetUseNumDiff(bool iset);

    /// is numeric differentiation enabled
    bool IsNumDiffEnabled(void);

// ----
    /// set kernel by name
    void SetKernel(const CSmallString& kernel);

    /// get kernel name
    const CSmallString GetKernelName(void);

// ----
    /// multiply of bin sizes
    void SetWFac(const CSmallString& spec);

    /// multiply of bin sizes
    void SetWFac(CSimpleVector<double>& wfac);

    /// multiply of bin sizes
    void SetWFac(size_t cvind, double value);

// ----
    /// setup kernel for calculation
    void SetupKernel(void);

// ----
    /// get kernel value
    double GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp);

    // kernel derivatives
    void   GetKernelDerI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder);
    void   GetKernelDerJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder);
    void   GetKernelDerIJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);

    // derivatives by wfac for hyper-parameter gradients
    double GetKernelValueWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv);
    void   GetKernelDerIWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder);
    void   GetKernelDerJWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder);
    void   GetKernelDerIJWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);

// section of protected data ---------------------------------------------------
protected:
    CPMFAccumulatorPtr      Accu;
    size_t                  NumOfCVs;
    size_t                  NumOfBins;
    EGPRKernel              Kernel;
    CSimpleVector<double>   WFac;

// section of private data ----------------------------------------------------
private:
    CSimpleVector<double>   CVLengths2;
    bool                    UseNumDiff;

// analytical derivatives
    void   GetKernelDerIAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder);
    void   GetKernelDerJAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder);
    void   GetKernelDerIJAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);

// numerical derivatives
    void   GetKernelDerINum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder);
    void   GetKernelDerJNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder);
    void   GetKernelDerIJNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock);

// derivatives by wfac fro hyperparameter gradients
    double GetKernelValueWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv);
    void   GetKernelDerIWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder);
    void   GetKernelDerJWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder);
    void   GetKernelDerIJWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);

    double GetKernelValueWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv);
    void   GetKernelDerIWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder);
    void   GetKernelDerJWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder);
    void   GetKernelDerIJWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock);
};

//------------------------------------------------------------------------------

#endif
