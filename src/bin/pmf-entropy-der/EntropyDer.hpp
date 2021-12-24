#ifndef EntropyDerH
#define EntropyDerH
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

#include "EntropyDerOptions.hpp"
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <PMFAccumulator.hpp>

enum EKSKernelType {
    EKSKT_PARABOLIC = 1,
    EKSKT_TRIWEIGHT = 2,
    EKSKT_GAUSSIAN  = 3,
};

//------------------------------------------------------------------------------

/// utility to extract enthalpy from  accumulator

class CEntropyDer {
public:
    CEntropyDer(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CEntropyDerOptions              Options;
    CStdIOFile                      InputFile;
    CStdIOFile                      OutputFile;
    CPMFAccumulatorPtr              InAccu;
    CPMFAccumulatorPtr              OutAccu;
    int                             State;

    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    size_t                  NSTLimit;
    size_t                  NumOfBins;
    size_t                  NumOfCVs;
    EKSKernelType           KSKernel;
    CSimpleVector<double>   KSWFac;

    void GetEtot(void);
    void GetICFBySGF(void);
    void GetICFByGPF(double wfac);

    void CalculatePPandPN(void);
    void CalculatePPandPN_KS(void);

    void CalculateKSWeights(CSimpleVector<double>& weights,CSimpleVector<double>& jpos);
    double GetKSKernelValue(double u2);

    void SetKSWFac(const CSmallString& spec);
    void SetKSKernel(const CSmallString& kskernel);
};

//------------------------------------------------------------------------------

#endif
