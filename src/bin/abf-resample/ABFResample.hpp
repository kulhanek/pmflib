#ifndef ABFResampleH
#define ABFResampleH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2022 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Martin Petrek, petrek@chemi.muni.cz
//                       Petr Kulhanek, kulhanek@enzim.hu
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

#include "ABFResampleOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <SmallTimeAndDate.hpp>
#include <PMFAccumulator.hpp>
#include <EnergyDerProxy.hpp>
#include <EnergySurface.hpp>
#include <IntegratorGPR.hpp>

//------------------------------------------------------------------------------

class CIntegratorGPR;
class CSmootherGPR;

//------------------------------------------------------------------------------

/// utility to integrate ABF accumulator

class CABFResample {
public:
    CABFResample(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CABFResampleOptions Options;
    CPMFAccumulatorPtr  InAccu;
    std::vector<int>    NBins;
    CPMFAccumulatorPtr  OutAccu;
    CEnergyDerProxyPtr  DerProxy;
    CEnergyDerProxyPtr  DerProxt;
    CEnergySurfacePtr   FES;
    CIntegratorGPR      Integrator;
    CPMFAccuDataPtr     NSAMPLES;
    CPMFAccuDataPtr     MICF;
    CPMFAccuDataPtr     M2ICF;
    CSmallTimeAndDate   StartTime;
    int                 State;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    void PrepResampledAccu(void);
    void ResampleAccu(void);
    void PrintSampledStat(void);
    bool Integrate(void);
    void LoadGPRHyprms(CIntegratorGPR& gpr);
    void DecodeIList(const CSmallString& spec, std::vector<int>& ilist,const CSmallString& optionname);
};

//------------------------------------------------------------------------------

#endif
