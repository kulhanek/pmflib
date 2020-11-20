#ifndef CABFDerivativesH
#define CABFDerivativesH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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

#include "ABFDeriOptions.hpp"
#include <ABFAccumulator.hpp>
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>

//------------------------------------------------------------------------------

/// utility to extract free energy derivatives from ABF accumulator

class CABFDerivatives {
public:
    CABFDerivatives(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CABFDeriOptions     Options;
    CStdIOFile          InputFile;
    CStdIOFile          OutputFile;
    CABFAccumulator     Accumulator;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// print derivatives
    bool PrintDG(void);     // mean force
};

//------------------------------------------------------------------------------

#endif
