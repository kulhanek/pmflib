#ifndef MTDEnergyH
#define MTDEnergyH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
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

#include <stdio.h>
#include <PMFAccumulator.hpp>
#include <SimpleVector.hpp>
#include <EnergySurface.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>
#include <EnergyProxy.hpp>
#include "MTDEneOptions.hpp"
#include <list>

//------------------------------------------------------------------------------

class CMTDEnergy {
public:
    // constructor
    CMTDEnergy(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CMTDEneOptions                  Options;            // program options
    CStdIOFile                      InputFile;
    CStdIOFile                      OutputFile;
    std::list<CPMFAccumulatorPtr>   Accus;
    CEnergySurfacePtr               FES;                // energy surface
    CTerminalStr                    Console;
    CVerboseStr                     vout;
};

//------------------------------------------------------------------------------

#endif
