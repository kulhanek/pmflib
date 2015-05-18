#ifndef RSTSamplesH
#define RSTSamplesH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "RSTSamplesOptions.hpp"
#include <RSTAccumulator.hpp>
#include <SimpleVector.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>

//------------------------------------------------------------------------------

class CRSTSamples {
public:
    // constructor
    CRSTSamples(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CRSTSamplesOptions  Options;            // program options
    FILE*               InputFile;          // input file
    bool                OwnInputFile;       // do we own input file handle?
    FILE*               OutputFile;         // output file
    bool                OwnOutputFile;      // do we own output file handle?
    CRSTAccumulator     Accumulator;        // RST accumulator
    CSimpleVector<int>  Point;              // current point in accumulator

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    // print samples - it is called recursivelly
    bool PrintSamples(int cv);
};

//------------------------------------------------------------------------------

#endif
