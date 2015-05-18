#ifndef BMIntegrateH
#define BMIntegrateH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
//                       Martin Petrek, petrek@chemi.muni.cz
//    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <DynamicArray.hpp>
#include "BMIntOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>

//------------------------------------------------------------------------------

class CBMIntegrate {
public:
    // constructor
    CBMIntegrate(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CBMIntOptions   Options;            // program options
    FILE*           InputFile;          // input file
    bool            OwnInputFile;       // do we own input file handle?
    FILE*           OutputFile;         // output file
    bool            OwnOutputFile;      // do we own output file handle?
    CSmallString    OutputFormat;       // output format
    CDynamicArray   X,Y,S,I,IE;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// read data from input stream
    bool ReadData(void);

    /// integrate data
    bool IntegrateData(void);

    /// print integrated data to output stream
    bool PrintData(void);
};

//------------------------------------------------------------------------------

#endif
