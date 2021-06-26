#ifndef PMFAccuInfoH
#define PMFAccuInfoH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "PMFAccuInfoOptions.hpp"
#include <PMFAccumulator.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <StdIOFile.hpp>

//------------------------------------------------------------------------------

/// utility to integrate ABF accumulator

class CPMFAccuInfo {
public:
    CPMFAccuInfo(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CPMFAccuInfoOptions Options;
    CStdIOFile          InputFile;
    CPMFAccumulatorPtr  Accu;

    // output ------------------------------------
    CTerminalStr    Console;
    CVerboseStr     vout;

    // data to print -----------------------------
    CSimpleVector<double>   Values;
    CSimpleVector<double>   Sigmas;
    CSimpleVector<double>   Errors;

// actions
    // print overall info
    void PrintInfo(void);

    // list data sections
    void ListSections(void);

    // get data from the section
    void GetSection(const CSmallString& name);

    // get MICF from ABF
    void GetMICF(void);

    // print header
    void PrintHeader(const CSmallString& sec_name, bool print_errors);

    // print data
    void PrintData(bool print_errors);

    // print data per CV
    void PrintDataPerCV(const CSmallString& sec_name);
};

//------------------------------------------------------------------------------

#endif
