#ifndef MWAServerH
#define MWAServerH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <ExtraServer.hpp>
#include <PrmFile.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <SimpleMutex.hpp>
#include "MWASrvOptions.hpp"
#include "MWAServerAccu.hpp"

//------------------------------------------------------------------------------

class CMWARegClient;

//------------------------------------------------------------------------------

class CMWAServer : public CExtraServer {
public:
    // constructor
    CMWAServer(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CMWASrvOptions      Options;        // program options
    CPrmFile            Controls;       // controls
    CSmallString        InputFileName;  // input file name
    CSmallString        OutputFileName; // output file name
    CSmallString        TrajFileName;   // output file name
    FILE*               TrajFile;
    bool                AutoRestart;
    int                 AutoSaveInterval;
    int                 SaveCounter;
    CSimpleMutex        AutoSaveMutex;

    // global data -------------------------------
    CMWAServerAccu      MWAAccumulator;     // global MWA accumulator

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// Ctrl+C signal handler
    static void CtrlCSignalHandler(int signal);

    /// process control file
    void ProcessFilesControl(void);

    /// process data for synchronous mode
    void ProcessSyncControl(void);

    /// autosave MWA accumulator
    void AutoSaveData(void);

    friend class CMWAProcessor;
    friend class CMWAServerAccu;
};

//------------------------------------------------------------------------------

extern CMWAServer MWAServer;

//------------------------------------------------------------------------------

#endif
