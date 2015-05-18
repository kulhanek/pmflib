#ifndef MTDServerH
#define MTDServerH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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
// ===============================================================================

#include <ExtraServer.hpp>
#include <PrmFile.hpp>
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <RegClientList.hpp>
#include "MTDSrvOptions.hpp"
#include "MTDServerHist.hpp"

//------------------------------------------------------------------------------

class CMTDServer : public CExtraServer {
public:
    // constructor
    CMTDServer(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data ----------------------------------------------------
private:
    CMTDSrvOptions      Options;        // program options
    CPrmFile            Controls;       // controls
    CSmallString        InputFileName;  // input file name
    CSmallString        OutputFileName; // output file name

    // global data -------------------------------
    CMTDServerHist      MTDHistory;         // global MTD history

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// Ctrl+C signal handler
    static void CtrlCSignalHandler(int signal);

    /// process control file
    bool ProcessFilesControl(void);

    friend class CMTDProcessor;
    friend class CMTDRegClient;
};

//------------------------------------------------------------------------------

extern CMTDServer MTDServer;

//------------------------------------------------------------------------------

#endif
