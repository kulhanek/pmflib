#ifndef StringServerH
#define StringServerH
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

#include <ExtraServer.hpp>
#include <PrmFile.hpp>
#include <RegClientList.hpp>
#include <SimpleCond.hpp>
#include "StringSrvOptions.hpp"
#include "BeadList.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include "Launcher.hpp"

//------------------------------------------------------------------------------

class CStringRegClient;

//------------------------------------------------------------------------------

class CStringServer : public CExtraServer {
public:
    // constructor
    CStringServer(void);

// main methods ---------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

    /// add data to global accumulator
    bool AddDataToStringAccu(CStringRegClient* p_client);

    /// terminate server - it redefines CServer::TerminateServer()
    void TerminateServer(void);

// section of private data ----------------------------------------------------
private:
    CStringSrvOptions   Options;        // program options
    CPrmFile            Controls;       // controls

    // control data ------------------------------
    CSmallString        ServerKeyName;
    int                 Port;

    // global data -------------------------------
    CRegClientList      RegClients;     // registered clients
    CBeadList           Beads;          // available beads
    CLauncher           Launcher;       // job luncher

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// Ctrl+C signal handler
    static void CtrlCSignalHandler(int signal);

    friend class CStringProcessor;
    friend class CLauncher;
};

//------------------------------------------------------------------------------

extern CStringServer StringServer;

//------------------------------------------------------------------------------

#endif
