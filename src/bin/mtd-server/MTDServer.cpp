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

#include <stdio.h>
#include <ErrorSystem.hpp>
#include <signal.h>
#include <SmallTime.hpp>
#include <CmdProcessorList.hpp>
#include <PMFOperation.hpp>
#include <ExtraOperation.hpp>
#include "MTDServer.hpp"
#include "MTDProcessor.hpp"
#include "MTDFactory.hpp"

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

CMTDServer MTDServer;
MAIN_ENTRY_OBJECT(MTDServer)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDServer::CMTDServer(void)
{
    SetProtocolName("mtd");
    OutputFileName = "_mtdserver.rst";
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDServer::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# mtd-server (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgControlFile() != "-") {
        vout << "# Control file name : " << Options.GetArgControlFile() << endl;
    } else {
        vout << "# Control file name : - (standard input)" << endl;
    }

// set SIGINT/SIGTERM hadler to cleanly shutdown server ----------
    signal(SIGINT,CtrlCSignalHandler);
    signal(SIGTERM,CtrlCSignalHandler);

// process control file ----------------------------------
    if(Options.GetArgControlFile() != "-") {
        if(Controls.Read(Options.GetArgControlFile()) == false) {
            ES_ERROR("unable to open and read control file");
            return(SO_USER_ERROR);
        }
    } else {
        if(Controls.Parse(stdin) == false) {
            ES_ERROR("unable to read control file from standard input");
            return(SO_USER_ERROR);
        }
    }

    if( ProcessServerControl(Controls,vout) == false ){
       return( SO_USER_ERROR );
       }

    if( ProcessFilesControl() == false ) {
        return(SO_USER_ERROR);
    }

    if( Controls.CountULines() > 0 ){
       ES_ERROR("unprocessed items found in control file");
       Controls.Dump(stderr,true);
       return( SO_USER_ERROR );
       }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CMTDServer::ProcessFilesControl(void)
{
    vout << endl;
    vout << "=== [files] ====================================================================" << endl;
    if(Controls.OpenSection("files") == false) {
        vout << "Input history potential (input)   = none                               (default)" << endl;
        vout << "Output history potential (output) = " << OutputFileName << "                     (default)" << endl;
        return(true);
    }

    if(Controls.GetStringByKey("input",InputFileName) == true) {
        vout << "Input history potential (input)   = " << InputFileName << endl;
    } else {
        vout << "Input history potential (input)   = none                               (default)" << endl;
    }

    if(Controls.GetStringByKey("output",OutputFileName) == true) {
        vout << "Output history potential (output) = " << OutputFileName << endl;
    } else {
        vout << "Output history potential (output) = " << OutputFileName << "                     (default)" << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDServer::Run(void)
{
    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Input data :::::::::::::::::::::::::::::::::" << endl;

    if(InputFileName != NULL) {
        vout << "Input MTD history file: " << InputFileName << endl;
        try {
            MTDHistory.Load(InputFileName);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to load input MTD history file",e);
            return(false);
        }
        // this is very important operation
        // data exchange between client and server is quantized, quantum is whole buffer
        // so no empty gap has to be presented in it
        MTDHistory.ReallocateAsSingleBuffer();
    } else {
        vout << "Input MTD history file: none" << endl;
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: MTD Server :::::::::::::::::::::::::::::::::" << endl;

// register operations
    CmdProcessorList.RegisterProcessor(Operation_GetServerInfo,&MTDFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_FlushServerData,&MTDFactory);
    CmdProcessorList.RegisterProcessor(Operation_RegisterClient,&MTDFactory);
    CmdProcessorList.RegisterProcessor(Operation_UnregisterClient,&MTDFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_GetInitialData,&MTDFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_ExchangeData,&MTDFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_GetData,&MTDFactory);

// start server listening
    if( StartServer() == false ){
        return(false);
    }

    vout << "Number of processed transactions: " << GetNumberOfTransactions() << endl;
    vout << "Number of illegal transactions  : " << GetNumberOfIllegalTransactions() << endl;

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;

// and now save all data
    if(MTDHistory.GetNumberOfCoords() > 0) {
        vout << "Output MTD history file: " << OutputFileName << endl;
        try{
            MTDHistory.Save(OutputFileName);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to save final MTD history file",e)
            return(false);
        }
    } else {
        vout << "Output MTD history file: " << OutputFileName << endl;
        vout << ">>> INFO: No data in MTD history." << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CMTDServer::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# mtd-server terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDServer::CtrlCSignalHandler(int signal)
{
    MTDServer.vout << endl;
    MTDServer.vout << "Ctrl+C/TERM signal recieved. ";
    MTDServer.vout << "Initiating server shutdown ... " << endl;
    MTDServer.vout << "Waiting for server finalization ... " << endl;
    MTDServer.TerminateServer();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
