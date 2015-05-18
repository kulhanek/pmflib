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

#include <stdio.h>
#include <signal.h>
#include <ErrorSystem.hpp>
#include <SmallTime.hpp>
#include <FileSystem.hpp>
#include <CmdProcessorList.hpp>
#include <PMFOperation.hpp>
#include <ExtraOperation.hpp>
#include "StringServer.hpp"
#include "StringProcessor.hpp"
#include "StringFactory.hpp"
#include <iomanip>
#include <PrmUtils.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

CStringServer StringServer;
MAIN_ENTRY_OBJECT(StringServer)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CStringServer::CStringServer(void)
{
    SetProtocolName("stm");
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CStringServer::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CStringIntOpts
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
    Beads.AttachVerboseStream(vout,Options.GetOptVerbose());

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# stm-server (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgControlFile() != "-") {
        vout << "# Control file name : " << Options.GetArgControlFile() << endl;
    } else {
        vout << "# Control file name : - (standard input)" << endl;
    }

// set SIGINT hadler to cleanly shutdown server ----------
    signal(SIGINT,CtrlCSignalHandler);

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

    if( ProcessServerControl(Controls) == false ){
       return( SO_USER_ERROR );
       }

    // additional setup
    Beads.ProcessSTMControl(Controls);
    Beads.ProcessFilesControl(Controls);
    Beads.ProcessIntervalsControl(Controls);

    // load luncher setup
    if( Launcher.ReadControl(Controls,vout) == false ){
        ES_TRACE_ERROR("unable to load launcher setup");
        return( SO_USER_ERROR );
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::::: {PATHS} ::::::::::::::::::::::::::::::::::" << endl;

    Beads.LoadPath(Controls);

    if( Controls.CountULines() > 0 ){
        vout << endl;
        vout << ":::::::::::::::::::::::::::::::: Unprocessed items :::::::::::::::::::::::::::::" << endl;
        ES_ERROR("unprocessed items found in control file");
        Controls.Dump(stderr,true);
        return( SO_USER_ERROR );
    }

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CStringServer::Run(void)
{
    vout << endl;
    Beads.PrintPathSummary(vout);

    // open trajectory
    if( Beads.OpenTrajectory() == false ){
        return(false);
    }

    vout << endl;
    vout << ":::::::::::::::::::::::::::::::::: String Server :::::::::::::::::::::::::::::::" << endl;

// number of worker threads for synchronous mode must be equal to number of beads
    if( Launcher.IsEnabled() == true ){
        if( Beads.IsAsynchronous() == false ){
            vout << endl;
            vout << ">>> WARNING: The job launcher is activated but the asynchronous mode is not selected!" << endl;
            vout << "             [stm]/async is set to true" << endl;
            Beads.SetAsynchronousMode(true);
        }
    }

    if( Beads.IsAsynchronous() == false ){
        if( Beads.GetNumOfBeads() > GetMaxNumOfWorkers() ){
            vout << endl;
            vout << ">>> WARNING: In the synchronous mode, number of worker threads [server]/maxthreads" << endl;
            vout << "             must be equal to or greater than the number of beads!" << endl;
            vout << "             [server]/maxthreads is set to " << Beads.GetNumOfBeads()+1 << endl;
            SetMaxNumOfWorkers(Beads.GetNumOfBeads()+1);
        }
    }
    if( Beads.IsAsynchronous() == true ){
        if( GetMaxNumOfWorkers() < 2 ){
            vout << endl;
            vout << ">>> WARNING: In the asynchronous mode, number of worker threads [server]/maxthreads" << endl;
            vout << "             must be equal or larger than two!" << endl;
            vout << "             [server]/maxthreads is set to two" << endl;
            SetMaxNumOfWorkers(2);
        }
    }
    if( Beads.IsAsynchronous() == true ){
        if( DoNotShutdown == false ){
            vout << endl;
            vout << ">>> WARNING: In the asynchronous mode, option [server]/donotshutdown must be set to true!" << endl;
            vout << "             [server]/donotshutdown is set to true" << endl;
            DoNotShutdown = true;
        }
    }

// register operations
    CmdProcessorList.RegisterProcessor(Operation_GetServerInfo,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_FlushServerData,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_RegisterClient,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_UnregisterClient,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_GetInitialData,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_ExchangeData,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_GetData,&StringFactory);
    CmdProcessorList.RegisterProcessor(Operation_ShutdownServer,&StringFactory);

// start job launcher
    if( Launcher.IsEnabled() ){
        if( Launcher.StartLauncher(vout) == false ) return(false);
    }

// start server listening
    if( StartServer(vout) == false ){
        return(false);
    }

// terminate job launcher
    if( Launcher.IsEnabled() ){
        Launcher.TerminateThread();
        Launcher.WaitForThread();
        vout << "Launcher server was terminated." << endl;
    }

    vout << endl;
    vout << "Number of processed transactions: " << GetNumberOfTransactions() << endl;
    vout << "Number of illegal transactions  : " << GetNumberOfIllegalTransactions() << endl;

    RegClients.PrintInfo(vout);

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;
    bool result = true;
    result &= Beads.SavePath();
    result &= Beads.SavePathSummary();
    vout << endl;
    Beads.PrintPathSummary(vout);

    return(result);
}

//------------------------------------------------------------------------------

void CStringServer::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# stm-server terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//------------------------------------------------------------------------------

void CStringServer::TerminateServer(void)
{
    // wake up all waiting replicas
    Beads.SetServerTerminated();
    // terminate server
    CServer::TerminateServer();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CStringServer::CtrlCSignalHandler(int signal)
{
    StringServer.vout << endl;
    StringServer.vout << "Ctrl+C signal recieved. ";
    StringServer.vout << "Initiating server shutdown ... " << endl;

    // stop launcher
    if(  StringServer.Launcher.IsEnabled() ) {
        StringServer.vout << "Waiting for launcher server finalization ... " << endl;
        StringServer.Launcher.TerminateThread();
        StringServer.Launcher.WaitForThread();
    }

    StringServer.vout << "Waiting for STM server finalization ... " << endl;

    // terminate server
    StringServer.TerminateServer();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
