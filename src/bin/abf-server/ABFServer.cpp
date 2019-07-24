// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <stdio.h>
#include <signal.h>
#include <ErrorSystem.hpp>
#include <SmallTime.hpp>
#include <CmdProcessorList.hpp>
#include <ExtraOperation.hpp>
#include <PMFOperation.hpp>
#include "ABFServer.hpp"
#include "ABFProcessor.hpp"
#include "ABFFactory.hpp"
#include <iomanip>
#include <sstream>
#include <PrmUtils.hpp>
#include <FileSystem.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

CABFServer ABFServer;
MAIN_ENTRY_OBJECT(ABFServer)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFServer::CABFServer(void)
{
    SetProtocolName("abf");
    OutputFileName = "_abfserver.rst";
    AutoRestart =  true;
    AutoSaveInterval = 500;
    SaveCounter = 0;
    TrajFile = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFServer::Init(int argc,char* argv[])
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
    vout << "# abf-server (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
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

    // read [files] section if any
    ProcessFilesControl();

    // read [sync] section if any
    ProcessSyncControl();

    if( Controls.CountULines() > 0 ){
        ES_TRACE_ERROR("unprocessed items found in control file");
        vout << low;
        vout << endl;
        vout << "<b><red> >>> ERROR: Unprocessed items found in the control file!</red></b>" << endl;
        vout << "<b><red>            Their list is provided below.</red></b>" << endl;
        vout << endl;
        Controls.Dump(vout,true);
        return( SO_USER_ERROR );
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

void CABFServer::ProcessFilesControl(void)
{
    vout << endl;
    vout << "=== [files] ====================================================================" << endl;
    if(Controls.OpenSection("files") == false) {
        vout << "Input ABF accumulator (input)             = none                       (default)" << endl;
        vout << "Output ABF accumulator (output)           = " << OutputFileName << "             (default)" << endl;
        vout << "Output ABF trajectory (trajectory)        = none                       (default)" << endl;
        vout << "Restart from the last run (autorestart)   = " << PrmFileOnOff(AutoRestart) << "                         (default)" << endl;
        vout << format("Autosave interval (saveinterval)          = %5d                      (default)\n") % AutoSaveInterval;
        return;
    }

    if(Controls.GetStringByKey("input",InputFileName) == true) {
        vout << "Input ABF accumulator (input)     = " << InputFileName << endl;
    } else {
        vout << "Input ABF accumulator (input)             = none                       (default)" << endl;
    }

    if(Controls.GetStringByKey("output",OutputFileName) == true) {
        vout << "Output ABF accumulator (output)           = " << OutputFileName << endl;
    } else {
        vout << "Output ABF accumulator (output)           = " << OutputFileName << "             (default)" << endl;
    }

    if(Controls.GetStringByKey("trajectory",TrajFileName) == true) {
        vout << "Output ABF trajectory (trajectory)        = " << TrajFileName << endl;
    } else {
        vout << "Output ABF trajectory (trajectory)        = none                       (default)" << endl;
    }

    if(Controls.GetLogicalByKey("autorestart",AutoRestart) == true) {
        vout << "Restart from the last run (autorestart)   = " << PrmFileOnOff(AutoRestart) << endl;
    } else {
        vout << "Restart from the last run (autorestart)   = " << PrmFileOnOff(AutoRestart) << "                         (default)" << endl;
    }

    if(Controls.GetIntegerByKey("saveinterval",AutoSaveInterval) == true) {
        vout << format("Autosave interval (saveinterval)          = %5d\n") % AutoSaveInterval;
    } else {
        vout << format("Autosave interval (saveinterval)          = %5d                      (default)\n") % AutoSaveInterval;
    }
}

//------------------------------------------------------------------------------

void CABFServer::ProcessSyncControl(void)
{
    vout << endl;
    vout << "=== [sync] =====================================================================" << endl;
    if(Controls.OpenSection("sync") == false) {
        vout << "Synchronous mode                  = off                                (default)" << endl;
        return;
    }

    int num_of_clients = 0;
    if(Controls.GetIntegerByKey("clients",num_of_clients) == true) {
        vout << "Number of clients (clients)       = " << num_of_clients << endl;
    } else {
        vout << "Number of clients (clients)       = " << setw(3) << num_of_clients << "                               (default)" << endl;
    }

    ABFAccumulator.SetupSynchronousMode(num_of_clients);

    if( ABFAccumulator.GetNumberOfClients() > 0 ) {
        vout << "Synchronous mode                  = on" << endl;
    } else{
        vout << "Synchronous mode                  = off" << endl;
    }

    // FIXME - Letif
    // do we need any further setup or parameters to mix data in CABFServerAccu::ExchangeDataSynchronousMode?
    // if so then we can read them here or in other more relavant section than in [sync]
}

//------------------------------------------------------------------------------

void CABFServer::AutoSaveData(void)
{
    // lock access
    AutoSaveMutex.Lock();

    SaveCounter++;
    if( AutoSaveInterval <= 0 ){
        AutoSaveMutex.Unlock();
        return; // disabled
    }
    if( SaveCounter % AutoSaveInterval != 0 ){
        AutoSaveMutex.Unlock();
        return; // too soon
    }

    vout << format(" --> %9d autosave to ABF accumulator: %s\n") % SaveCounter % OutputFileName;
    try {
        ABFAccumulator.Save(OutputFileName);
        if( TrajFile ){
            ABFAccumulator.Save(TrajFile);
        }
    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to save ABF output accumulator",e);
        AutoSaveMutex.Unlock();
        return; // skip error
    }

    // unlock access
    AutoSaveMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFServer::Run(void)
{
    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Input data :::::::::::::::::::::::::::::::::" << endl;

    if(InputFileName != NULL) {
        vout << "Input ABF accumulator: " << InputFileName << endl;
        try {
            ABFAccumulator.Load(InputFileName);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to load input ABF accumulator",e);
            return(false);
        }
    } else {
            vout << "Input ABF accumulator: none" << endl;
        if( AutoRestart ){
            vout << "Autorestart is on, searching for the last ABF accumulator to restart from ..." << endl;
            if( CFileSystem::IsFile(OutputFileName) ){
                vout << " == Last ABF accumulator found as:         " << OutputFileName << endl;
                try {
                    ABFAccumulator.Load(OutputFileName);
                } catch(std::exception& e) {
                    ES_ERROR_FROM_EXCEPTION("unable to load restart ABF accumulator",e);
                    return(false);
                }
                CSmallTimeAndDate dt;
                dt.GetActualTimeAndDate();
                std::stringstream backup_name;
                backup_name << format("%s.backuped_%04d-%02d-%02d_%02d:%02d:%02d")
                               % OutputFileName % dt.GetYear() % dt.GetMonth() % dt.GetDay()
                               % dt.GetHour() % dt.GetMinute() % dt.GetSecond();
                vout << " == Backuping the last ABF accumulator as: " << backup_name.str() << endl;
                try {
                    ABFAccumulator.Save(backup_name.str());
                } catch(std::exception& e) {
                    ES_ERROR_FROM_EXCEPTION("unable to save backup of the ABF output accumulator",e);
                    return(false);
                }
            } else {
                vout << " == No last ABF accumulator found, fresh start ..." << endl;
            }
        }
    }

    if(TrajFileName != NULL) {
        vout << "Output ABF trajectory: " << TrajFileName << endl;
        TrajFile = fopen(TrajFileName,"w");
        if( TrajFile == NULL ){
            ES_ERROR("unable to open ABF trajectory file for writing");
            return(false);
        }
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: ABF Server :::::::::::::::::::::::::::::::::" << endl;

// number of worker threads for synchronous mode must be equal to number of clients+1
    if( ABFAccumulator.GetNumberOfClients() > 0 ){
        if( ABFAccumulator.GetNumberOfClients() >= GetMaxNumOfWorkers() ){
            vout << endl;
            vout << ">>> WARNING: In the synchronous mode, number of worker threads [server]/maxthreads" << endl;
            vout << "             must be equal to or greater than the number of clients!" << endl;
            vout << "             [server]/maxthreads is set to " << ABFAccumulator.GetNumberOfClients()+1 << endl;
            SetMaxNumOfWorkers(ABFAccumulator.GetNumberOfClients()+1);
        }
    }

// register operations
    CmdProcessorList.RegisterProcessor(Operation_GetServerInfo,&ABFFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_FlushServerData,&ABFFactory);
    CmdProcessorList.RegisterProcessor(Operation_RegisterClient,&ABFFactory);
    CmdProcessorList.RegisterProcessor(Operation_UnregisterClient,&ABFFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_GetInitialData,&ABFFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_ExchangeData,&ABFFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_GetData,&ABFFactory);

// start server listening
    if( StartServer(vout) == false ){
        return(false);
    }

    // closing trajectory if opened
    if( TrajFile != NULL ){
        fclose(TrajFile);
    }

// print server statistics
    vout << endl;
    vout << "=== Short server summary =======================================================" << endl;
    vout << endl;

    vout << "Number of processed transactions: " << GetNumberOfTransactions() << endl;
    vout << "Number of illegal transactions  : " << GetNumberOfIllegalTransactions() << endl;

    RegClients.PrintInfo(vout);

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;

// and now save all data
    if(ABFAccumulator.GetNumberOfCoords() > 0) {
        vout << "Output ABF accumulator: " << OutputFileName <<  endl;
        try {
            ABFAccumulator.Save(OutputFileName);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to save ABF output accumulator",e);
            return(false);
        }
    } else {
        vout << "Output ABF accumulator: " << OutputFileName <<  endl;
        vout << ">>> INFO: No data in ABF accumulator." << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFServer::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-server terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFServer::CtrlCSignalHandler(int signal)
{
    ABFServer.vout << endl;
    ABFServer.vout << "Ctrl+C/TERM signal recieved. ";
    ABFServer.vout << "Initiating server shutdown ... " << endl;
    ABFServer.vout << "Waiting for server finalization ... " << endl;
    ABFServer.ABFAccumulator.SetServerTerminated();
    ABFServer.TerminateServer();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
