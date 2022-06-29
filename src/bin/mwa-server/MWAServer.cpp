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

#include <stdio.h>
#include <signal.h>
#include <ErrorSystem.hpp>
#include <SmallTime.hpp>
#include <CmdProcessorList.hpp>
#include <ExtraOperation.hpp>
#include <PMFOperation.hpp>
#include "MWAServer.hpp"
#include "MWAProcessor.hpp"
#include "MWAFactory.hpp"
#include <iomanip>
#include <sstream>
#include <PrmUtils.hpp>
#include <FileSystem.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

CMWAServer MWAServer;
MAIN_ENTRY_OBJECT(MWAServer)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMWAServer::CMWAServer(void)
{
    SetProtocolName("mwa");
    OutputFileName = "_mwaserver.rst";
    AutoRestart =  true;
    AutoSaveInterval = 50;
    SaveCounter = 0;
    TrajFile = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMWAServer::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CMWAIntOpts
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
    vout << "# mwa-server (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
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

void CMWAServer::ProcessFilesControl(void)
{
    vout << endl;
    vout << "=== [files] ====================================================================" << endl;
    if(Controls.OpenSection("files") == false) {
        vout << "Input MWA accumulator (input)             = none                       (default)" << endl;
        vout << "Output MWA accumulator (output)           = " << OutputFileName << "             (default)" << endl;
        vout << "Output MWA trajectory (trajectory)        = none                       (default)" << endl;
        vout << "Restart from the last run (autorestart)   = " << PrmFileOnOff(AutoRestart) << "                         (default)" << endl;
        vout << format("Autosave interval (saveinterval)          = %5d                      (default)\n") % AutoSaveInterval;
        return;
    }

    if(Controls.GetStringByKey("input",InputFileName) == true) {
        vout << "Input MWA accumulator (input)     = " << InputFileName << endl;
    } else {
        vout << "Input MWA accumulator (input)             = none                       (default)" << endl;
    }

    if(Controls.GetStringByKey("output",OutputFileName) == true) {
        vout << "Output MWA accumulator (output)           = " << OutputFileName << endl;
    } else {
        vout << "Output MWA accumulator (output)           = " << OutputFileName << "             (default)" << endl;
    }

    if(Controls.GetStringByKey("trajectory",TrajFileName) == true) {
        vout << "Output MWA trajectory (trajectory)        = " << TrajFileName << endl;
    } else {
        vout << "Output MWA trajectory (trajectory)        = none                       (default)" << endl;
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

void CMWAServer::ProcessSyncControl(void)
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

    MWAAccumulator.SetupSynchronousMode(num_of_clients);

    if( MWAAccumulator.GetNumberOfClients() > 0 ) {
        vout << "Synchronous mode                  = on" << endl;
    } else{
        vout << "Synchronous mode                  = off" << endl;
    }

    // FIXME - Letif
    // do we need any further setup or parameters to mix data in CMWAServerAccu::ExchangeDataSynchronousMode?
    // if so then we can read them here or in other more relavant section than in [sync]
}

//------------------------------------------------------------------------------

void CMWAServer::AutoSaveData(void)
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

    vout << format(" --> %9d autosave to MWA accumulator: %s\n") % SaveCounter % OutputFileName;
    try {
        MWAAccumulator.AccuMutex.Lock();
        MWAAccumulator.Save(OutputFileName);
        if( TrajFile ){
            MWAAccumulator.Save(TrajFile);
        }
        MWAAccumulator.AccuMutex.Unlock();
    } catch(std::exception& e) {
        MWAAccumulator.AccuMutex.Unlock();
        AutoSaveMutex.Unlock();
        ES_ERROR_FROM_EXCEPTION("unable to save MWA output accumulator",e);
        return; // skip error
    }

    // unlock access
    AutoSaveMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMWAServer::Run(void)
{
    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Input data :::::::::::::::::::::::::::::::::" << endl;

    if(InputFileName != NULL) {
        vout << "Input MWA accumulator: " << InputFileName << endl;
        try {
            MWAAccumulator.Load(InputFileName);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to load input MWA accumulator",e);
            return(false);
        }
    } else {
            vout << "Input MWA accumulator: none" << endl;
        if( AutoRestart ){
            vout << "Autorestart is on, searching for the last MWA accumulator to restart from ..." << endl;
            if( CFileSystem::IsFile(OutputFileName) ){
                vout << " == Last MWA accumulator found as:         " << OutputFileName << endl;
                try {
                    MWAAccumulator.Load(OutputFileName);
                } catch(std::exception& e) {
                    ES_ERROR_FROM_EXCEPTION("unable to load restart MWA accumulator",e);
                    return(false);
                }
                CSmallTimeAndDate dt;
                dt.GetActualTimeAndDate();
                std::stringstream backup_name;
                backup_name << format("%s.backuped_%04d-%02d-%02d_%02d:%02d:%02d")
                               % OutputFileName % dt.GetYear() % dt.GetMonth() % dt.GetDay()
                               % dt.GetHour() % dt.GetMinute() % dt.GetSecond();
                vout << " == Backuping the last MWA accumulator as: " << backup_name.str() << endl;
                try {
                    MWAAccumulator.Save(backup_name.str());
                } catch(std::exception& e) {
                    ES_ERROR_FROM_EXCEPTION("unable to save backup of the MWA output accumulator",e);
                    return(false);
                }
            } else {
                vout << " == No last MWA accumulator found, fresh start ..." << endl;
            }
        }
    }

    if( MWAAccumulator.GetNumOfCVs() > 0 ){
        MWAAccumulator.PrintAccuInfo(vout);
        MWAAccumulator.PrintCVSInfo(vout);
    }

    if(TrajFileName != NULL) {
        vout << "Output MWA trajectory: " << TrajFileName << endl;
        TrajFile = fopen(TrajFileName,"w");
        if( TrajFile == NULL ){
            ES_ERROR("unable to open MWA trajectory file for writing");
            return(false);
        }
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: MWA Server :::::::::::::::::::::::::::::::::" << endl;

// number of worker threads for synchronous mode must be equal to number of clients+1
    if( MWAAccumulator.GetNumberOfClients() > 0 ){
        if( MWAAccumulator.GetNumberOfClients() >= GetMaxNumOfWorkers() ){
            vout << endl;
            vout << ">>> WARNING: In the synchronous mode, number of worker threads [server]/maxthreads" << endl;
            vout << "             must be equal to or greater than the number of clients!" << endl;
            vout << "             [server]/maxthreads is set to " << MWAAccumulator.GetNumberOfClients()+1 << endl;
            SetMaxNumOfWorkers(MWAAccumulator.GetNumberOfClients()+1);
        }
    }

// register operations
    CmdProcessorList.RegisterProcessor(Operation_GetServerInfo,&MWAFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_FlushServerData,&MWAFactory);
    CmdProcessorList.RegisterProcessor(Operation_RegisterClient,&MWAFactory);
    CmdProcessorList.RegisterProcessor(Operation_UnregisterClient,&MWAFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_GetInitialData,&MWAFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_ExchangeData,&MWAFactory);
    CmdProcessorList.RegisterProcessor(OperationPMF_GetData,&MWAFactory);

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
    if(MWAAccumulator.GetNumOfCVs() > 0) {
        vout << "Output MWA accumulator: " << OutputFileName <<  endl;
        try {
            MWAAccumulator.Save(OutputFileName);
        } catch(std::exception& e) {
            ES_ERROR_FROM_EXCEPTION("unable to save MWA output accumulator",e);
            return(false);
        }
    } else {
        vout << "Output MWA accumulator: " << OutputFileName <<  endl;
        vout << ">>> INFO: No data in MWA accumulator." << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CMWAServer::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# mwa-server terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMWAServer::CtrlCSignalHandler(int signal)
{
    MWAServer.vout << endl;
    MWAServer.vout << "Ctrl+C/TERM signal recieved. ";
    MWAServer.vout << "Initiating server shutdown ... " << endl;
    MWAServer.vout << "Waiting for server finalization ... " << endl;
    MWAServer.MWAAccumulator.SetServerTerminated();
    MWAServer.TerminateServer();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
