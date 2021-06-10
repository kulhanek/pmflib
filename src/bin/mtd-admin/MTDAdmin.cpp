// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
// ===============================================================================

#include <stdio.h>
#include "MTDAdmin.hpp"
#include <ErrorSystem.hpp>
#include <MTDAccumulator.hpp>
#include <PMFOperation.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CMTDAdmin)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDAdmin::CMTDAdmin(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDAdmin::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to terminal stream and set desired verbosity level
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
    vout << "# mtd-admin (PMFLib utility) started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# Server request (by user)   : " << Options.GetArgCommand() << endl;

// set server description - there should no be any problem
// syntax was already checked
    ActionRequest.SetProtocolName("mtd");
    ActionRequest.SetQualifiedName(Options.GetArgCommand());

    if( Options.IsOptServerKeySet() || ActionRequest.IsServerKey() ) {
        ActionRequest.ReadServerKey(Options.GetOptServerKey());
    }

// print header --------------------------------------------------------------
    vout << "# Server request (processed) : " << ActionRequest.GetQualifiedName() << endl;

    vout << debug;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << "# Server name                : " << ActionRequest.GetName() << endl;
    vout << "# Server IP                  : " << ActionRequest.GetIP() << endl;
    vout << "# Server port                : " << ActionRequest.GetPort() << endl;
    vout << "# Command                    : " << ActionRequest.GetAction() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    CSmallString passfile;
    bool         passfile_set;
    passfile = Options.GetOptPassword(passfile_set);
    if(passfile_set == true) {
        vout << "# Password                   : from file '" << passfile << "'" << endl;
    } else {
        vout << "# Password                   : from command line or server key file" << endl;
    }
    vout << high;

    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDAdmin::Run(void)
{
    if( ! ( Options.IsOptServerKeySet() || ActionRequest.IsServerKey() )  ){
       if( ReadPassword(Options.IsOptPasswordSet(),Options.GetOptPassword(),
                        Options.GetOptVerbose()) == false ) return(false);
       }

    vout << debug;
    vout << endl;
    vout << "Sending request ... " << ActionRequest.GetAction() << endl << endl;
    vout << high;

    bool result;

    if(ActionRequest.GetAction() == "info") {
        result = GetServerInfo(vout);
    } else if(ActionRequest.GetAction() == "shutdown") {
        result = ShutdownServer(vout);
    } else if(ActionRequest.GetAction() == "flush") {
        result = FlushServerData();
    } else if(ActionRequest.GetAction() == "errors") {
        result = GetServerErrors(vout);
    } else if(ActionRequest.GetAction() == "get") {
        result = GetMTDHistory();
    } else {
        CSmallString error;
        error << "action " << ActionRequest.GetAction() << " is not implemented";
        ES_ERROR(error);
        return(false);
    }

    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDAdmin::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# mtd-admin terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDAdmin::GetMTDHistory(void)
{
    CClientCommand      cmd;
    CMTDAccumulator     accu;

    try{

        // init command
        InitCommand(&cmd,OperationPMF_GetData);

        // execute command
        ExecuteCommand(&cmd);

        // process response
        CXMLElement* p_rele = cmd.GetRootResultElement();
        if(p_rele == NULL) {
            LOGIC_ERROR("unable to open RESULTS element");
        }

        // read cvs info
        accu.LoadCVSInfo(p_rele);
        accu.ReadMTDData(p_rele);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;

    CSmallString file_output;

    if(ActionRequest.GetParameterKeyValue("file",file_output) == false) {
        file_output = "_mtdserver.rst";
    }

// and now save all data
    if(accu.GetNumOfCVs() > 0) {
        vout << "Output MTD accumulator: " << file_output << endl;
        try {
            accu.Save(file_output);
        } catch(...) {
            ES_ERROR("unable to save MTD accumulator file");
        }
    } else {
        vout << "Output MTD accumulator: " << file_output << endl;
        vout << ">>> INFO: No data in MTD accumulator." << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

