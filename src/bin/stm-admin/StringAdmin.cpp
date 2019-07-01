// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "StringAdmin.hpp"
#include <ExtraOperation.hpp>
#include <BeadList.hpp>
#include <PMFOperation.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CStringAdmin)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CStringAdmin::CStringAdmin(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CStringAdmin::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CStringIntOpts
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
    vout << "# stm-admin (PMFLib utility) started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# Server request (by user)   : " << Options.GetArgCommand() << endl;

// set server description - there should no be any problem
// syntax was already checked
    ActionRequest.SetProtocolName("stm");
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

bool CStringAdmin::Run(void)
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
        result = GetStringPath();
    } else if(ActionRequest.GetAction() == "release") {
        result = ReleaseBead();
    } else if(ActionRequest.GetAction() == "terminate") {
        result = TerminateServer();
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

void CStringAdmin::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# stm-admin terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//------------------------------------------------------------------------------

bool CStringAdmin::GetServerInfo(std::ostream& vout)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_GetServerInfo);

        // execute command
        ExecuteCommand(&cmd);

        // print response
        CBeadList beads;
        beads.LoadInfo(cmd.GetRootResultElement());

        vout << endl;
        beads.PrintPathSummary(vout);

        vout << endl;
        beads.PrintPathUpdate(vout);

        GetShortServerInfo(&cmd,vout);
        GetLongServerInfo(&cmd,vout);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CStringAdmin::GetStringPath(void)
{
    CSmallString file_output;

    if(ActionRequest.GetParameterKeyValue("file",file_output) == false) {
        file_output = "_stm.results";
    }

    bool result = true;

    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_GetData);

        // execute command
        ExecuteCommand(&cmd);

        // print response
        CBeadList beads;
        beads.LoadInfo(cmd.GetRootResultElement());

        vout << endl;
        vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;
        vout << "Output STM path results : " << file_output <<  endl;
        try {
            beads.SavePathSummary(file_output);
        } catch(std::exception& e) {
            result = false;
            ES_ERROR_FROM_EXCEPTION("unable to save STM path",e);
        }

        vout << endl;
        beads.PrintPathSummary(vout);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(result);
}

//------------------------------------------------------------------------------

bool CStringAdmin::ReleaseBead(void)
{
    int bead_id;

    if(ActionRequest.GetParameterKeyValue("id",bead_id) == false) {
        ES_ERROR("bead id must be specified");
        return(false);
    }

    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_UnregisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("bead_id",bead_id);
        p_ele->SetAttribute("release",true);

        // execute command
        ExecuteCommand(&cmd);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

#define PASS_LENGTH 5

bool CStringAdmin::TerminateServer(void)
{
    CSmallString passphrase;

    // set random generator
    time_t now = time(NULL);
    srand(now);
    passphrase.SetLength(PASS_LENGTH);
    int pos = 0;
    int attempts = 0;
    do {
        attempts++;
        if( attempts % 1000 == 0 ){
            time_t now = time(NULL);
            srand(now);
        }
        char c = rand()%(126-33) + 33;
        if( (c == '*') || (c == '#') || (c == '!') || (c == '"') || (c == '\'') ) continue;
        passphrase[pos++] = c;
    } while(pos < PASS_LENGTH);

    vout << endl;
    vout << "Generated security passphrase : " << passphrase << endl;
    vout << "Rewrite passphrase to softly terminate the server: ";
    string user;
    cin >> user;
    if( user != string(passphrase) ){
        vout << "Passphrases do not match - exiting ..." << endl;
        return(false);
    } else {
        vout << "Passphrases match - sending soft termination signal to the server ..." << endl;
    }

    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_ShutdownServer);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("soft",true);

        // execute command
        ExecuteCommand(&cmd);

        // print response
        GetShortServerInfo(&cmd,vout);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

