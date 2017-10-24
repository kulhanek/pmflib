// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
#include <ErrorSystem.hpp>
#include <ABFAccumulator.hpp>
#include "ABFAdmin.hpp"
#include <PMFOperation.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFAdmin)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFAdmin::CABFAdmin(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAdmin::Init(int argc,char* argv[])
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

// print header --------------------------------------------------------------
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-admin (PMFLib utility) started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# Server request (by user)   : " << Options.GetArgCommand() << endl;

// set server description - there should no be any problem
// syntax was already checked
    ActionRequest.SetProtocolName("abf");
    ActionRequest.SetQualifiedName(Options.GetArgCommand());

    if(Options.IsOptServerKeySet()) {
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

bool CABFAdmin::Run(void)
{
    if( ! Options.IsOptServerKeySet() ){
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
        result = ShutdownServer(vout,Options.GetOptForce());
    } else if(ActionRequest.GetAction() == "flush") {
        result = FlushServerData();
    } else if(ActionRequest.GetAction() == "errors") {
        result = GetServerErrors(vout);
    } else if(ActionRequest.GetAction() == "get") {
        result = GetABFAccumulator();
    }  else {
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

void CABFAdmin::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-admin terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }
    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFAdmin::GetABFAccumulator(void)
{
    CClientCommand  cmd;
    CABFAccumulator accu;

    try{

        // init command
        InitCommand(&cmd,OperationPMF_GetData);

        // execute command
        ExecuteCommand(&cmd);

        // process response
        CXMLElement* p_cvs = cmd.GetResultElementByPath("CVS",false);
        if(p_cvs == NULL) {
            LOGIC_ERROR("unable to open CVS element");
        }

        CXMLElement* p_data = cmd.GetResultElementByPath("DATA",false);
        if(p_data == NULL) {
            LOGIC_ERROR("unable to open DATA element");
        }

        // read cvs info
        accu.LoadCVSInfo(p_cvs);
        accu.ReadABFData(p_data);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;

    CSmallString file_output;

    if(ActionRequest.GetParameterKeyValue("file",file_output) == false) {
        file_output = "_abfserver.rst";
    }

// and now save all data
    if(accu.GetNumberOfCoords() > 0) {
        vout << "Output ABF accumulator: " << file_output << endl;
        try {
            accu.Save(file_output);
        } catch(...) {
            ES_ERROR("unable to save ABF output accumulator");
        }
    } else {
        vout << "Output ABF accumulator: " << file_output << endl;
        vout << ">>> INFO: No data in ABF accumulator." << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

