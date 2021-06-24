// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301  USA
// =============================================================================

#include <MTDClient.hpp>
#include <PMFOperation.hpp>
#include <ClientCommand.hpp>
#include <ErrorSystem.hpp>
#include <ColVariable.hpp>
#include <string.h>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <XMLPrinter.hpp>
#include <ExtraOperation.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDClient      MTDClient;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDClient::CMTDClient(void)
{
    ClientID = 0;
}

//------------------------------------------------------------------------------

CMTDClient::~CMTDClient(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDClient::GetClientID(void)
{
    return(ClientID);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDClient::RegisterClient(void)
{
    ClientID = -1;

    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_RegisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        CSmallString job_id = getenv("INF_JOB_ID");
        if( job_id != NULL ) {
            p_ele->SetAttribute("job_id",job_id);
        }
        p_ele = cmd.GetCommandElementByPath("CVS",true);
        SaveCVSInfo(p_ele);

        // execute command
        ExecuteCommand(&cmd);

        // process response
        p_ele = cmd.GetRootResultElement();
        if( p_ele->GetAttribute("client_id",ClientID) == false ) {
            LOGIC_ERROR("unable to get client_id");
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(-1);
    }

    return(ClientID);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDClient::UnregisterClient(void)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_UnregisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);

        // execute command
        ExecuteCommand(&cmd);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDClient::GetInitialData(void)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_GetInitialData);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);

        // execute command
        ExecuteCommand(&cmd);

        // process results
        CXMLElement* p_aele = cmd.GetResultElementByPath("ACCU",false);
        if( p_aele == NULL ){
            throw("ACCU element not found");
        }
        Load(p_aele);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDClient::ExchangeData(void)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_ExchangeData);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);

        CXMLElement* p_aele = cmd.GetResultElementByPath("ACCU",true);
        Save(p_aele);

        // execute command
        ExecuteCommand(&cmd);

        // process results
        p_aele = cmd.GetResultElementByPath("ACCU",false);
        if( p_ele == NULL ){
            throw("ACCU element not found");
        }
        Load(p_aele);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

