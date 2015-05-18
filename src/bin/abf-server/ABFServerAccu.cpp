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

#include <string.h>
#include "ABFServerAccu.hpp"
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <ServerCommand.hpp>
#include <CmdProcessor.hpp>
#include "ABFServer.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFServerAccu::CABFServerAccu(void)
{
    SynchronousMode = false;
    NumOfClients = 0;
    ServerTerminated = false;
}

//------------------------------------------------------------------------------

void CABFServerAccu::SetupSynchronousMode(int numofclients)
{
    if( numofclients > 0 ){
        SynchronousMode = true;
    } else {
        SynchronousMode = false;
    }
    NumOfClients = numofclients;
}

//------------------------------------------------------------------------------

int CABFServerAccu::GetNumberOfClients(void)
{
    return(NumOfClients);
}

//------------------------------------------------------------------------------

void CABFServerAccu::SetServerTerminated(void)
{
    ServerTerminated = true;
    RendezvousCond.BroadcastSignal();
}

//------------------------------------------------------------------------------

bool CABFServerAccu::RegisterOrCheckCoords(CServerCommand* p_cmd)
{
    CXMLElement* p_cele = p_cmd->GetRootCommandElement();
    CXMLElement* p_cvs = p_cele->GetFirstChildElement("CVS");
    if(p_cvs == NULL) {
        LOGIC_ERROR("unable to open CVS element");
    }

    bool result = true;

    AccuMutex.Lock();

    if(GetNumberOfCoords() == 0) {
        try {
            // load coordinates from first client
            LoadCVSInfo(p_cvs);

            // print info about CVs
            PrintCVSInfo(ABFServer.vout);

        } catch(std::exception& e) {
            CMD_ERROR_FROM_EXCEPTION(p_cmd,"unable to load coordinates",e);
            result = false;
        }
    } else {
        // coodinates have been loaded - so check if other clients are the same
        result = CheckCVSInfo(p_cvs);
        if( result == false ){
            // copy last reason to command error stack
            CMD_ERROR(p_cmd,ErrorSystem.GetLastError());
        }
    }

    AccuMutex.Unlock();

    return(result);
}

//------------------------------------------------------------------------------

void CABFServerAccu::GetInitialData(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_rele = p_command->GetRootResultElement();

    AccuMutex.Lock();
        try {
            WriteABFData(p_rele);
        } catch(...) {
            AccuMutex.Unlock();
            throw;
        }
    AccuMutex.Unlock();
}

//------------------------------------------------------------------------------

void CABFServerAccu::ExchangeDataWithClient(CRegClient* p_client)
{
    if(p_client == NULL) {
        INVALID_ARGUMENT("p_client is NULL");
    }

    if( ServerTerminated ){
        // server was terminated prematurely
        return;
    }

    // check if the data storage are ready
    // FIXME - for Letif add weights array check
    if( (GetNSamplesArray() == NULL)
        || (GetABFForceSumArray() == NULL)
        || (GetABFForceSquareSumArray() == NULL)) {
        LOGIC_ERROR("data array(s) is(are) NULL");
    }

    if( SynchronousMode ){
        ExchangeDataSynchronousMode(p_client);
    } else {
        ExchangeDataAsynchronousMode(p_client);
    }
}

//------------------------------------------------------------------------------

void CABFServerAccu::ExchangeDataAsynchronousMode(CRegClient* p_client)
{
    if(p_client == NULL) {
        INVALID_ARGUMENT("p_client is NULL");
    }
    CServerCommand* p_command = p_client->GetCommand();
    if(p_command == NULL) {
        RUNTIME_ERROR("p_commad is NULL");
    }

    CXMLElement* p_cele = p_command->GetRootCommandElement();
    CXMLElement* p_rele = p_command->GetRootResultElement();

    AccuMutex.Lock();
    try {
        AddABFData(p_cele);
        WriteABFData(p_rele);
    } catch(...) {
        AccuMutex.Unlock();
        throw;
    }
    AccuMutex.Unlock();

}

//------------------------------------------------------------------------------

void CABFServerAccu::ExchangeDataSynchronousMode(CRegClient* p_client)
{
    if(p_client == NULL) {
        INVALID_ARGUMENT("p_client is NULL");
    }
    CServerCommand* p_command = p_client->GetCommand();
    if(p_command == NULL) {
        RUNTIME_ERROR("p_commad is NULL");
    }

    RendezvousMutex.Lock();
        CServerCommand* p_cmd = p_client->GetCommand();

        // add command into command queue
        Commands.InsertToEnd(p_cmd);

        // do we need to wait for other clients?
        if( Commands.NumOfMembers() != NumOfClients ){
            p_client->SetExtraStatus("W");
            // we need to wait
            RendezvousCond.WaitForSignal(RendezvousMutex);
            // return status back
            p_client->SetExtraStatus("");
        } else {
            // this is the last client that will do all work for all clients
            try {
                AccuMutex.Lock();
                // FIXME - Letif
                // do all necessary operations on all clients
                // client command list is accessible in Commands list

                // this code performs convential ABF
                CSimpleIterator<CServerCommand> I(Commands);

                // for each client add data to central accumulator
                while( (p_cmd = I.Current()) != NULL ){
                    CXMLElement* p_cele = p_cmd->GetRootCommandElement();
                    AddABFData(p_cele);
                    I++;
                }

                // then for each client return back the whole ABF accumulator
                I.SetToBegin();
                while( (p_cmd = I.Current()) != NULL ){
                    CXMLElement* p_rele = p_cmd->GetRootResultElement();
                    WriteABFData(p_rele);
                    I++;
                }

                AccuMutex.Unlock();
            } catch(...) {
                Commands.RemoveAll();
                RendezvousCond.BroadcastSignal();
                AccuMutex.Unlock();
                RendezvousMutex.Unlock();
                throw;
            }
            // unblock all other clients
            RendezvousCond.BroadcastSignal();
            Commands.RemoveAll();
        }
    RendezvousMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFServerAccu::GetData(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_cvs = p_command->GetResultElementByPath("CVS",true);
    CXMLElement* p_data = p_command->GetResultElementByPath("DATA",true);

    AccuMutex.Lock();
    try {
        SaveCVSInfo(p_cvs);
        if( GetNumberOfCoords() > 0 ) {
            WriteABFData(p_data);
        }
    } catch(...) {
        AccuMutex.Unlock();
        throw;
    }
    AccuMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFServerAccu::FlushData(const CSmallString& name)
{
    AccuMutex.Lock();
        try {
            Save(name);
        } catch(...){
            AccuMutex.Unlock();
            throw;
        }
    AccuMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

