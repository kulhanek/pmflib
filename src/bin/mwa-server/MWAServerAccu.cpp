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

#include <string.h>
#include "MWAServerAccu.hpp"
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <ServerCommand.hpp>
#include <CmdProcessor.hpp>
#include "MWAServer.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMWAServerAccu::CMWAServerAccu(void)
{
    SynchronousMode = false;
    NumOfClients = 0;
    ServerTerminated = false;
}

//------------------------------------------------------------------------------

void CMWAServerAccu::SetupSynchronousMode(int numofclients)
{
    if( numofclients > 0 ){
        SynchronousMode = true;
    } else {
        SynchronousMode = false;
    }
    NumOfClients = numofclients;
}

//------------------------------------------------------------------------------

int CMWAServerAccu::GetNumberOfClients(void)
{
    return(NumOfClients);
}

//------------------------------------------------------------------------------

void CMWAServerAccu::SetServerTerminated(void)
{
    ServerTerminated = true;
    RendezvousCond.BroadcastSignal();
}

//------------------------------------------------------------------------------

bool CMWAServerAccu::RegisterOrCheckCoords(CServerCommand* p_cmd)
{
    CXMLElement* p_accu = p_cmd->GetCommandElementByPath("PMFLIB-V6",false);
    if(p_accu == NULL) {
        LOGIC_ERROR("unable to open PMFLIB-V6 element");
    }
    bool result = true;

    AccuMutex.Lock();

    if(GetNumOfCVs() == 0) {
        try {
            // load coordinates from first client
            Load(p_accu);

            // print info about PMF accumulator and CVs
            PrintAccuInfo(MWAServer.vout);
            PrintCVSInfo(MWAServer.vout);
            MWAServer.vout << std::endl;

        } catch(std::exception& e) {
            CMD_ERROR_FROM_EXCEPTION(p_cmd,"unable to load coordinates",e);
            result = false;
        }
    } else {
        try {
            CPMFAccumulatorPtr accu = CPMFAccumulatorPtr(new CPMFAccumulator);
            accu->Load(p_accu);
            // coordinates have been loaded - so check if other clients are the same
            result = CheckCVSInfo(accu);
        } catch(std::exception& e) {
            CMD_ERROR_FROM_EXCEPTION(p_cmd,"unable to check coordinates",e);
            result = false;
        }
        if( result == false ){
            // copy last reason to command error stack
            CMD_ERROR(p_cmd,ErrorSystem.GetLastError());
        }
    }

    AccuMutex.Unlock();

    return(result);
}

//------------------------------------------------------------------------------

void CMWAServerAccu::GetInitialData(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_root = p_command->GetRootResultElement();

    AccuMutex.Lock();
        try {
            Save(p_root);
        } catch(...) {
            AccuMutex.Unlock();
            throw;
        }
    AccuMutex.Unlock();
}

//------------------------------------------------------------------------------

void CMWAServerAccu::ExchangeDataWithClient(CRegClient* p_client)
{
    if(p_client == NULL) {
        INVALID_ARGUMENT("p_client is NULL");
    }

    if( ServerTerminated ){
        // server was terminated prematurely
        return;
    }

    if( SynchronousMode ){
        ExchangeDataSynchronousMode(p_client);
    } else {
        ExchangeDataAsynchronousMode(p_client);
    }
}

//------------------------------------------------------------------------------

void CMWAServerAccu::ExchangeDataAsynchronousMode(CRegClient* p_client)
{
    if(p_client == NULL) {
        INVALID_ARGUMENT("p_client is NULL");
    }
    CServerCommand* p_command = p_client->GetCommand();
    if(p_command == NULL) {
        RUNTIME_ERROR("p_commad is NULL");
    }

    CXMLElement* p_cele = p_command->GetCommandElementByPath("PMFLIB-V6",false);
    if(p_cele == NULL) {
        LOGIC_ERROR("unable to open PMFLIB-V6 element");
    }
    CXMLElement* p_rele = p_command->GetRootResultElement();

    AccuMutex.Lock();
    try {
        Combine(p_cele);
        Save(p_rele);
    } catch(...) {
        AccuMutex.Unlock();
        throw;
    }
    AccuMutex.Unlock();

}

//------------------------------------------------------------------------------

void CMWAServerAccu::ExchangeDataSynchronousMode(CRegClient* p_client)
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

                // this code performs conventional MWA
                CSimpleIterator<CServerCommand> I(Commands);

                // for each client add data to central accumulator
                while( (p_cmd = I.Current()) != NULL ){
                    CXMLElement* p_cele = p_cmd->GetCommandElementByPath("PMFLIB-V6",false);
                    if( p_cele ){
                        Combine(p_cele);
                    }
                    I++;
                }

                // then for each client return back the whole MWA accumulator
                I.SetToBegin();
                while( (p_cmd = I.Current()) != NULL ){
                    CXMLElement* p_rele = p_cmd->GetResultElementByPath("PMFLIB-V6",true);
                    Save(p_rele);
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

void CMWAServerAccu::GetData(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_accu = p_command->GetRootResultElement();

    AccuMutex.Lock();
    try {
        Save(p_accu);
    } catch(...) {
        AccuMutex.Unlock();
        throw;
    }
    AccuMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMWAServerAccu::FlushData(const CSmallString& name)
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

