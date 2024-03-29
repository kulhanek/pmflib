// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <ErrorSystem.hpp>
#include "MWAProcessor.hpp"
#include "MWAServer.hpp"
#include "MWAServerAccu.hpp"
#include <XMLElement.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMWAProcessor::ExchangeData(void)
{
    int client_id = -1;

    if(CommandElement->GetAttribute("client_id",client_id) == false) {
        LOGIC_ERROR("unable to get client_id");
    }

    CRegClient* p_client = MWAServer.RegClients.FindClient(client_id);
    if(p_client == NULL) {
        RUNTIME_ERROR("unable to find client");
    }

    // assign command to client
    if( p_client->GetCommand() != NULL ){
        CSmallString error;
        error << "client '" << client_id << "' has already assigned command";
        RUNTIME_ERROR(error)
    }
    p_client->SetCommand(Command);

    // exchange data
    try{
        MWAServer.MWAAccumulator.ExchangeDataWithClient(p_client);
        p_client->RegisterOperation();
    } catch(...){
        p_client->SetCommand(NULL);
        throw;
    }

    // autosave data
    MWAServer.AutoSaveData();

    // release command
    p_client->SetCommand(NULL);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

