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
#include <ErrorSystem.hpp>
#include "StringProcessor.hpp"
#include "StringServer.hpp"
#include <XMLElement.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CStringProcessor::ExchangeData(void)
{
    int client_id = -1;

// get client ID
    if(CommandElement->GetAttribute("client_id",client_id) == false) {
        RUNTIME_ERROR("unable to get client_id");
    }

    CRegClient* p_client = StringServer.RegClients.FindClient(client_id);
    if(p_client == NULL) {
        RUNTIME_ERROR("unable to find client");
    }

// exchange data
    StringServer.Beads.ExchangeData(CommandElement,ResultElement);

    p_client->RegisterOperation();

    if( StringServer.Beads.IsAsynchronous() == false ){
        if( (StringServer.Beads.GetSTMStatus() == ESTMS_MAX_STEPS_REACHED) ||
            (StringServer.Beads.GetSTMStatus() == ESTMS_COMPLETED) ) {
            StringServer.TerminateServer();
        }
    }
    // termination in asynchronous mode is done by Launcher
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

