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
#include <XMLElement.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMWAProcessor::UnregisterClient(void)
{
    int client_id = -1;

// get client ID
    if(CommandElement->GetAttribute("client_id",client_id) == false) {
        LOGIC_ERROR("unable to get client_id");
    }

//  and unregister it from list
    MWAServer.RegClients.UnregisterClient(client_id);

// automatically shutdown server if last client is unregistered
    if(MWAServer.RegClients.GetNumberOfActiveRegistration() == 0) {
        if(MWAServer.DoNotShutdown == false) MWAServer.TerminateServer();
    }

// did we reach the target number of processed clients?
    if( MWAServer.TargetRegs > 0 ){
        if( (MWAServer.RegClients.GetNumberOfClients() >= MWAServer.TargetRegs) &&
            (MWAServer.RegClients.GetNumberOfActiveRegistration() == 0) ){
            MWAServer.TerminateServer();
        }
    }

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
