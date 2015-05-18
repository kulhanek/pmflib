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

void CStringProcessor::UnregisterClient(void)
{
    int client_id = -1;
    int bead_id = -1;

    CBead* p_bead = NULL;
    if(CommandElement->GetAttribute("bead_id",bead_id) == true) {
        // request to release the bead
        p_bead = StringServer.Beads.GetBead(bead_id);
        client_id = p_bead->GetClientID();
    } else {
        // request to unregister the client
        // get client ID
        if(CommandElement->GetAttribute("client_id",client_id) == false) {
            LOGIC_ERROR("unable to get client_id");
        }
        // remove client registration from the bead
        p_bead = StringServer.Beads.GetBeadByClientID(client_id);
    }

    bool release = false;
    CommandElement->GetAttribute("release",release);

    if( release ){
        p_bead->ReleaseBead();
    }
    p_bead->SetClientID(-1);

//  and unregister it from list
    StringServer.RegClients.UnregisterClient(client_id);

// automatically shutdown server if last client is unregistered
    if(StringServer.RegClients.GetNumberOfActiveRegistration() == 0) {
        if(StringServer.DoNotShutdown == false) StringServer.TerminateServer();
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
