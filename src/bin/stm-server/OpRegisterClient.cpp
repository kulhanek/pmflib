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
#include <RegClient.hpp>
#include <ServerCommand.hpp>
#include "StringProcessor.hpp"
#include "StringServer.hpp"
#include <XMLElement.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CStringProcessor::RegisterClient(void)
{
// check if client is eligible to connect
    if( StringServer.Beads.CheckClient(CommandElement) == false ){
        RUNTIME_ERROR("client is not eligible to connect");
    }

// register client
    CSmallString job_id;
    CommandElement->GetAttribute("job_id",job_id);
    CRegClient* p_rc = new CRegClient(Command->GetClientName(),job_id);
    StringServer.RegClients.RegisterClient(p_rc);

// find bead and connect client
    int bead_id = -1;
    CommandElement->GetAttribute("bead_id",bead_id);
    StringServer.Beads.RegisterBead(bead_id,p_rc->GetClientID());

// write response
    ResultElement->SetAttribute("client_id",p_rc->GetClientID());   
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
