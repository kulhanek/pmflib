// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "ABPProcessor.hpp"
#include "ABPServer.hpp"
#include <XMLElement.hpp>
#include <ServerCommand.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPProcessor::RegisterClient(void)
{
// set or check coordinates
    if(ABPServer.ABPAccumulator.RegisterOrCheckCoords(Command) == false) {
        LOGIC_ERROR("unable to create or check coordinates");
    }

// in synchronous mode only the same number of clients can register
    if( ABPServer.ABPAccumulator.GetNumberOfClients() > 0 ){
        // how many clients we have now?
        int nregcl = ABPServer.RegClients.GetNumberOfRegClients();
        if( nregcl >= ABPServer.ABPAccumulator.GetNumberOfClients() ){
            CSmallString error;
            error << "in synchronous mode, only '" << ABPServer.ABPAccumulator.GetNumberOfClients() << "' clients are allowed (as requested in [sync] section)";
            LOGIC_ERROR(error);
        }
    }

// register client
    CSmallString job_id;
    CommandElement->GetAttribute("job_id",job_id);
    CRegClient* p_rc = new CRegClient(Command->GetClientName(),job_id);
    ABPServer.RegClients.RegisterClient(p_rc);

// write response
    ResultElement->SetAttribute("client_id",p_rc->GetClientID());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
