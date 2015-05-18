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
#include "ABFProcessor.hpp"
#include "ABFServer.hpp"
#include <XMLElement.hpp>
#include <ServerCommand.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFProcessor::RegisterClient(void)
{
// set or check coordinates
    if(ABFServer.ABFAccumulator.RegisterOrCheckCoords(Command) == false) {
        LOGIC_ERROR("unable to create or check coordinates");
    }

// in synchronous mode only the same number of clients can register
    if( ABFServer.ABFAccumulator.GetNumberOfClients() > 0 ){
        // how many clients we have now?
        int nregcl = ABFServer.RegClients.GetNumberOfRegClients();
        if( nregcl >= ABFServer.ABFAccumulator.GetNumberOfClients() ){
            CSmallString error;
            error << "in synchronous mode, only '" << ABFServer.ABFAccumulator.GetNumberOfClients() << "' clients are allowed (as requested in [sync] section)";
            LOGIC_ERROR(error);
        }
    }

// register client
    CSmallString job_id;
    CommandElement->GetAttribute("job_id",job_id);
    CRegClient* p_rc = new CRegClient(Command->GetClientName(),job_id);
    ABFServer.RegClients.RegisterClient(p_rc);

// write response
    ResultElement->SetAttribute("client_id",p_rc->GetClientID());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
