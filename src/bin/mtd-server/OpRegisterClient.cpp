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

#include <stdio.h>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <RegClient.hpp>
#include <ServerCommand.hpp>
#include "MTDProcessor.hpp"
#include "MTDServer.hpp"
#include "MTDRegClient.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDProcessor::RegisterClient(void)
{
// set or check coordinates
    if(MTDServer.MTDHistory.RegisterOrCheckCoords(Command) == false) {
        LOGIC_ERROR("unable to create or check coordinates");
    }

// register client
    CSmallString job_id;
    CommandElement->GetAttribute("job_id",job_id);
    CMTDRegClient* p_rc = new CMTDRegClient(Command->GetClientName(),job_id);
    MTDServer.RegClients.RegisterClient(p_rc);

// now register MTDBuffers that are new for client
    MTDServer.MTDHistory.RegisterAllBuffersAsNew(p_rc);

// write response
    ResultElement->SetAttribute("client_id",p_rc->GetClientID());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
