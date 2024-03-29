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

void CMWAProcessor::GetInitialData(void)
{
    int client_id = -1;

    if(CommandElement->GetAttribute("client_id",client_id) == false) {
        RUNTIME_ERROR("unable to get client_id");
    }

    CRegClient* p_client = MWAServer.RegClients.FindClient(client_id);
    if(p_client == NULL) {
        RUNTIME_ERROR("unable to find client");
    }

    MWAServer.MWAAccumulator.GetInitialData(Command);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
