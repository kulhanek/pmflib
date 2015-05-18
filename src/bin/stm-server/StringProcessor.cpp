// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
// ===============================================================================

#include <stdio.h>
#include <ErrorSystem.hpp>
#include <PMFOperation.hpp>
#include <ExtraOperation.hpp>
#include "StringProcessor.hpp"
#include "StringServer.hpp"
#include <ServerCommand.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CStringProcessor::CStringProcessor(CServerCommand* p_cmd)
    : CCmdProcessor(p_cmd)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CStringProcessor::ProcessCommand(void)
{
    try {
        if(Operation == Operation_GetServerInfo) {
            GetServerInfo();
            return(true);
        }
        if(Operation == Operation_FlushServerData) {
            FlushServerData();
            return(true);
        }
        if(Operation == Operation_RegisterClient) {
            RegisterClient();
            return(true);
        }
        if(Operation == Operation_UnregisterClient) {
            UnregisterClient();
            return(true);
        }
        if(Operation == Operation_ExchangeData) {
            ExchangeData();
            return(true);
        }
        if(Operation == Operation_GetData) {
            GetData();
            return(true);
        }
        if(Operation == Operation_ShutdownServer) {
            ShutdownTerminateServer();
            return(true);
        }
    } catch(std::exception& e) {
        CMD_ERROR_FROM_EXCEPTION(Command,"unable process the request",e);
        return(false);
    }

    CSmallString error;
    error << "operation " << Operation.GetStringForm() << " is not implemented";
    CMD_ERROR(Command,error);

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
