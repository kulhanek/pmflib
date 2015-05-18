// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
// ===============================================================================

#include <ErrorSystem.hpp>
#include <ServerCommand.hpp>
#include "MTDRegClient.hpp"
#include "MTDServer.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDRegClient::CMTDRegClient(const CSmallString& client_name,const CSmallString& job_id)
    : CRegClient(client_name,job_id)
{

}

//------------------------------------------------------------------------------

CMTDRegClient::~CMTDRegClient(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDRegClient::RegisterBuffer(CMTDBuffer* p_buffer)
{
    ClientMutex.Lock();
        try {
            NewBuffers.InsertToEnd(p_buffer);
        } catch(...){
            ClientMutex.Unlock();
            throw;
        }
    ClientMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDRegClient::GetInitialData(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_rele = p_command->GetRootResultElement();

    ClientMutex.Lock();

        try {
            // write new data
            MTDServer.MTDHistory.WriteMTDData(p_rele,NewBuffers);

            // remove registered buffers
            NewBuffers.RemoveAll();
        } catch(...) {
            ClientMutex.Unlock();
            throw;
        }

    ClientMutex.Unlock();
}

//------------------------------------------------------------------------------

void CMTDRegClient::ExchangeDataWithClient(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_cele = p_command->GetRootCommandElement();
    CXMLElement* p_rele = p_command->GetRootResultElement();

    ClientMutex.Lock();
        try {
            // add received data as new buffer to the history
            CMTDBuffer* p_buffer = MTDServer.MTDHistory.AddBuffer(p_cele);

            p_buffer->SetLevel(GetClientID());

            // now register this buffer for all other clients
            int numofclients = MTDServer.RegClients.GetNumberOfRegClients();

            for(int i=0; i < numofclients; i++) {
                CMTDRegClient* p_client = dynamic_cast<CMTDRegClient*>(MTDServer.RegClients.GetClient(i));
                if((p_client != NULL) && (p_client != this)) {
                    p_client->RegisterBuffer(p_buffer);
                }
            }

            // write new data from other clients
            MTDServer.MTDHistory.WriteMTDData(p_rele,NewBuffers);

            // remove registered buffers
            NewBuffers.RemoveAll();

        } catch(...) {
            ClientMutex.Unlock();
            throw;
        }

    ClientMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

