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

#include <string.h>
#include <ErrorSystem.hpp>
#include <ServerCommand.hpp>
#include "MTDServerHist.hpp"
#include "MTDRegClient.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDServerHist::RegisterOrCheckCoords(CServerCommand* p_cmd)
{
    CXMLElement* p_cele = p_cmd->GetRootCommandElement();
    CXMLElement* p_cvs = p_cele->GetFirstChildElement("CVS");
    if(p_cvs == NULL) {
        LOGIC_ERROR("unable to open CVS element");
    }

    bool result = true;

    HistMutex.Lock();

        if(GetNumOfCVs() == 0) {
            try {
                // load coordinates from first client
                LoadCVSInfo(p_cvs);
            } catch(std::exception& e) {
                CMD_ERROR_FROM_EXCEPTION(p_cmd,"unable to load coordinates",e);
                result = false;
            }
        } else {
            // coodinates have been loaded - so check if other clients are the same
            result = CheckCVSInfo(p_cvs);
            if( result == false ){
                // copy last reason to command error stack
                CMD_ERROR(p_cmd,ErrorSystem.GetLastError());
            }
        }

    HistMutex.Unlock();

    return(result);
}

//------------------------------------------------------------------------------

void CMTDServerHist::GetData(CServerCommand* p_command)
{
    if(p_command == NULL) {
        INVALID_ARGUMENT("p_command is NULL");
    }

    CXMLElement* p_result = p_command->GetRootResultElement();

    HistMutex.Lock();

        try {
            SaveCVSInfo(p_result);
            if(GetNumOfCVs() > 0) {
                WriteMTDData(p_result);
            }
        } catch(...){
            HistMutex.Unlock();
            throw;
        }

    HistMutex.Unlock();
}

//------------------------------------------------------------------------------

CMTDBuffer* CMTDServerHist::AddBuffer(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CMTDBuffer* p_buffer;

    HistMutex.Lock();
        try {
            p_buffer = AddMTDDataAsSingleBuffer(p_ele);
        } catch(...) {
            HistMutex.Unlock();
            throw;
        }
    HistMutex.Unlock();

    return(p_buffer);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDServerHist::FlushData(const CSmallString& name)
{
    HistMutex.Lock();
        try {
            Save(name);
        } catch(...) {
            HistMutex.Unlock();
            throw;
        }
    HistMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDServerHist::RegisterAllBuffersAsNew(CMTDRegClient* p_rc)
{
    if(p_rc == NULL) {
        INVALID_ARGUMENT("p_rc is NULL");
    }

    HistMutex.Lock();

        try {
            CSimpleIterator<CMTDBuffer>     I(Buffers);
            CMTDBuffer* p_buffer;

            while((p_buffer = I.Current()) != NULL) {
                p_rc->RegisterBuffer(p_buffer);
                I++;
            }
        } catch(...) {
            HistMutex.Unlock();
            throw;
        }

    HistMutex.Unlock();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

