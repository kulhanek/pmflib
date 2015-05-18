#ifndef MTDServerHistH
#define MTDServerHistH
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

#include <Server.hpp>
#include <MTDHistory.hpp>
#include <SimpleMutex.hpp>

//------------------------------------------------------------------------------

class CServerCommand;
class CMTDRegClient;

//------------------------------------------------------------------------------

class CMTDServerHist : public CMTDHistory {
public:

    /// send initial data to client
    bool RegisterOrCheckCoords(CServerCommand* p_cmd);

    /// send data to admin utility
    void GetData(CServerCommand* p_command);

    /// flush data
    void FlushData(const CSmallString& name);

    /// set new buffers for client
    void RegisterAllBuffersAsNew(CMTDRegClient* p_rc);

    /// add new buffer to history
    CMTDBuffer* AddBuffer(CXMLElement* p_ele);

// section of private data ----------------------------------------------------
private:
    CSimpleMutex        HistMutex;  // access mutex
};

//------------------------------------------------------------------------------

#endif
