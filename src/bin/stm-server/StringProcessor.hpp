#ifndef StringProcessorH
#define StringProcessorH
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

#include <CmdProcessor.hpp>

//------------------------------------------------------------------------------

class CStringProcessor : public CCmdProcessor {
public:
    // constructor
    CStringProcessor(CServerCommand* p_cmd);

// section of private data -----------------------------------------------------
private:
    virtual bool ProcessCommand(void);

// implemented operations ------------------------------------------------------
    void GetServerInfo(void);
    void FlushServerData(void);
    void RegisterClient(void);
    void UnregisterClient(void);
    void ExchangeData(void);
    void GetData(void);
    void ShutdownTerminateServer(void);
};

//------------------------------------------------------------------------------

#endif
