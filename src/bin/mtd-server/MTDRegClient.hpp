#ifndef MTDRegClientH
#define MTDRegClientH
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

#include <RegClient.hpp>
#include <MTDBuffer.hpp>
#include <SimpleMutex.hpp>
#include <SimpleList.hpp>

//------------------------------------------------------------------------------

class CServerCommand;

//------------------------------------------------------------------------------

/*! \brief MTD registered client
 *
 */

class CMTDRegClient : public CRegClient {
public:
    CMTDRegClient(const CSmallString& client_name,const CSmallString& job_id);
    virtual ~CMTDRegClient(void);

    /// register buffer for client
    void RegisterBuffer(CMTDBuffer* p_buffer);

    /// send initial data to client
    void GetInitialData(CServerCommand* p_command);

    /// exchange data with client
    void ExchangeDataWithClient(CServerCommand* p_command);

// section of private data ----------------------------------------------------
private:
    CSimpleMutex                ClientMutex;
    CSimpleList<CMTDBuffer>     NewBuffers;
};

//------------------------------------------------------------------------------

#endif
