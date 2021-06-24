#ifndef MWAServerAccuH
#define MWAServerAccuH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <PMFAccumulator.hpp>
#include <SimpleMutex.hpp>
#include <ServerCommand.hpp>

//------------------------------------------------------------------------------

class CRegClient;

//------------------------------------------------------------------------------

/// total MWA accumulator for MWA approach

class CMWAServerAccu : public CPMFAccumulator {
public:
    /// constructor
    CMWAServerAccu(void);

    /// send initial data to client
    bool RegisterOrCheckCoords(CServerCommand* p_cmd);

    /// send initial data to client
    void GetInitialData(CServerCommand* p_command);

    /// exchange data with client
    void ExchangeDataWithClient(CRegClient* p_client);

    /// send data to admin utility
    void GetData(CServerCommand* p_command);

    /// flush data
    void FlushData(const CSmallString& name);

    /// setup synchronous mode, if numofclients > 0
    void SetupSynchronousMode(int numofclients);

    /// get number of clients in synchronous mode
    int GetNumberOfClients(void);

    /// server is terminated - unblock waiting beads
    void SetServerTerminated(void);

    /// set temperature and energy unit
    void SetMainHeader(CServerCommand* p_command);

// section of private data ----------------------------------------------------
private:
    CSimpleMutex    AccuMutex;          // access mutex

    // synchronous exchange
    bool                        SynchronousMode;        // should we wait for all clients?
    bool                        ServerTerminated;
    int                         NumOfClients;           // number of clients
    CSimpleMutex                RendezvousMutex;
    CSimpleCond                 RendezvousCond;
    CSimpleList<CServerCommand> Commands;               // waiting client data

    /// echange data in asynchronous mode
    void ExchangeDataAsynchronousMode(CRegClient* p_client);

    /// exchange data in synchronous mode
    void ExchangeDataSynchronousMode(CRegClient* p_client);
};

//------------------------------------------------------------------------------

#endif
