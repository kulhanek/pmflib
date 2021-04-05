#ifndef REMDClientH
#define REMDClientH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301  USA
// =============================================================================

#include <PMFMainHeader.hpp>
#include <Client.hpp>
#include <Snapshot.hpp>

//------------------------------------------------------------------------------

class CRCoordinate;
class CXMLElement;

//------------------------------------------------------------------------------

/*! \brief REMD client for fpmf library
 *
 */

class PMF_PACKAGE CREMDClient : public CClient {
public:
    CREMDClient(void);
    ~CREMDClient(void);

// commands -------------------------------------------------------------------
    /// register client on server side
    int RegisterClient(int numofatoms);

    /// unregister client on server side
    bool UnregisterClient(void);

    /// get initial temperature
    bool GetInitialData(int& mode,int& period,int& bath_id,double& temp);

    /// set snapshot data
    void SetSnapshotData(double* crds,
                         double* vels,
                         double* box_abc,
                         double* box_angles);

    /// set snapshot data
    void GetSnapshotData(double* crds,
                         double* vels,
                         double* box_abc,
                         double* box_angles);

    /// exchange data with server
    bool ExchangeData(double inc_epotsum,int& bath_id,double& ctemp,double& otemp);

// section of private data ----------------------------------------------------
private:
    int             ReplicaID;      // replica ID
    int             BathID;         // bath ID
    int             Mode;           // mode - 0 - temp swap, 1 - phase point swap
    int             NumOfAtoms;     // finger print
    CSnapshot       Snapshot;       // system snapshot
};

//------------------------------------------------------------------------------

extern CREMDClient      REMDClient;

//------------------------------------------------------------------------------

#endif
