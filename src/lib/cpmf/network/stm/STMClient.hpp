#ifndef STMClientH
#define STMClientH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//      Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//      Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <SimpleVector.hpp>

//------------------------------------------------------------------------------

class CColVariable;
class CXMLElement;

//------------------------------------------------------------------------------

/*! \brief STM client for fpmf library
 *
 */

class PMF_PACKAGE CSTMClient : public CClient {
public:
    CSTMClient(void);
    ~CSTMClient(void);

// access methods -------------------------------------------------------------
    /// set number of cvs
    void SetNumberOfItems(int nitems);

    /// set cv
    void SetCoord(int id,const CSmallString& name,const CSmallString& type);

// commands -------------------------------------------------------------------
    /// register client on server side
    bool RegisterClient(int& cid,int& bid);

    /// unregister client on server side
    bool UnregisterClient(void);

    /// exchange data with server
    bool ExchangeData(int& mode,int& isteps,double* bpos,double* rpmf,double* rfz);

// section of private data ----------------------------------------------------
private:
    int                             ClientID;       // client ID
    int                             BeadID;         // bead ID
    int                             NumOfCVs;       // number of CVs
    CSimpleVector<CColVariable>     CVs;            // list of CVs

    /// save cvs into XML
    void SaveCVSInfo(CXMLElement* p_tele);
};

//------------------------------------------------------------------------------

extern CSTMClient      STMClient;

//------------------------------------------------------------------------------

#endif
