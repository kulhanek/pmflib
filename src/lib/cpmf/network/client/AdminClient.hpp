#ifndef AdminClientH
#define AdminClientH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2006      Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <ExtraClient.hpp>
#include <iostream>

//------------------------------------------------------------------------------

class PMF_PACKAGE CAdminClient : protected CExtraClient {
    public:
        // constructor
        CAdminClient(void);

// supported operations -------------------------------------------------------
    /// get status information from server
    bool GetServerInfo(std::ostream& vout);

    /// flush all data to disk
    bool FlushServerData(void);

    /// server info - helper method
    void GetCVServerInfo(CClientCommand* p_command,std::ostream& sout);
};

//------------------------------------------------------------------------------

#endif
