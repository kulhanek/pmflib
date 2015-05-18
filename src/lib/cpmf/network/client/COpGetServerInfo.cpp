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

#include <AdminClient.hpp>
#include <ExtraOperation.hpp>
#include <ErrorSystem.hpp>
#include <ClientCommand.hpp>
#include <SmallTime.hpp>
#include <RegClient.hpp>
#include <ColVariable.hpp>
#include <XMLElement.hpp>
#include <iostream>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAdminClient::GetServerInfo(std::ostream& vout)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_GetServerInfo);

        // execute command
        ExecuteCommand(&cmd);

        // print response
        GetCVServerInfo(&cmd,vout);
        GetShortServerInfo(&cmd,vout);
        GetLongServerInfo(&cmd,vout);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CAdminClient::GetCVServerInfo(CClientCommand* p_command,std::ostream& vout)
{
    if( p_command == NULL ){
        INVALID_ARGUMENT("p_command is NULL");
    }

    vout << endl;
    vout << "=== Collective Variables =======================================================" << endl;
    vout << endl;
    vout << "ID  Type      Name                          Min value       Max value     NBins " << endl;
    vout << "-- ---------- -------------------------- --------------- --------------- -------" << endl;

    // write cv summary if it is present
    CXMLElement* p_ele  = p_command->GetRootResultElement();
    CXMLElement* p_icel = p_ele->GetChildElementByPath("CVS/COORD");
    int          counter = 1;

    while( p_icel != NULL ){
        CColVariable cv;

        cv.LoadInfo(p_icel);
        cv.PrintInfo(vout);
        counter++;

        p_icel = p_icel->GetNextSiblingElement("COORD");
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

