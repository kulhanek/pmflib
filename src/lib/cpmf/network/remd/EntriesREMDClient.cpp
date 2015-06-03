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

#include <stdio.h>
#include <SmallString.hpp>
#include <ErrorSystem.hpp>
#include <REMDClient.hpp>

extern "C" {

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_remd_client_reg_by_name_(int* natoms,char* fserver,char* fpassword,
                                   int* replica_id,
                                   unsigned int fserver_len,
                                   unsigned int fpassword_len)
{
// setup info about server
    CSmallString   server_name;
    CSmallString   server_password;

    try {

        // server name and protocol
        server_name.SetFromFortran(fserver,fserver_len);
        REMDClient.ActionRequest.SetProtocolName("remd");
        REMDClient.ActionRequest.SetQualifiedName(server_name);

        // password
        server_password.SetFromFortran(fpassword,fpassword_len);
        REMDClient.ActionRequest.SetPassword(server_password);

    } catch(...) {
        *replica_id = -1;
        return;
    }

// register client
    *replica_id = REMDClient.RegisterClient(*natoms);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_remd_client_reg_by_key_(int* natoms,char* fserverkey,char* fserver,
                                  int* replica_id,
                                  unsigned int fserverkey_len,
                                  unsigned int fserver_len)
{
    // setup info about server
        CSmallString   serverkey_name;
        CSmallString   server_name;

        try {

            // server name and protocol
            serverkey_name.SetFromFortran(fserverkey,fserverkey_len);
            REMDClient.ActionRequest.SetProtocolName("remd");
            REMDClient.ActionRequest.ReadServerKey(serverkey_name);

            // get server name
            server_name = REMDClient.ActionRequest.GetQualifiedName();
            server_name.ConvertToFortran(fserver,fserver_len);

        } catch(...) {
            *replica_id = -1;
            return;
        }

    // register client
        *replica_id = REMDClient.RegisterClient(*natoms);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_remd_client_initial_data_(int* ret_st,
                                    int* mode,
                                    int* period,
                                    int* bath_id,
                                    double* temp)
{

    if(REMDClient.GetInitialData(*mode,*period,*bath_id,*temp) == false) {
        fprintf(stderr,"\n>>> ERROR: Unable to get initial data!\n");
        *ret_st = 1;
        return;
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_remd_client_set_data_(int*    ret_st,
                                double* crds,
                                double* vels,
                                double* box_abc,
                                double* box_angles)
{
    REMDClient.SetSnapshotData(crds,vels,box_abc,box_angles);
    *ret_st = 0;
}

//------------------------------------------------------------------------------

void cpmf_remd_client_get_data_(int*    ret_st,
                                double* crds,
                                double* vels,
                                double* box_abc,
                                double* box_angles)
{
    REMDClient.GetSnapshotData(crds,vels,box_abc,box_angles);
    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_remd_client_exchange_data_(int*   ret_st,
                                     double* iepot,
                                     int*    bath_id,
                                     double* ctemp,
                                     double* otemp)
{
    if(REMDClient.ExchangeData(*iepot,*bath_id,*ctemp,*otemp) == false) {
        fprintf(stderr,"\n>>> ERROR: Unable to exchange data!\n");
        *ret_st = 1;
        return;
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_remd_client_unregister_(void)
{
    REMDClient.UnregisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

}
