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

#include <stdio.h>
#include <SmallString.hpp>
#include <ErrorSystem.hpp>
#include <STMClient.hpp>

extern "C" {

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// set number items in cpmf part

void PMF_PACKAGE cpmf_stm_client_set_header_(int* ret_st,int* nitems)
{
    int l_nitems = *nitems;

    try {
        STMClient.SetNumberOfItems(l_nitems);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to set the number of items",e);
        *ret_st = 1;
        return;
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_stm_client_set_coord_(int* ret_st,
                                int* id,
                                char* type,
                                char* name,
                                unsigned int type_len,
                                unsigned int name_len)
{
    int            l_id = *id - 1;
    CSmallString   l_type;
    CSmallString   l_name;

    try {

        l_type.SetFromFortran(type,type_len);
        l_name.SetFromFortran(name,name_len);

        STMClient.SetCoord(l_id,l_name,l_type);

    } catch(std::exception& e) {
        CSmallString error;
        error << "unable to set the cv id: " << l_id;
        ES_ERROR_FROM_EXCEPTION(error,e);
        *ret_st = 1;
        return;
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_stm_client_reg_by_name_(char* fserver,char* fpassword,
                                  int* client_id,
                                  int* bead_id,
                                  unsigned int fserver_len,
                                  unsigned int fpassword_len)
{
// setup info about server
    CSmallString   server_name;
    CSmallString   server_password;

    try {

        // server name and protocol
        server_name.SetFromFortran(fserver,fserver_len);
        STMClient.ActionRequest.SetProtocolName("stm");
        STMClient.ActionRequest.SetQualifiedName(server_name);

        // password
        server_password.SetFromFortran(fpassword,fpassword_len);
        STMClient.ActionRequest.SetPassword(server_password);

    } catch(...) {
        *client_id = -1;
        *bead_id = -1;
        return;
    }

// register client
    STMClient.RegisterClient(*client_id,*bead_id);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_stm_client_reg_by_key_(char* fserverkey,char* fserver,
                                 int* client_id,
                                 int* bead_id,
                                 unsigned int fserverkey_len,
                                 unsigned int fserver_len)
{
// setup info about server
    CSmallString   serverkey_name;
    CSmallString   server_name;

    try {
        // server name and protocol
        serverkey_name.SetFromFortran(fserverkey,fserverkey_len);
        STMClient.ActionRequest.SetProtocolName("stm");
        STMClient.ActionRequest.ReadServerKey(serverkey_name);

        // get server name
        server_name = STMClient.ActionRequest.GetQualifiedName();
        server_name.ConvertToFortran(fserver,fserver_len);
    } catch(...) {
        // no client assigment
        *client_id = -1;
        *bead_id = -1;
        return;
    }

// register client
    STMClient.RegisterClient(*client_id,*bead_id);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

//integer         :: ret_st
//integer         :: mode
//integer         :: isteps
//real(8)         :: bpos(*)
//real(8)         :: rpmf(*)
//real(8)         :: rfz(*)

void PMF_PACKAGE cpmf_stm_client_exchange_data_(int* ret_st,
                                    int* mode,
                                    int* isteps,
                                    double* bpos,
                                    double* rpmf,
                                    double* rfz)
{

    if(STMClient.ExchangeData(*mode,*isteps,bpos,rpmf,rfz) == false) {
        ES_ERROR("unable to exchange data");
        *ret_st = 1;
        return;
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_stm_client_unregister_(void)
{
    STMClient.UnregisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

}
