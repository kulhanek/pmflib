// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
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
#include <MTDClient.hpp>
#include <PMFMainHeader.hpp>
#include <SimpleVector.hpp>

extern "C" {

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// set number items in cpmf part

void PMF_PACKAGE cpmf_mtd_client_set_header_(FTINT* ret_st,FTINT* nitems)
{
    int l_nitems = *nitems;

    try {
        MTDClient.SetNumOfCVs(l_nitems);
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

void PMF_PACKAGE cpmf_mtd_client_set_coord_(FTINT* ret_st,
                                FTINT* id,
                                char* type,
                                char* name,
                                double* min_value,
                                double* max_value,
                                FTINT* nbins,
                                UFTINT type_len,
                                UFTINT name_len
                                )
{
    int            l_id = *id - 1;
    CSmallString   l_type;
    CSmallString   l_name;
    int            l_nbins = *nbins;
    double         l_min_value = *min_value;
    double         l_max_value = *max_value;

    try {
        l_type.SetFromFortran(type,type_len);
        l_name.SetFromFortran(name,name_len);

        MTDClient.SetCV(l_id,l_name,l_type,l_min_value,l_max_value,l_nbins);

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

void PMF_PACKAGE cpmf_mtd_client_reg_by_name_(char* fserver,char* fpassword,
                                  FTINT* client_id,
                                  UFTINT fserver_len,
                                  UFTINT fpassword_len)
{
// setup info about server
    CSmallString   server_name;
    CSmallString   server_password;

    try {

        // server name and protocol
        server_name.SetFromFortran(fserver,fserver_len);
        MTDClient.ActionRequest.SetProtocolName("mtd");
        MTDClient.ActionRequest.SetQualifiedName(server_name);

        // password
        server_password.SetFromFortran(fpassword,fpassword_len);
        MTDClient.ActionRequest.SetPassword(server_password);

    } catch(...) {
        *client_id = -1;
        return;
    }

// register client
    *client_id = MTDClient.RegisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_reg_by_key_(char* fserverkey,char* fserver,
                                 FTINT* client_id,
                                 UFTINT fserverkey_len,
                                 UFTINT fserver_len)
{
// setup info about server
    CSmallString   serverkey_name;
    CSmallString   server_name;

    try {

        // server name and protocol
        serverkey_name.SetFromFortran(fserverkey,fserverkey_len);
        MTDClient.ActionRequest.SetProtocolName("mtd");
        MTDClient.ActionRequest.ReadServerKey(serverkey_name);

        // get server name
        server_name = MTDClient.ActionRequest.GetQualifiedName();
        server_name.ConvertToFortran(fserver,fserver_len);
    } catch(...) {
        *client_id = -1;
        return;
    }

    // register client
    *client_id = MTDClient.RegisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_initial_data_(FTINT* ret_st)
{
    if(MTDClient.GetInitialData() == false) {
        *ret_st = 1;
        return;
    }
    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_exchange_data_(FTINT* ret_st)
{
    if(MTDClient.ExchangeData() == false) {
        *ret_st = 1;
        return;
    }
    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_get_buffer_info_(FTINT* ret_st,
                                      FTINT* id,
                                      FTINT* level,
                                      FTINT* start,
                                      FTINT* num_of_hills)
{
    *level = 0;
    *start = 0;
    *num_of_hills = 0;

// is there a buffer with index id?
    CMTDBuffer* p_buffer = MTDClient.GetBuffer(*id-1);
    if(p_buffer == NULL) {
        *ret_st = 1;
        return;
    }

    *level = p_buffer->GetLevel();
    *start = p_buffer->GetStart();
    *num_of_hills = p_buffer->GetNumOfHills();

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_add_buffer_data_(FTINT* ret_st,
                                      FTINT* level,
                                      FTINT* start,
                                      FTINT* num_of_hills,
                                      double* data)
{
    int l_num_of_hills = *num_of_hills;

// create new buffer
    CMTDBuffer* p_buffer = MTDClient.GetNewBuffer(l_num_of_hills);
    if(p_buffer == NULL) {
        *ret_st = 1;
        return;
    }

    p_buffer->SetLevel(*level);
    p_buffer->SetStart(*start);

// set data
    double* src = data;

    for(int i=0; i < l_num_of_hills; i++) {
        p_buffer->SetHeight(i,*src++);
        for(int j=0; j < p_buffer->GetNumOfCVs(); j++) {
            p_buffer->SetValue(i,j,*src++);
            p_buffer->SetWidth(i,j,*src++);
        }
        p_buffer->IncNumberOfHills();
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_get_buffer_data_(FTINT* ret_st,
                                      FTINT* id,
                                      double* data)
{
// is there a buffer with index id?
    CMTDBuffer* p_buffer = MTDClient.GetBuffer(*id-1);
    if(p_buffer == NULL) {
        *ret_st = 1;
        return;
    }

// copy all data
    double* dst = data;

// copy data
    for(int i=0; i < p_buffer->GetNumOfHills(); i++) {
        *dst++ = p_buffer->GetHeight(i);
        for(int j=0; j < p_buffer->GetNumOfCVs(); j++) {
            *dst++ = p_buffer->GetValue(i,j);
            *dst++ = p_buffer->GetWidth(i,j);
        }
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_clear_buf_list_(void)
{
    MTDClient.Deallocate();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_unregister_(void)
{
    MTDClient.UnregisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

}
