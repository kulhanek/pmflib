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
#include <ABFClient.hpp>
#include <PMFMainHeader.hpp>
#include <SimpleVector.hpp>

extern "C" {

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// set number items in cpmf part

void PMF_PACKAGE cpmf_abf_client_set_header_(FTINT* ret_st,
                                 FTINT* nitems,
                                 FTINT* tot_nbins,
                                 char* energy_unit,
                                 double* energy_unit_fac,
                                 UFTINT energy_unit_len
                                 )
{
    int l_nitems = *nitems;
    int l_tot_nbins = *tot_nbins;

    try {
        ABFClient.SetNumberOfItems(l_nitems,l_tot_nbins);

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

void PMF_PACKAGE cpmf_abf_client_set_coord_(FTINT* ret_st,
                                FTINT* id,
                                char* name,
                                char* type,
                                double* min_value,
                                double* max_value,
                                FTINT* nbins,
                                UFTINT name_len,
                                UFTINT type_len
)
{
    int            l_id = *id - 1;
    CSmallString   l_name;
    CSmallString   l_type;
    double         l_min_value = *min_value;
    double         l_max_value = *max_value;
    int            l_nbins = *nbins;

    try {

        l_name.SetFromFortran(name,name_len);
        l_type.SetFromFortran(type,type_len);

        ABFClient.SetCoord(l_id,l_name,l_type,l_min_value,l_max_value,l_nbins);

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

void PMF_PACKAGE cpmf_abf_client_reg_by_name_(char* fserver,char* fpassword,
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
        ABFClient.ActionRequest.SetProtocolName("abf");
        ABFClient.ActionRequest.SetQualifiedName(server_name);

        // password
        server_password.SetFromFortran(fpassword,fpassword_len);
        ABFClient.ActionRequest.SetPassword(server_password);

    } catch(...) {
        *client_id = -1;
        return;
    }

// register client
    *client_id = ABFClient.RegisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_abf_client_reg_by_key_(char* fserverkey,char* fserver,
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
        ABFClient.ActionRequest.SetProtocolName("abf");
        ABFClient.ActionRequest.ReadServerKey(serverkey_name);

        // get server name
        server_name = ABFClient.ActionRequest.GetQualifiedName();
        server_name.ConvertToFortran(fserver,fserver_len);
    } catch(...) {
        *client_id = -1;
        return;
    }

// register client
    *client_id = ABFClient.RegisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_abf_client_initial_data_(FTINT* ret_st,
                                   FTINT* nisamples,
                                   double* iabfforce,
                                   double* iabfforce2)
{
    CSimpleVector<int> isamples;

    if(ABFClient.GetInitialData(isamples,iabfforce,iabfforce2) == false) {
        ES_ERROR("unable to get initial data");
        *ret_st = 1;
        return;
    }

    for(int i = 0; i < ABFClient.GetNumberOfBins(); i++){
        nisamples[i] = isamples[i];
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_abf_client_exchange_data_(FTINT* ret_st,
                                    FTINT* nisamples,
                                    double* iabfforce,
                                    double* iabfforce2)
{
    CSimpleVector<int> isamples;

    isamples.CreateVector(ABFClient.GetNumberOfBins());
    for(int i = 0; i < ABFClient.GetNumberOfBins(); i++){
        isamples[i] = nisamples[i];
    }

    if(ABFClient.ExchangeData(isamples,iabfforce,iabfforce2) == false) {
        ES_ERROR("unable to exchange data");
        *ret_st = 1;
        return;
    }

    for(int i = 0; i < ABFClient.GetNumberOfBins(); i++){
        nisamples[i] = isamples[i];
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_abf_client_unregister_(void)
{
    ABFClient.UnregisterClient();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

}
