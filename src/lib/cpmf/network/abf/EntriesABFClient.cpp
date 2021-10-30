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

void PMF_PACKAGE cpmf_abf_client_set_header_(FTINT* ret_st,
                                 FTINT*     ncvs,
                                 FTINT*     nbins,
                                 char*      version,
                                 char*      driver,
                                 double*    temp,
                                 char*      temp_unit,
                                 double*    temp_fconv,
                                 char*      ene_unit,
                                 double*    ene_fconv,
                                 FTINT*     enthalpy_enabled,
                                 FTINT*     mwa_mode,
                                 UFTINT     version_len,
                                 UFTINT     driver_len,
                                 UFTINT     temp_unit_len,
                                 UFTINT     ene_unit_len
                                 )
{
    int             l_ncvs          = *ncvs;
    int             l_nbins         = *nbins;
    double          l_temp          = *temp;
    double          l_temp_fconv    = *temp_fconv;
    double          l_ene_fconv     = *ene_fconv;
    int             l_mwa_mode      = *mwa_mode;
    CSmallString    l_version;
    CSmallString    l_driver;
    CSmallString    l_temp_unit;
    CSmallString    l_ene_unit;

    bool            l_enthalpy_enabled = false;
    if( *enthalpy_enabled > 0 ) l_enthalpy_enabled = true;

    try {
        l_version.SetFromFortran(version,version_len);
        l_driver.SetFromFortran(driver,driver_len);
        l_temp_unit.SetFromFortran(temp_unit,temp_unit_len);
        l_ene_unit.SetFromFortran(ene_unit,ene_unit_len);

        ABFClient.NumOfCVs = l_ncvs;
        ABFClient.NumOfBins = l_nbins;
        ABFClient.EnthalpyEnabled = l_enthalpy_enabled;
        ABFClient.MWAMode = l_mwa_mode;

        ABFClient.Accu->SetNumOfCVs(l_ncvs);
        ABFClient.Accu->SetHeaders("ABF",l_version,l_driver,l_temp,l_temp_unit,l_temp_fconv,l_ene_unit,l_ene_fconv);

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
                                FTINT*  id,
                                char*   name,
                                char*   type,
                                double* min_value,
                                double* max_value,
                                FTINT*  nbins,
                                double* fconv,
                                char*   unit,
                                UFTINT  name_len,
                                UFTINT  type_len,
                                UFTINT  unit_len
)
{
    int            l_id = *id - 1;
    CSmallString   l_name;
    CSmallString   l_type;
    CSmallString   l_unit;
    double         l_min_value = *min_value;
    double         l_max_value = *max_value;
    double         l_fconv = *fconv;
    int            l_nbins = *nbins;

    try {

        l_name.SetFromFortran(name,name_len);
        l_type.SetFromFortran(type,type_len);
        l_unit.SetFromFortran(unit,unit_len);

        ABFClient.Accu->SetCV(l_id,l_name,l_type,l_min_value,l_max_value,l_nbins,l_fconv,l_unit);

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

void PMF_PACKAGE cpmf_abf_client_reg_by_key_(char* fserverkey,char* fserver,
                                 FTINT* client_id,
                                 UFTINT fserverkey_len,
                                 UFTINT fserver_len)
{
// setup info about server
    CSmallString   serverkey_name;
    CSmallString   server_name;

    try {

        // finalize accumulator after CVs update
        ABFClient.Accu->UpdateNumOfBins();

        // server name and protocol
        serverkey_name.SetFromFortran(fserverkey,fserverkey_len);
        ABFClient.ActionRequest.SetProtocolName("mwa");
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

void PMF_PACKAGE cpmf_abf_client_initial_data_( FTINT* ret_st,
                                                double* inc_nsamples,
                                                double* inc_micf,
                                                double* inc_m2icf,
                                                double* inc_mepot,
                                                double* inc_m2epot)
{
    if(ABFClient.GetInitialData(inc_nsamples,inc_micf,inc_m2icf,inc_mepot,inc_m2epot) == false) {
        ES_ERROR("unable to get initial data - core");
        *ret_st = 1;
        return;
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_abf_client_exchange_data_(FTINT* ret_st,
                                                double* inc_nsamples,
                                                double* inc_micf,
                                                double* inc_m2icf,
                                                double* inc_mepot,
                                                double* inc_m2epot)
{
    if(ABFClient.ExchangeData(inc_nsamples,inc_micf,inc_m2icf,
                              inc_mepot,inc_m2epot) == false) {
        ES_ERROR("unable to exchange data");
        *ret_st = 1;
        return;
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
