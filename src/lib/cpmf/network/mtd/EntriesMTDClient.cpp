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

void PMF_PACKAGE cpmf_mtd_client_set_header_(FTINT* ret_st,
                                 FTINT*     ncvs,
                                 FTINT*     nbins,
                                 char*      version,
                                 char*      driver,
                                 double*    temp,
                                 char*      temp_unit,
                                 double*    temp_fconv,
                                 char*      ene_unit,
                                 double*    ene_fconv,
                                 double*    widths,
                                 double*    wt_temp,
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
    double          l_wt_temp       = *wt_temp;
    CSmallString    l_version;
    CSmallString    l_driver;
    CSmallString    l_temp_unit;
    CSmallString    l_ene_unit;

    try {
        l_version.SetFromFortran(version,version_len);
        l_driver.SetFromFortran(driver,driver_len);
        l_temp_unit.SetFromFortran(temp_unit,temp_unit_len);
        l_ene_unit.SetFromFortran(ene_unit,ene_unit_len);

        MTDClient.NumOfCVs = l_ncvs;
        MTDClient.NumOfBins = l_nbins;

        MTDClient.Accu->SetNumOfCVs(l_ncvs);
        MTDClient.Accu->SetHeaders("MTD",l_version,l_driver,l_temp,l_temp_unit,l_temp_fconv,l_ene_unit,l_ene_fconv);
        MTDClient.WTemp = l_wt_temp;

        MTDClient.Widths.CreateVector(l_ncvs);
        for(int i=0; i < l_ncvs; i++){
            MTDClient.Widths[i] = widths[i];
        }

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

        MTDClient.Accu->SetCV(l_id,l_name,l_type,l_min_value,l_max_value,l_nbins,l_fconv,l_unit);

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
        MTDClient.ActionRequest.SetProtocolName("mwa");
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

void PMF_PACKAGE cpmf_mtd_client_initial_data_(FTINT* ret_st,
                                    FTINT*  nsamples,
                                    double* mtdforce,
                                    double* mtdpot)
{
    CSimpleVector<double> l_nsamples;

    l_nsamples.CreateVector(MTDClient.NumOfBins);

    if(MTDClient.GetInitialData(l_nsamples,mtdforce,mtdpot) == false) {
        ES_ERROR("unable to get initial data");
        *ret_st = 1;
        return;
    }

    for(int i = 0; i < MTDClient.NumOfBins; i++){
        nsamples[i] = l_nsamples[i];
    }

    *ret_st = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE cpmf_mtd_client_exchange_data_(FTINT* ret_st,
                                    FTINT*  inc_nsamples,
                                    double* inc_mtdforce,
                                    double* inc_mtdpot)
{
    CSimpleVector<double> l_inc_nsamples;

    l_inc_nsamples.CreateVector(MTDClient.NumOfBins);

    for(int i = 0; i < MTDClient.NumOfBins; i++){
        l_inc_nsamples[i] = inc_nsamples[i];
    }

    if(MTDClient.ExchangeData(l_inc_nsamples,inc_mtdforce,inc_mtdpot) == false) {
        ES_ERROR("unable to exchange data");
        *ret_st = 1;
        return;
    }

    for(int i = 0; i < MTDClient.NumOfBins; i++){
        inc_nsamples[i] = l_inc_nsamples[i];
    }

    *ret_st = 0;
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
