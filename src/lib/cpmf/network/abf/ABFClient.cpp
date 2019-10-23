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

#include <ABFClient.hpp>
#include <PMFOperation.hpp>
#include <ClientCommand.hpp>
#include <ErrorSystem.hpp>
#include <ColVariable.hpp>
#include <string.h>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <ExtraOperation.hpp>
#include <cstdlib>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFClient      ABFClient;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFClient::CABFClient(void)
{
    ClientID = -1;
    NItems   = 0;
    NTotBins = 0;
}

//------------------------------------------------------------------------------

CABFClient::~CABFClient(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFClient::SetNumberOfItems(int nitems,int ntotbins)
{
    if(nitems <= 0) {
        INVALID_ARGUMENT("number of items has to be grater than zero");
    }

    if(ntotbins <= 0) {
        INVALID_ARGUMENT("number of bins has to be grater than zero");
    }

    NItems = 0;
    Coords.resize(nitems);
    NItems = nitems;
    NTotBins = ntotbins;
}

//------------------------------------------------------------------------------

int CABFClient::GetNumberOfBins(void)
{
    return(NTotBins);
}

//------------------------------------------------------------------------------

void CABFClient::SetCoord(int id,const CSmallString& name,const CSmallString& type,
                          double min_value,double max_value,int nbins)
{
    if((id < 0) || (id >= NItems)) {
        INVALID_ARGUMENT("cv id is out of range");
    }

    Coords[id].SetCoord(id,name,type,min_value,max_value,nbins);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFClient::RegisterClient(void)
{
    ClientID = -1;

    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_RegisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        CSmallString job_id = getenv("INF_JOB_ID");
        if( job_id != NULL ) {
            p_ele->SetAttribute("job_id",job_id);
        }
        p_ele = cmd.GetCommandElementByPath("CVS",true);
        SaveCVSInfo(p_ele);

        // execute command
        ExecuteCommand(&cmd);

        // process response
        p_ele = cmd.GetRootResultElement();
        if( p_ele->GetAttribute("client_id",ClientID) == false ) {
            LOGIC_ERROR("unable to get client_id");
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(-1);
    }

    return(ClientID);
}

//------------------------------------------------------------------------------

bool CABFClient::UnregisterClient(void)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_UnregisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);

        // execute command
        ExecuteCommand(&cmd);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFClient::GetInitialData(int* nisamples,
                                double* iabfforce,
                                double* iabfforce2)
{
    ClearExchangeData(nisamples,iabfforce,iabfforce2);

    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_GetInitialData);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);

        // execute command
        ExecuteCommand(&cmd);

        // process results
        p_ele = cmd.GetRootResultElement();
        ReadExchangeData(p_ele,nisamples,iabfforce,iabfforce2);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFClient::ExchangeData(int* nisamples,
                              double* iabfforce,
                              double* iabfforce2)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_ExchangeData);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);
        WriteExchangeData(p_ele,nisamples,iabfforce,iabfforce2);

        // execute command
        ExecuteCommand(&cmd);

        // process results
        p_ele = cmd.GetRootResultElement();
        ReadExchangeData(p_ele,nisamples,iabfforce,iabfforce2);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFClient::SaveCVSInfo(CXMLElement* p_tele)
{
    if(p_tele == NULL) {
        INVALID_ARGUMENT("p_tele is NULL");
    }

    CXMLElement* p_ele = p_tele->CreateChildElement("CVS");

    p_ele->SetAttribute("NCoords",NItems);

    for(int i=0; i < NItems; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("COORD");
        Coords[i].SaveInfo(p_iele);
    }
}


//------------------------------------------------------------------------------

void CABFClient::WriteExchangeData(CXMLElement* p_cele,
                                   int* nisamples,
                                   double* iabfforce,
                                   double* iabfforce2)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_cele is NULL");
    }

    int nisamples_size = NTotBins*sizeof(int);
    int iabfforce_size = NTotBins*NItems*sizeof(double);

    if((nisamples_size == 0) || (iabfforce_size == 0)) {
        LOGIC_ERROR("size(s) is(are) zero");
    }

// write new data -------------------------------
    CXMLBinData* p_nisamples = p_cele->CreateChildBinData("NISAMPLES");
    p_nisamples->SetData(nisamples,nisamples_size);

    CXMLBinData* p_iabfforce = p_cele->CreateChildBinData("IABFFORCE");
    p_iabfforce->SetData(iabfforce,iabfforce_size);

    CXMLBinData* p_iabfforce2 = p_cele->CreateChildBinData("IABFFORCE2");
    p_iabfforce2->SetData(iabfforce2,iabfforce_size);
}

//------------------------------------------------------------------------------

void CABFClient::ReadExchangeData(CXMLElement* p_rele,
                                  int* nisamples,
                                  double* iabfforce,
                                  double* iabfforce2)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    unsigned int nisamples_size = NTotBins*sizeof(int);
    unsigned int iabfforce_size = NTotBins*NItems*sizeof(double);

    if((nisamples_size == 0) || (iabfforce_size == 0)) {
        LOGIC_ERROR("size(s) is(are) zero");
    }

    CXMLBinData* p_nisamples = p_rele->GetFirstChildBinData("NISAMPLES");
    if(p_nisamples == NULL) {
        LOGIC_ERROR("unable to open NISAMPLES element");
    }

    void* p_data = p_nisamples->GetData();
    if((p_data == NULL) || (p_nisamples->GetLength() != nisamples_size)) {
        LOGIC_ERROR("inconsistent NISAMPLES element data");
    }
    memcpy(nisamples,p_data,nisamples_size);


    CXMLBinData* p_iabfforce = p_rele->GetFirstChildBinData("IABFFORCE");
    if(p_iabfforce == NULL) {
        LOGIC_ERROR("unable to open IABFFORCE element");
    }

    p_data = p_iabfforce->GetData();
    if((p_data == NULL) || (p_iabfforce->GetLength() != iabfforce_size)) {
        LOGIC_ERROR("inconsistent IABFFORCE element data");
    }
    memcpy(iabfforce,p_data,iabfforce_size);


    CXMLBinData* p_iabfforce2 = p_rele->GetFirstChildBinData("IABFFORCE2");
    if(p_iabfforce2 == NULL) {
        LOGIC_ERROR("unable to open IABFFORCE2 element");
    }

    p_data = p_iabfforce2->GetData();
    if((p_data == NULL) || (p_iabfforce2->GetLength() != iabfforce_size)) {
        LOGIC_ERROR("inconsistent IABFFORCE2 element data");
    }
    memcpy(iabfforce2,p_data,iabfforce_size);
}

//------------------------------------------------------------------------------

void CABFClient::ClearExchangeData(int* nisamples,
                                   double* iabfforce,
                                   double* iabfforce2)
{
    int nisamples_size = NTotBins;
    int iabfforce_size = NTotBins*NItems;

    for(int i=0; i < nisamples_size; i++) {
        nisamples[i] = 0;
    }

    for(int i=0; i < iabfforce_size; i++) {
        iabfforce[i] = 0.0;
        iabfforce2[i] = 0.0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

