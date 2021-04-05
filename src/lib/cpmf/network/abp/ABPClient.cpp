// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <ABPClient.hpp>
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

CABPClient      ABPClient;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABPClient::CABPClient(void)
{
    ClientID = -1;
    NCVs   = 0;
    NTotBins = 0;
}

//------------------------------------------------------------------------------

CABPClient::~CABPClient(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPClient::SetNumberOfItems(int NCVs,int ntotbins)
{
    if(NCVs <= 0) {
        INVALID_ARGUMENT("number of items has to be grater than zero");
    }

    if(ntotbins <= 0) {
        INVALID_ARGUMENT("number of bins has to be grater than zero");
    }

    NCVs = 0;
    CVs.resize(NCVs);
    NCVs = NCVs;
    NTotBins = ntotbins;
}

//------------------------------------------------------------------------------

int CABPClient::GetNumOfBins(void)
{
    return(NTotBins);
}

//------------------------------------------------------------------------------

void CABPClient::SetCV(int id,const CSmallString& name,const CSmallString& type,
                       double min_value,double max_value,int nbins,double alpha,double fconv,const CSmallString& unit)
{
    if((id < 0) || (id >= NCVs)) {
        INVALID_ARGUMENT("cv id is out of range");
    }

    CVs[id].SetCV(id,name,type,min_value,max_value,nbins,alpha,fconv,unit);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPClient::RegisterClient(void)
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

bool CABPClient::UnregisterClient(void)
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

bool CABPClient::GetInitialData(int* nisamples,
                                double* idpop,
                                double* ipop)
{
    ClearExchangeData(nisamples,idpop,ipop);

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
        ReadExchangeData(p_ele,nisamples,idpop,ipop);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABPClient::ExchangeData(int* nisamples,
                                double* idpop,
                                double* ipop)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_ExchangeData);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);
        WriteExchangeData(p_ele,nisamples,idpop,ipop);

        // execute command
        ExecuteCommand(&cmd);

        // process results
        p_ele = cmd.GetRootResultElement();
        ReadExchangeData(p_ele,nisamples,idpop,ipop);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPClient::SaveCVSInfo(CXMLElement* p_tele)
{
    if(p_tele == NULL) {
        INVALID_ARGUMENT("p_tele is NULL");
    }

    CXMLElement* p_ele = p_tele->CreateChildElement("CVS");

    p_ele->SetAttribute("NCVs",NCVs);

    for(int i=0; i < NCVs; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("COORD");
        CVs[i].SaveInfo(p_iele);
    }
}

//------------------------------------------------------------------------------

void CABPClient::WriteExchangeData(CXMLElement* p_cele,
                                    int* nisamples,
                                    double* idpop,
                                    double* ipop)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_cele is NULL");
    }

    int nisamples_size  = NTotBins*sizeof(int);
    int idpop_size      = NTotBins*NCVs*sizeof(double);
    int ipop_size       = NTotBins*sizeof(double);

    if( (nisamples_size == 0) || (idpop_size == 0) || (ipop_size == 0) ) {
        LOGIC_ERROR("size(s) is(are) zero");
    }

// write new data -------------------------------
    CXMLBinData* p_nisamples = p_cele->CreateChildBinData("NISAMPLES");
    p_nisamples->SetData(nisamples,nisamples_size);

    CXMLBinData* p_idpop = p_cele->CreateChildBinData("IDPOP");
    p_idpop->SetData(idpop,idpop_size);

    CXMLBinData* p_ipop = p_cele->CreateChildBinData("IPOP");
    p_ipop->SetData(ipop,ipop_size);
}

//------------------------------------------------------------------------------

void CABPClient::ReadExchangeData(CXMLElement* p_rele,
                                    int* nisamples,
                                    double* idpop,
                                    double* ipop)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    unsigned int nisamples_size = NTotBins*sizeof(int);
    unsigned int idpop_size     = NTotBins*NCVs*sizeof(double);
    unsigned int ipop_size      = NTotBins*sizeof(double);

    if( (nisamples_size == 0) || (idpop_size == 0) || (ipop_size == 0) ) {
        LOGIC_ERROR("size(s) is(are) zero");
    }

// --------------------------
    CXMLBinData* p_nisamples = p_rele->GetFirstChildBinData("NISAMPLES");
    if(p_nisamples == NULL) {
        LOGIC_ERROR("unable to open NISAMPLES element");
    }

    void* p_data = p_nisamples->GetData();
    if((p_data == NULL) || (p_nisamples->GetLength() != nisamples_size)) {
        LOGIC_ERROR("inconsistent NISAMPLES element data");
    }
    memcpy(nisamples,p_data,nisamples_size);

// --------------------------
    CXMLBinData* p_idpop= p_rele->GetFirstChildBinData("IDPOP");
    if(p_idpop == NULL) {
        LOGIC_ERROR("unable to open IDPOP element");
    }

    p_data = p_idpop->GetData();
    if((p_data == NULL) || (p_idpop->GetLength() != idpop_size)) {
        LOGIC_ERROR("inconsistent IDPOP element data");
    }
    memcpy(idpop,p_data,idpop_size);

// --------------------------
    CXMLBinData* p_ipop = p_rele->GetFirstChildBinData("IPOP");
    if(p_ipop == NULL) {
        LOGIC_ERROR("unable to open IPOP element");
    }

    p_data = p_ipop->GetData();
    if((p_data == NULL) || (p_ipop->GetLength() != ipop_size)) {
        LOGIC_ERROR("inconsistent IPOP element data");
    }
    memcpy(ipop,p_data,ipop_size);
}

//------------------------------------------------------------------------------

void CABPClient::ClearExchangeData(int* nisamples,
                                    double* idpop,
                                    double* ipop)
{
    int nisamples_size = NTotBins;
    int idpop_size = NTotBins*NCVs;
    int ipop_size = NTotBins;

    for(int i=0; i < nisamples_size; i++) {
        nisamples[i] = 0;
    }

    for(int i=0; i < idpop_size; i++) {
        idpop[i] = 0.0;
    }

    for(int i=0; i < ipop_size; i++) {
        ipop[i] = 0.0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

