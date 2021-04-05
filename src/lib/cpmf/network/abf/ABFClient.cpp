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
    NCVs   = 0;
    NTotBins = 0;
    EpotEnabled = false;
    Temperature = 300.0;
    EnergyFConv = 1.0;
    EnergyUnit = "kcal mol-1";
}

//------------------------------------------------------------------------------

CABFClient::~CABFClient(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFClient::SetNumberOfItems(int ncvs,int ntotbins)
{
    if(ncvs <= 0) {
        INVALID_ARGUMENT("number of CVs has to be grater than zero");
    }

    if(ntotbins <= 0) {
        INVALID_ARGUMENT("number of bins has to be grater than zero");
    }

    NCVs = ncvs;
    CVs.resize(ncvs);
    NTotBins = ntotbins;
}

//------------------------------------------------------------------------------

void CABFClient::SetTemperature(double temp)
{
    Temperature = temp;
}

//------------------------------------------------------------------------------

void CABFClient::SetEnergyUnit(double fconv,const CSmallString& unit)
{
    EnergyFConv = fconv;
    EnergyUnit = unit;
}

//------------------------------------------------------------------------------

void CABFClient::SetEpotEnabled(bool enabled)
{
    EpotEnabled = enabled;
}

//------------------------------------------------------------------------------

int CABFClient::GetNumOfBins(void)
{
    return(NTotBins);
}

//------------------------------------------------------------------------------

void CABFClient::SetCV(int id,const CSmallString& name,const CSmallString& type,
                          double min_value,double max_value,int nbins,double fconv,const CSmallString& unit)
{
    if((id < 0) || (id >= NCVs)) {
        INVALID_ARGUMENT("cv id is out of range");
    }

    CVs[id].SetCV(id,name,type,min_value,max_value,nbins,fconv,unit);
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
        p_ele->SetAttribute("temperature",Temperature);
        p_ele->SetAttribute("epot",EpotEnabled);

        p_ele = cmd.GetCommandElementByPath("ENERGY",true);
        p_ele->SetAttribute("fconv",EnergyFConv);
        p_ele->SetAttribute("unit",EnergyUnit);

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
                                double* inc_icfsum,
                                double* inc_icfsum2,
                                double* inc_epotsum,
                                double* inc_epotsum2,
                                double* inc_icfepotsum,
                                double* inc_icfepotsum2)
{
    ClearExchangeData(nisamples,inc_icfsum,inc_icfsum2,inc_epotsum,inc_epotsum2,inc_icfepotsum,inc_icfepotsum2);

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
        ReadExchangeData(p_ele,nisamples,inc_icfsum,inc_icfsum2,inc_epotsum,inc_epotsum2,inc_icfepotsum,inc_icfepotsum2);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFClient::ExchangeData(int* nisamples,
                                double* inc_icfsum,
                                double* inc_icfsum2,
                                double* inc_epotsum,
                                double* inc_epotsum2,
                                double* inc_icfepotsum,
                                double* inc_icfepotsum2)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_ExchangeData);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);
        WriteExchangeData(p_ele,nisamples,inc_icfsum,inc_icfsum2,inc_epotsum,inc_epotsum2,inc_icfepotsum,inc_icfepotsum2);

        // execute command
        ExecuteCommand(&cmd);

        // process results
        p_ele = cmd.GetRootResultElement();
        ReadExchangeData(p_ele,nisamples,inc_icfsum,inc_icfsum2,inc_epotsum,inc_epotsum2,inc_icfepotsum,inc_icfepotsum2);

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

    p_ele->SetAttribute("ncvs",NCVs);

    for(int i=0; i < NCVs; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("CV");
        CVs[i].SaveInfo(p_iele);
    }
}

//------------------------------------------------------------------------------

void CABFClient::WriteExchangeData(CXMLElement* p_cele,
                                    int*    inc_nsamples,
                                    double* inc_icfsum,
                                    double* inc_icfsum2,
                                    double* inc_epotsum,
                                    double* inc_epotsum2,
                                    double* inc_icfepotsum,
                                    double* inc_icfepotsum2)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_cele is NULL");
    }

    size_t nsamples_size   = NTotBins*sizeof(int);
    size_t icfsum_size     = NTotBins*NCVs*sizeof(double);
    size_t epotsum_size    = NTotBins*sizeof(double);
    size_t icfepotsum_size = NTotBins*NCVs*sizeof(double);

    if( (nsamples_size == 0) || (icfsum_size == 0) || (epotsum_size == 0) || (icfepotsum_size == 0) ) {
        LOGIC_ERROR("size(s) is(are) zero");
    }

// write new data -------------------------------
    CXMLBinData* p_nisamples = p_cele->CreateChildBinData("NSAMPLES");
    p_nisamples->SetData(inc_nsamples,nsamples_size);

    CXMLBinData* p_inc_icfsum = p_cele->CreateChildBinData("ICFSUM");
    p_inc_icfsum->SetData(inc_icfsum,icfsum_size);

    CXMLBinData* p_inc_icfsum2 = p_cele->CreateChildBinData("ICFSUM2");
    p_inc_icfsum2->SetData(inc_icfsum2,icfsum_size);

    CXMLBinData* p_inc_epotsum = p_cele->CreateChildBinData("EPOTSUM");
    p_inc_epotsum->SetData(inc_epotsum,epotsum_size);

    CXMLBinData* p_inc_epotsum2 = p_cele->CreateChildBinData("EPOTSUM2");
    p_inc_epotsum2->SetData(inc_epotsum2,epotsum_size);

    CXMLBinData* p_inc_icfepotsum = p_cele->CreateChildBinData("ICFEPOTSUM");
    p_inc_icfepotsum->SetData(inc_icfepotsum,icfepotsum_size);

    CXMLBinData* p_inc_icfepotsum2 = p_cele->CreateChildBinData("ICFEPOTSUM2");
    p_inc_icfepotsum2->SetData(inc_icfepotsum2,icfepotsum_size);
}

//------------------------------------------------------------------------------

void CABFClient::ReadExchangeData(CXMLElement* p_rele,
                                    int*    inc_nsamples,
                                    double* inc_icfsum,
                                    double* inc_icfsum2,
                                    double* inc_epotsum,
                                    double* inc_epotsum2,
                                    double* inc_icfepotsum,
                                    double* inc_icfepotsum2)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    size_t nsamples_size   = NTotBins*sizeof(int);
    size_t icfsum_size     = NTotBins*NCVs*sizeof(double);
    size_t epotsum_size    = NTotBins*sizeof(double);
    size_t icfepotsum_size = NTotBins*NCVs*sizeof(double);

    if( (nsamples_size == 0) || (icfsum_size == 0) || (epotsum_size == 0) || (icfepotsum_size == 0) ) {
        LOGIC_ERROR("size(s) is(are) zero");
    }

// --------------------------
    CXMLBinData* p_inc_nsamples = p_rele->GetFirstChildBinData("NSAMPLES");
    if(p_inc_nsamples == NULL) {
        LOGIC_ERROR("unable to open NSAMPLES element");
    }

    void* p_data = p_inc_nsamples->GetData();
    if((p_data == NULL) || (p_inc_nsamples->GetLength() != nsamples_size)) {
        LOGIC_ERROR("inconsistent NISAMPLES element data");
    }
    memcpy(inc_nsamples,p_data,nsamples_size);

// --------------------------
    CXMLBinData* p_inc_icfsum = p_rele->GetFirstChildBinData("ICFSUM");
    if(p_inc_icfsum == NULL) {
        LOGIC_ERROR("unable to open ICFSUM element");
    }

    p_data = p_inc_icfsum->GetData();
    if((p_data == NULL) || (p_inc_icfsum->GetLength() != icfsum_size)) {
        LOGIC_ERROR("inconsistent ICFSUM element data");
    }
    memcpy(inc_icfsum,p_data,icfsum_size);

// --------------------------
    CXMLBinData* p_inc_icfsum2 = p_rele->GetFirstChildBinData("ICFSUM2");
    if(p_inc_icfsum2 == NULL) {
        LOGIC_ERROR("unable to open ICFSUM2 element");
    }

    p_data = p_inc_icfsum2->GetData();
    if((p_data == NULL) || (p_inc_icfsum2->GetLength() != icfsum_size)) {
        LOGIC_ERROR("inconsistent ICFSUM2 element data");
    }
    memcpy(inc_icfsum2,p_data,icfsum_size);

// --------------------------
    CXMLBinData* p_inc_epotsum = p_rele->GetFirstChildBinData("EPOTSUM");
    if(p_inc_epotsum == NULL) {
        LOGIC_ERROR("unable to open EPOTSUM element");
    }

    p_data = p_inc_epotsum->GetData();
    if((p_data == NULL) || (p_inc_epotsum->GetLength() != epotsum_size)) {
        LOGIC_ERROR("inconsistent EPOTSUM element data");
    }
    memcpy(inc_epotsum,p_data,epotsum_size);

// --------------------------
    CXMLBinData* p_inc_epotsum2 = p_rele->GetFirstChildBinData("EPOTSUM2");
    if(p_inc_epotsum2 == NULL) {
        LOGIC_ERROR("unable to open EPOTSUM2 element");
    }

    p_data = p_inc_epotsum2->GetData();
    if((p_data == NULL) || (p_inc_epotsum2->GetLength() != epotsum_size)) {
        LOGIC_ERROR("inconsistent EPOTSUM2 element data");
    }
    memcpy(inc_epotsum2,p_data,epotsum_size);

// --------------------------
    CXMLBinData* p_inc_icfepotsum = p_rele->GetFirstChildBinData("ICFEPOTSUM");
    if(p_inc_icfepotsum == NULL) {
        LOGIC_ERROR("unable to open ICFEPOTSUM element");
    }

    p_data = p_inc_icfepotsum->GetData();
    if((p_data == NULL) || (p_inc_icfepotsum->GetLength() != icfepotsum_size)) {
        LOGIC_ERROR("inconsistent ICFEPOTSUM element data");
    }
    memcpy(inc_icfepotsum,p_data,icfepotsum_size);

// --------------------------
    CXMLBinData* p_inc_icfepotsum2 = p_rele->GetFirstChildBinData("ICFEPOTSUM2");
    if(p_inc_icfepotsum2 == NULL) {
        LOGIC_ERROR("unable to open ICFEPOTSUM2 element");
    }

    p_data = p_inc_icfepotsum2->GetData();
    if((p_data == NULL) || (p_inc_icfepotsum2->GetLength() != icfepotsum_size)) {
        LOGIC_ERROR("inconsistent ICFEPOTSUM2 element data");
    }
    memcpy(inc_icfepotsum2,p_data,icfepotsum_size);
}

//------------------------------------------------------------------------------

void CABFClient::ClearExchangeData(int* nisamples,
                                    double* inc_icfsum,
                                    double* inc_icfsum2,
                                    double* inc_epotsum,
                                    double* inc_epotsum2,
                                    double* inc_icfepotsum,
                                    double* inc_icfepotsum2)
{
    size_t nsamples_size   = NTotBins*sizeof(int);
    size_t icfsum_size     = NTotBins*NCVs*sizeof(double);
    size_t epotsum_size    = NTotBins*sizeof(double);
    size_t icfepotsum_size = NTotBins*NCVs*sizeof(double);

    for(size_t i=0; i < nsamples_size; i++) {
        nisamples[i] = 0;
    }

    for(size_t i=0; i < icfsum_size; i++) {
        inc_icfsum[i] = 0.0;
        inc_icfsum2[i] = 0.0;
    }

    for(size_t i=0; i < epotsum_size; i++) {
        inc_epotsum[i] = 0.0;
        inc_epotsum2[i] = 0.0;
    }

    for(size_t i=0; i < icfepotsum_size; i++) {
        inc_icfepotsum[i] = 0.0;
        inc_icfepotsum2[i] = 0.0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

