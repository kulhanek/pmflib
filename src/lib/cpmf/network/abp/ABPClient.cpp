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
    ClientID        = -1;
    Accu            = CPMFAccumulatorPtr(new CPMFAccumulator);
    NumOfCVs        = 0;
    NumOfBins       = 0;
}

//------------------------------------------------------------------------------

CABPClient::~CABPClient(void)
{
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

        // create sections
        Accu->CreateSectionData("NSAMPLES", "AD","I","B");
        Accu->CreateSectionData("MICF",     "WA","R","M");
        Accu->CreateSectionData("M2ICF",    "AD","R","M");

        // FIXME - CV widths

        Accu->Save(p_ele);

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

bool CABPClient::GetInitialData(double* nsamples,
                                double* dpop,
                                double* pop)
{
    ClearExchangeData(nsamples,dpop,pop);

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
        p_ele = cmd.GetResultElementByPath("PMFLIB-V6",false);
        if( p_ele != NULL ){

            CPMFAccumulatorPtr accu = CPMFAccumulatorPtr(new CPMFAccumulator);
            accu->Load(p_ele);

            // check dimensions
            if( (NumOfBins != accu->GetNumOfBins()) ||
                (NumOfCVs  != accu->GetNumOfCVs()) ) {
                RUNTIME_ERROR("size inconsistency");
            }

            if( Accu->CheckCVSInfo(accu) == false ){
                RUNTIME_ERROR("inconsistent CVs");
            }

            // core data
            CPMFAccuDataPtr data;
            if( accu->HasSectionData("NSAMPLES") ){
                data = accu->GetSectionData("NSAMPLES");
                data->GetDataBlob(nsamples);
            }

            if( accu->HasSectionData("DPOP") ){
                data = accu->GetSectionData("DPOP");
                data->GetDataBlob(dpop);
            }

            if( accu->HasSectionData("POP") ){
                data = accu->GetSectionData("POP");
                data->GetDataBlob(pop);
            }
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABPClient::ExchangeData(  double* inc_nsamples,
                                double* inc_dpop,
                                double* inc_pop)
{
    CClientCommand cmd;
    try{

    // init command
        InitCommand(&cmd,OperationPMF_ExchangeData);

    // prepare input data
        CXMLElement* p_ele;

        p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);

        CPMFAccuDataPtr data;
        data = Accu->GetSectionData("NSAMPLES");
        data->SetDataBlob(inc_nsamples);

        data = Accu->GetSectionData("DPOP");
        data->SetDataBlob(inc_dpop);

        data = Accu->GetSectionData("POP");
        data->SetDataBlob(inc_pop);

        Accu->Save(p_ele);

    // execute command
        ExecuteCommand(&cmd);

    // process results
        p_ele = cmd.GetResultElementByPath("PMFLIB-V6",false);
        if( p_ele != NULL ){

            CPMFAccumulatorPtr accu = CPMFAccumulatorPtr(new CPMFAccumulator);
            accu->Load(p_ele);

            if( (NumOfBins != accu->GetNumOfBins()) ||
                (NumOfCVs  != accu->GetNumOfCVs()) ) {
                RUNTIME_ERROR("size inconsistency");
            }

            if( Accu->CheckCVSInfo(accu) == false ){
                RUNTIME_ERROR("inconsistent CVs");
            }

            // core data
            CPMFAccuDataPtr data;
            if( accu->HasSectionData("NSAMPLES") ){
                data = accu->GetSectionData("NSAMPLES");
                data->GetDataBlob(inc_nsamples);
            }

            if( accu->HasSectionData("DPOP") ){
                data = accu->GetSectionData("DPOP");
                data->GetDataBlob(inc_dpop);
            }

            if( accu->HasSectionData("POP") ){
                data = accu->GetSectionData("POP");
                data->GetDataBlob(inc_pop);
            }
        }


    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABPClient::ClearExchangeData( double* inc_nsamples,
                                    double* inc_dpop,
                                    double* inc_pop)
{
    int inc_nsamples_size = NumOfBins;
    int inc_dpop_size = NumOfBins*NumOfCVs;
    int inc_pop_size = NumOfBins;

    for(int i=0; i < inc_nsamples_size; i++) {
        inc_nsamples[i] = 0;
    }

    for(int i=0; i < inc_dpop_size; i++) {
        inc_dpop[i] = 0.0;
    }

    for(int i=0; i < inc_pop_size; i++) {
        inc_pop[i] = 0.0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

