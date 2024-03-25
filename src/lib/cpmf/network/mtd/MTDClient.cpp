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

#include <MTDClient.hpp>
#include <PMFOperation.hpp>
#include <ClientCommand.hpp>
#include <ErrorSystem.hpp>
#include <ColVariable.hpp>
#include <string.h>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <XMLPrinter.hpp>
#include <ExtraOperation.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDClient      MTDClient;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDClient::CMTDClient(void)
{
    ClientID        = -1;
    Accu            = CPMFAccumulatorPtr(new CPMFAccumulator);
    NumOfCVs        = 0;
    NumOfBins       = 0;
    WTemp           = 0.0;
}

//------------------------------------------------------------------------------

CMTDClient::~CMTDClient(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDClient::RegisterClient(void)
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

        Accu->UpdateNumOfBins();

        // create sections
        Accu->CreateSectionData("NSAMPLES", "AD","I","B");
        Accu->CreateSectionData("MTDPOT",   "AD","R","B");
        Accu->CreateSectionData("MTDFORCE", "AD","R","M");

        // set widths
        Accu->CreateSectionData("WIDTHS","SA","R","C");
        for(int i=0; i < NumOfCVs; i++){
           Accu->SetData("WIDTHS",i,Widths[i]);
        }

        if( WTemp > 0.0 ){
            Accu->CreateSectionData("MTD-WT","SA","R","D",1);
            Accu->SetData("MTD-WT",0,WTemp);
        }

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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMTDClient::UnregisterClient(void)
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

bool CMTDClient::GetInitialData(double* nsamples,
                                double* mtdforce,
                                double* mtdpot)
{
    ClearExchangeData(nsamples,mtdforce,mtdpot);

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

            if( accu->HasSectionData("MTDFORCE") ){
                data = accu->GetSectionData("MTDFORCE");
                data->GetDataBlob(mtdforce);
            }

            if( accu->HasSectionData("MTDPOT") ){
                data = accu->GetSectionData("MTDPOT");
                data->GetDataBlob(mtdpot);
            }
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CMTDClient::ExchangeData(  double* inc_nsamples,
                                double* inc_mtdforce,
                                double* inc_mtdpot)
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

        data = Accu->GetSectionData("MTDFORCE");
        data->SetDataBlob(inc_mtdforce);

        data = Accu->GetSectionData("MTDPOT");
        data->SetDataBlob(inc_mtdpot);

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

            if( accu->HasSectionData("MTDFORCE") ){
                data = accu->GetSectionData("MTDFORCE");
                data->GetDataBlob(inc_mtdforce);
            }

            if( accu->HasSectionData("MTDPOT") ){
                data = accu->GetSectionData("MTDPOT");
                data->GetDataBlob(inc_mtdpot);
            }
        }


    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CMTDClient::ClearExchangeData( double* inc_nsamples,
                                    double* inc_mtdforce,
                                    double* inc_mtdpot)
{
    int inc_nsamples_size   = NumOfBins;
    int inc_mtdforce_size   = NumOfBins*NumOfCVs;
    int inc_mtdpot_size     = NumOfBins;

    for(int i=0; i < inc_nsamples_size; i++) {
        inc_nsamples[i] = 0;
    }

    for(int i=0; i < inc_mtdforce_size; i++) {
        inc_mtdforce[i] = 0.0;
    }

    for(int i=0; i < inc_mtdpot_size; i++) {
        inc_mtdpot[i] = 0.0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

