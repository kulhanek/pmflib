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
    ClientID        = -1;
    Accu            = CPMFAccumulatorPtr(new CPMFAccumulator);
    EnthalpyEnabled = false;
    NumOfCVs        = 0;
    NumOfBins       = 0;
    MWAMode         = 0;
}

//------------------------------------------------------------------------------

CABFClient::~CABFClient(void)
{
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

        Accu->UpdateNumOfBins();

        // create sections
        Accu->CreateSectionData("NSAMPLES", "AD","R","B");
        Accu->CreateSectionData("MICF",     "WA","R","M", "NSAMPLES");
        Accu->CreateSectionData("M2ICF",    "M2","R","M", "NSAMPLES","MICF");

        if( MWAMode == 0 ) {
            if( EnthalpyEnabled ){
                Accu->CreateSectionData("MEPOT",    "WA","R","B", "NSAMPLES");
                Accu->CreateSectionData("M2EPOT",   "M2","R","B", "NSAMPLES","MEPOT");
            }
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

bool CABFClient::GetInitialData(double* nsamples,
                                double* micf,
                                double* m2icf,
                                double* mepot,
                                double* m2epot)
{
    ClearExchangeData(nsamples,micf,m2icf,mepot,m2epot);

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

            if( accu->HasSectionData("MICF") ){
                data = accu->GetSectionData("MICF");
                data->GetDataBlob(micf);
            }

            if( accu->HasSectionData("M2ICF") ){
                data = accu->GetSectionData("M2ICF");
                data->GetDataBlob(m2icf);
            }

            if( EnthalpyEnabled ) {
                if( accu->HasSectionData("MEPOT") ){
                    data = accu->GetSectionData("MEPOT");
                    data->GetDataBlob(mepot);
                }

                if( accu->HasSectionData("M2EPOT") ){
                    data = accu->GetSectionData("M2EPOT");
                    data->GetDataBlob(m2epot);
                }
            }
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFClient::ExchangeData(  double* inc_nsamples,
                                double* inc_micf,
                                double* inc_m2icf,
                                double* inc_mepot,
                                double* inc_m2epot)
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

        data = Accu->GetSectionData("MICF");
        data->SetDataBlob(inc_micf);

        data = Accu->GetSectionData("M2ICF");
        data->SetDataBlob(inc_m2icf);

        if( MWAMode == 0 ) {
            // enthalpy
            if( EnthalpyEnabled ){
                data = Accu->GetSectionData("MEPOT");
                data->SetDataBlob(inc_mepot);

                data = Accu->GetSectionData("M2EPOT");
                data->SetDataBlob(inc_m2epot);
            }
        }

        Accu->Save(p_ele);

// execute command
        ExecuteCommand(&cmd);

        ClearExchangeData(inc_nsamples,inc_micf,inc_m2icf,inc_mepot,inc_m2epot);

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

            if( accu->HasSectionData("MICF") ){
                data = accu->GetSectionData("MICF");
                data->GetDataBlob(inc_micf);
            }

            if( accu->HasSectionData("M2ICF") ){
                data = accu->GetSectionData("M2ICF");
                data->GetDataBlob(inc_m2icf);
            }

            if( MWAMode == 0 ) {
                if( EnthalpyEnabled ) {
                    if( accu->HasSectionData("MEPOT") ){
                        data = accu->GetSectionData("MEPOT");
                        data->GetDataBlob(inc_mepot);
                    }

                    if( accu->HasSectionData("M2EPOT") ){
                        data = accu->GetSectionData("M2EPOT");
                        data->GetDataBlob(inc_m2epot);
                    }
                }
            }
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFClient::ClearExchangeData( double* inc_nsamples,
                                    double* inc_micf,
                                    double* inc_m2icf,
                                    double* inc_mepot,
                                    double* inc_m2epot)
{
    size_t nsamples_size   = NumOfBins;
    size_t micf_size       = NumOfBins*NumOfCVs;

    for(size_t i=0; i < nsamples_size; i++) {
        inc_nsamples[i] = 0.0;
    }

    for(size_t i=0; i < micf_size; i++) {
        inc_micf[i]  = 0.0;
        inc_m2icf[i] = 0.0;
    }

// ----------------------------------------------

    if( EnthalpyEnabled ){
        size_t mepot_size      = NumOfBins;

        for(size_t i=0; i < mepot_size; i++) {
            inc_mepot[i]  = 0.0;
            inc_m2epot[i] = 0.0;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

