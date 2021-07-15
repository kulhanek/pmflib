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
    EntropyEnabled  = false;
    NumOfCVs        = 0;
    NumOfBins       = 0;
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

        // create sections
        Accu->CreateSectionData("NSAMPLES", "AD","I","B");
        Accu->CreateSectionData("MICF",     "WA","R","M");
        Accu->CreateSectionData("M2ICF",    "AD","R","M");

        if( EnthalpyEnabled ){
            Accu->CreateSectionData("MEPOT",    "WA","R","B");
            Accu->CreateSectionData("M2EPOT",   "AD","R","B");
        }

        if( EnthalpyEnabled ){
            Accu->CreateSectionData("METOT",    "WA","R","B");
            Accu->CreateSectionData("M2ETOT",   "AD","R","B");
            Accu->CreateSectionData("C11HH",    "AD","R","M");
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
                                double* m2epot,
                                double* metot,
                                double* m2etot,
                                double* c11hh)
{
    ClearExchangeData(nsamples,micf,m2icf,mepot,m2epot,metot,m2etot,c11hh);

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

            if( EntropyEnabled ) {
                if( accu->HasSectionData("METOT") ){
                    data = accu->GetSectionData("METOT");
                    data->GetDataBlob(metot);
                }

                if( accu->HasSectionData("M2ETOT") ){
                    data = accu->GetSectionData("M2ETOT");
                    data->GetDataBlob(m2etot);
                }

                if( accu->HasSectionData("C11HH") ){
                    data = accu->GetSectionData("C11HH");
                    data->GetDataBlob(c11hh);
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
                                double* inc_m2epot,
                                double* inc_metot,
                                double* inc_m2etot,
                                double* inc_c11hh)
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

        // enthalpy
        if( EnthalpyEnabled ){
            data = Accu->GetSectionData("MEPOT");
            data->SetDataBlob(inc_mepot);

            data = Accu->GetSectionData("M2EPOT");
            data->SetDataBlob(inc_m2epot);
        }

        // entropy
        if( EntropyEnabled ){
            data = Accu->GetSectionData("METOT");
            data->SetDataBlob(inc_metot);

            data = Accu->GetSectionData("M2ETOT");
            data->SetDataBlob(inc_m2etot);

            data = Accu->GetSectionData("C11HH");
            data->SetDataBlob(inc_c11hh);
        }

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

            if( accu->HasSectionData("MICF") ){
                data = accu->GetSectionData("MICF");
                data->GetDataBlob(inc_micf);
            }

            if( accu->HasSectionData("M2ICF") ){
                data = accu->GetSectionData("M2ICF");
                data->GetDataBlob(inc_m2icf);
            }

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

            if( EntropyEnabled ) {
                if( accu->HasSectionData("METOT") ){
                    data = accu->GetSectionData("METOT");
                    data->GetDataBlob(inc_metot);
                }

                if( accu->HasSectionData("M2ETOT") ){
                    data = accu->GetSectionData("M2ETOT");
                    data->GetDataBlob(inc_m2etot);
                }

                if( accu->HasSectionData("C11HH") ){
                    data = accu->GetSectionData("C11HH");
                    data->GetDataBlob(inc_c11hh);
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
                                    double* inc_m2epot,
                                    double* inc_metot,
                                    double* inc_m2etot,
                                    double* inc_c11hh)
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

// ----------------------------------------------

    if( EntropyEnabled ) {
        size_t metot_size      = NumOfBins;
        size_t c11hh_size      = NumOfBins*NumOfCVs;

        for(size_t i=0; i < metot_size; i++) {
            inc_metot[i]  = 0.0;
            inc_m2etot[i] = 0.0;
        }
        for(size_t i=0; i < c11hh_size; i++) {
            inc_c11hh[i]  = 0.0;
        }
    }

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

