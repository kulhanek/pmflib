// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//      Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//      Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <STMClient.hpp>
#include <PMFOperation.hpp>
#include <ClientCommand.hpp>
#include <ErrorSystem.hpp>
#include <ColVariable.hpp>
#include <string.h>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <ExtraOperation.hpp>
#include <cstdlib>
#include <Bead.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSTMClient      STMClient;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSTMClient::CSTMClient(void)
{
    ClientID = -1;
    BeadID = -1;
    NumOfCVs   = 0;
}

//------------------------------------------------------------------------------

CSTMClient::~CSTMClient(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMClient::SetNumberOfItems(int nitems)
{
    if(nitems <= 0) {
        INVALID_ARGUMENT("number of items has to be grater than zero");
    }

    CVs.CreateVector(nitems);
    NumOfCVs = nitems;
}

//------------------------------------------------------------------------------

void CSTMClient::SetCoord(int id,const CSmallString& name,const CSmallString& type)
{
    if((id < 0) || (id >= NumOfCVs)) {
        INVALID_ARGUMENT("cv id is out of range");
    }

    CVs[id].SetCoord(id,name,type);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSTMClient::RegisterClient(int& cid,int& bid)
{
    ClientID = -1;
    BeadID = bid;
    cid = ClientID;

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
        p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("bead_id",BeadID);

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
        bid = -1;
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    cid = ClientID;

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMClient::UnregisterClient(void)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_UnregisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);
        p_ele->SetAttribute("bead_id",BeadID);

        // execute command
        ExecuteCommand(&cmd);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CSTMClient::ExchangeData(int& mode,int& isteps,double* bpos,double* rpmf,double* rfz)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,OperationPMF_ExchangeData);

        long unsigned int bpos_size = NumOfCVs*sizeof(double);
        long unsigned int rpmf_size = NumOfCVs*sizeof(double);
        long unsigned int rfz_size  = NumOfCVs*NumOfCVs*sizeof(double);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("client_id",ClientID);
        p_ele->SetAttribute("bead_id",BeadID);
        p_ele->SetAttribute("mode",mode);

        if( (mode == BMO_ACCUMULATION) || (mode == BMO_PRODUCTION) ) {
            CXMLBinData* p_rpmfele = p_ele->CreateChildBinData("PMF");
            p_rpmfele->SetData(rpmf,rpmf_size,false,EXBDT_DOUBLE);

            CXMLBinData* p_rfzele = p_ele->CreateChildBinData("MTZ");
            p_rfzele->SetData(rfz,rfz_size,false,EXBDT_DOUBLE);
            p_rfzele->SetAttribute("rows",NumOfCVs);
            p_rfzele->SetAttribute("columns",NumOfCVs);
        }

        // execute command
        ExecuteCommand(&cmd);

        // process results
        p_ele = cmd.GetRootResultElement();

        p_ele->GetAttribute("mode",mode);
        p_ele->GetAttribute("steps",isteps);

        if( mode != BMO_TERMINATE ){
            CXMLBinData* p_bpos = p_ele->GetFirstChildBinData("BPOS");
            if(p_bpos == NULL) {
                LOGIC_ERROR("unable to open BPOS element");
            }

            void* p_data = p_bpos->GetData();
            if((p_data == NULL) || (p_bpos->GetLength() != bpos_size)) {
                LOGIC_ERROR("inconsistent BPOS element data");
            }
            memcpy(bpos,p_data,bpos_size);
        }
//        } else {
//            // do not change bead prosition!
//        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSTMClient::SaveCVSInfo(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_tele is NULL");
    }

    p_ele->SetAttribute("ncvs",NumOfCVs);

    for(int i=0; i < NumOfCVs; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("COORD");
        CVs[i].SaveInfo(p_iele);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

