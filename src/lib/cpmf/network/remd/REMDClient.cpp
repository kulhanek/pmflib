// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <REMDClient.hpp>
#include <PMFOperation.hpp>
#include <ClientCommand.hpp>
#include <ErrorSystem.hpp>
#include <string.h>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <ExtraOperation.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CREMDClient      REMDClient;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CREMDClient::CREMDClient(void)
{
    ReplicaID = -1;
    BathID = -1;
    Mode = -1;
    NumOfAtoms = 0;
}

//------------------------------------------------------------------------------

CREMDClient::~CREMDClient(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CREMDClient::RegisterClient(int numofatoms)
{
    ReplicaID = -1;

    CClientCommand cmd;
    try{
        // prepare snapshot
        Snapshot.SetNumOfAtoms(numofatoms);

        // init command
        InitCommand(&cmd,Operation_RegisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        CSmallString job_id = getenv("INF_JOB_ID");
        if( job_id != NULL ) {
            p_ele->SetAttribute("job_id",job_id);
        }
        p_ele->SetAttribute("natoms",NumOfAtoms);

        // execute command
        ExecuteCommand(&cmd);

        // process response
        p_ele = cmd.GetRootResultElement();
        if( p_ele->GetAttribute("replica_id",ReplicaID) == false ) {
            LOGIC_ERROR("unable to get replica_id");
        }

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(-1);
    }

    return(ReplicaID);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CREMDClient::UnregisterClient(void)
{
    CClientCommand cmd;
    try{

        // init command
        InitCommand(&cmd,Operation_UnregisterClient);

        // prepare input data
        CXMLElement* p_ele = cmd.GetRootCommandElement();
        p_ele->SetAttribute("replica_id",ReplicaID);

        // execute command
        ExecuteCommand(&cmd);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CREMDClient::GetInitialData(int& mode,int& period,int& bath_id,double& temp)
{
//// create command
//    CClientCommand* p_command = NULL;
//    try{
//        p_command = CreateCommand(Operation_GetInitialData);
//    } catch(...) {
//        ES_ERROR("unable to create command");
//        return(false);
//    }

//// set replica ID
//    try{
//        CXMLElement* p_ele = p_command->GetRootCommandElement();
//        p_ele->SetAttribute("replica_id",ReplicaID);
//    } catch(...){
//        ES_ERROR("unable to prepare input data");
//        delete p_command;
//        return(false);
//    }

//// execute command
//    bool result = ExecuteCommand(p_command);
//    if( result == false ) {
//        ES_ERROR("unable to execute command");
//        delete p_command;
//        return(false);
//    }

//// read new data
//    try{
//        CXMLElement* p_ele = p_command->GetRootResultElement();
//        result = true;
//        result &= p_ele->GetAttribute("mode",mode);
//        result &= p_ele->GetAttribute("exch",period);
//        result &= p_ele->GetAttribute("bath_id",bath_id);
//        result &= p_ele->GetAttribute("temp",temp);

//        if(result == false) {
//            ES_ERROR("unable to get mode or temp attribute");
//            delete p_command;
//            return(false);
//        }

//        Mode = mode;
//        BathID = bath_id;

//    } catch(...){
//        ES_ERROR("unable to read initial data");
//        delete p_command;
//        return(false);
//    }

//// release command
//    delete p_command;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CREMDClient::SetSnapshotData(double* crds,
                                  double* vels,
                                  double* box_abc,
                                  double* box_angles)
{
    memcpy(Snapshot.Crds,crds,NumOfAtoms*sizeof(double));
    memcpy(Snapshot.Vels,vels,NumOfAtoms*sizeof(double));
    Snapshot.BoxA = box_abc[0];
    Snapshot.BoxB = box_abc[1];
    Snapshot.BoxC = box_abc[2];
    Snapshot.BoxAlpha = box_angles[0];
    Snapshot.BoxBeta  = box_angles[1];
    Snapshot.BoxGamma = box_angles[2];
}

//------------------------------------------------------------------------------

void CREMDClient::GetSnapshotData(double* crds,
                                  double* vels,
                                  double* box_abc,
                                  double* box_angles)
{
    memcpy(crds,Snapshot.Crds,NumOfAtoms*sizeof(double));
    memcpy(vels,Snapshot.Vels,NumOfAtoms*sizeof(double));
    box_abc[0] = Snapshot.BoxA;
    box_abc[1] = Snapshot.BoxB;
    box_abc[2] = Snapshot.BoxC;
    box_angles[0] = Snapshot.BoxAlpha;
    box_angles[1] = Snapshot.BoxBeta;
    box_angles[2] = Snapshot.BoxGamma;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CREMDClient::ExchangeData(double  iepot,
                               int&    bath_id,
                               double& ctemp,
                               double& otemp)
{
    //double temp = ctemp;

//// create command
//    CClientCommand* p_command = NULL;
//    try{
//        p_command = CreateCommand(Operation_ExchangeData);
//    } catch(...) {
//        ES_ERROR("unable to create command");
//        return(false);
//    }

//// set data
//    try{
//        CXMLElement* p_ele = p_command->GetRootCommandElement();
//        p_ele->SetAttribute("replica_id",ReplicaID);
//        p_ele->SetAttribute("bath_id",BathID);
//        p_ele->SetAttribute("epot",iepot);
//    } catch(...){
//        ES_ERROR("unable to prepare input data");
//        delete p_command;
//        return(false);
//    }

//// execute command
//    bool result = ExecuteCommand(p_command);
//    if( result == false ) {
//        ES_ERROR("unable to execute command");
//        delete p_command;
//        return(false);
//    }

//// read new data
//    try{
//        CXMLElement* p_ele = p_command->GetRootResultElement();
//        result = true;
//        result &= p_ele->GetAttribute("bath_id",BathID);
//        result &= p_ele->GetAttribute("temp",temp);

//        if(result == false) {
//            ES_ERROR("unable to get attributes");
//            delete p_command;
//            return(false);
//        }

//    } catch(...){
//        ES_ERROR("unable to read initial data");
//        delete p_command;
//        return(false);
//    }

//// release command
//    delete p_command;

//    bath_id = BathID;
//    switch(Mode) {
//        case 0:
//            otemp = ctemp;
//            ctemp  = temp;
//            break;
//        case 1:
//            otemp  = temp;
//            break;
//        default:
//            return(false);
//            break;
//    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

