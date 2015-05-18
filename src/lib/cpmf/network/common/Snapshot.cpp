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

#include <Snapshot.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSnapshot::CSnapshot(void)
{
    BoxA = 0;
    BoxB = 0;
    BoxC = 0;
    BoxAlpha = 0;
    BoxBeta = 0;
    BoxGamma = 0;
}

//------------------------------------------------------------------------------

CSnapshot::~CSnapshot(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSnapshot::SetNumOfAtoms(unsigned int numofatoms)
{
    Crds.CreateVector(numofatoms);
    Vels.CreateVector(numofatoms);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CSnapshot::Load(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        ES_ERROR("p_ele is NULL");
        return(false);
    }

    // load coordinates
    CXMLBinData* p_crd = p_ele->GetFirstChildBinData("CRDS");
    if(p_crd == NULL) {
        ES_ERROR("unable to find CRDS element");
        return(false);
    }
    Crds.Load(p_crd);

    // load velocities
    CXMLBinData* p_vel = p_ele->GetFirstChildBinData("VELS");
    if(p_vel == NULL) {
        ES_ERROR("unable to find VELS element");
        return(false);
    }
    Vels.Load(p_vel);

    // load box data
    CXMLElement* p_box = p_ele->GetFirstChildElement("BOX");
    if(p_vel == NULL) {
        ES_ERROR("unable to find BOX element");
        return(false);
    }

    bool result = true;
    result &= p_box->GetAttribute("a",BoxA);
    result &= p_box->GetAttribute("b",BoxB);
    result &= p_box->GetAttribute("c",BoxC);
    result &= p_box->GetAttribute("alpha",BoxAlpha);
    result &= p_box->GetAttribute("beta",BoxBeta);
    result &= p_box->GetAttribute("gamma",BoxGamma);

    if(result == false) {
        ES_ERROR("unable to get some box parameter");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CSnapshot::Save(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    // save coordinates
    CXMLBinData* p_crd = p_ele->CreateChildBinData("CRDS");
    Crds.Save(p_crd);

    // save velocities
    CXMLBinData* p_vel = p_ele->CreateChildBinData("VELS");
    Vels.Save(p_vel);

    // save box data
    CXMLElement* p_box = p_ele->CreateChildElement("BOX");
    p_box->SetAttribute("a",BoxA);
    p_box->SetAttribute("b",BoxB);
    p_box->SetAttribute("c",BoxC);
    p_box->SetAttribute("alpha",BoxAlpha);
    p_box->SetAttribute("beta",BoxBeta);
    p_box->SetAttribute("gamma",BoxGamma);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

