// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// ===============================================================================

#include <stdio.h>
#include <ErrorSystem.hpp>
#include "MWAAdmin.hpp"
#include <PMFOperation.hpp>
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CMWAAdmin::GetPMFAccumulator(void)
{
    CClientCommand  cmd;
    CPMFAccumulator accu;

    try{

        // init command
        InitCommand(&cmd,OperationPMF_GetData);

        // execute command
        ExecuteCommand(&cmd);

        // process response
        CXMLElement* p_accu = cmd.GetResultElementByPath("PMFLIB-V6",false);
        if(p_accu == NULL) {
            LOGIC_ERROR("unable to open PMFLIB-V6 element");
        }
        accu.Load(p_accu);

    } catch(std::exception& e) {
        ES_ERROR_FROM_EXCEPTION("unable to process command",e);
        return(false);
    }

    vout << endl;
    vout << "::::::::::::::::::::::::::::::::::: Output data ::::::::::::::::::::::::::::::::" << endl;

    CSmallString file_output;

    if(ActionRequest.GetParameterKeyValue("file",file_output) == false) {
        file_output = "_mwaserver.rst";
    }

// and now save all data
    if(accu.GetNumOfCVs() > 0) {
        vout << "Output PMF accumulator: " << file_output << endl;
        try {
            accu.Save(file_output);
        } catch(...) {
            ES_ERROR("unable to save PMF output accumulator");
        }
    } else {
        vout << "Output PMF accumulator: " << file_output << endl;
        vout << ">>> INFO: No data in PMF accumulator." << endl;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

