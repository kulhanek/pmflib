// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
// ===============================================================================

#include <PMFMainHeader.hpp>
#include <ErrorSystem.hpp>
#include <PMFMainHeaderConfig.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const char* LibBuildVersion_PMF = PMF_VERSION " (" PMF_DATE ")";

//------------------------------------------------------------------------------

extern "C" {

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

    void PMF_PACKAGE cpmf_print_errors_(void)
    {
        ErrorSystem.PrintErrors();
    }

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

#if defined _WIN32 || defined __CYGWIN__

#include <windows.h>

BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
    return(TRUE);
}

#endif

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
