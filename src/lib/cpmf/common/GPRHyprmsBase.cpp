// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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
// =============================================================================

#include <GPRHyprmsBase.hpp>


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGPRHyprmsBase::CGPRHyprmsBase(void)
{
}

//------------------------------------------------------------------------------

CGPRHyprmsBase::~CGPRHyprmsBase(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGPRHyprmsBase::GetLogML(void)
{
    double ml = 0.0;
    return(ml);
}

//------------------------------------------------------------------------------

double CGPRHyprmsBase::GetLogPL(void)
{
    double loo = 0.0;
    return(loo);
}

//------------------------------------------------------------------------------

void CGPRHyprmsBase::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
}

//------------------------------------------------------------------------------

void CGPRHyprmsBase::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
}

//------------------------------------------------------------------------------
