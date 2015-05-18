// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2009 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <PMFResidue.hpp>
#include <string.h>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFResidue::CPMFResidue(void)
{
    Index = -1;
    FirstAtomIndex = -1;
    NumOfAtoms = 0;
}

//------------------------------------------------------------------------------

CPMFResidue::~CPMFResidue(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFResidue::GetIndex(void) const
{
    return(Index);
}

//------------------------------------------------------------------------------

const char* CPMFResidue::GetName(void) const
{
    return(Name);
}

//------------------------------------------------------------------------------

int CPMFResidue::GetFirstAtomIndex(void) const
{
    return(FirstAtomIndex);
}

//------------------------------------------------------------------------------

int CPMFResidue::GetNumberOfAtoms(void) const
{
    return(NumOfAtoms);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
