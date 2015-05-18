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

#include <PMFAtom.hpp>
#include <string.h>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAtom::CPMFAtom(void)
{
    Index = -1;
    Residue = NULL;
    Mass = 0.0;
}

//------------------------------------------------------------------------------

CPMFAtom::~CPMFAtom(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFAtom::GetIndex(void) const
{
    return(Index);
}

//------------------------------------------------------------------------------

const char* CPMFAtom::GetName(void) const
{
    return(Name);
}

//------------------------------------------------------------------------------

const char* CPMFAtom::GetType(void) const
{
    return(Type);
}

//------------------------------------------------------------------------------

double CPMFAtom::GetMass(void) const
{
    return(Mass);
}

//------------------------------------------------------------------------------

const CPoint& CPMFAtom::GetPosition(void) const
{
    return(Position);
}

//------------------------------------------------------------------------------

CPMFResidue* CPMFAtom::GetResidue(void) const
{
    return(Residue);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
