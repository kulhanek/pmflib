#ifndef PMFMaskH
#define PMFMaskH
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

#include <SmallString.hpp>
#include <PMFMainHeader.hpp>

//---------------------------------------------------------------------------

class CPMFTopology;
class CPMFAtom;
class CPMFMaskSelection;

//---------------------------------------------------------------------------

/// mask support class
/*!
 mask support is now fully compatible with AMBER 9.0
*/

class PMF_PACKAGE CPMFMask {
public:
    CPMFMask(void);
    ~CPMFMask(void);

    // topology assigment ------------------------------------------------------
    bool            AssignTopology(CPMFTopology* p_top);
    CPMFTopology*   GetTopology(void) const;

    // mask setup --------------------------------------------------------------
    /// select all atoms
    bool SelectAllAtoms(void);

    /// set mask from string specification
    bool SetMask(const CSmallString& mask);

    /// return current mask specification
    const CSmallString&     GetMask(void) const;

    // results -----------------------------------------------------------------
    int             GetNumberOfTopologyAtoms(void);
    int             GetNumberOfSelectedAtoms(void);
    CPMFAtom*       GetSelectedAtom(int index);
    bool            IsAtomSelected(int index);

// section of private data ----------------------------------------------------
private:
    CPMFTopology*           Topology;
    CSmallString            Mask;
    CPMFMaskSelection*      Selection;
};

//---------------------------------------------------------------------------

#endif
