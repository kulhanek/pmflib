#ifndef PMFResidueH
#define PMFResidueH
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

//------------------------------------------------------------------------------

/** \brief simple topology residue
*/

class PMF_PACKAGE CPMFResidue {
public:
// constructor and destructor -------------------------------------------------
    CPMFResidue(void);
    ~CPMFResidue(void);

// information methods --------------------------------------------------------
    /// get residue index - counted from zero
    int         GetIndex(void) const;

    /// get residue name
    const char* GetName(void) const;

    /// return index of the first atom
    /*! index starts from zero and has regular counting order
    */
    int         GetFirstAtomIndex(void) const;

    /// return number of atoms which belongs to residue
    int         GetNumberOfAtoms(void) const;

// section of private data ----------------------------------------------------
private:
    CSmallString    Name;               /// residue name
    int             Index;              /// residue index - counted from zero
    int             FirstAtomIndex;     /// index of first atom in residue - counted from zero
    int             NumOfAtoms;         /// number of atoms in residues

    friend class CPMFTopology;
};

//------------------------------------------------------------------------------

#endif
