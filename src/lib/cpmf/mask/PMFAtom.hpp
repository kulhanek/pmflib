#ifndef PMFAtomH
#define PMFAtomH
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
#include <Point.hpp>

//---------------------------------------------------------------------------

class CPMFResidue;

//---------------------------------------------------------------------------

/** \brief simple topology atom
*/

class PMF_PACKAGE CPMFAtom {
public:
// constructor and destructor -------------------------------------------------
    CPMFAtom(void);
    ~CPMFAtom(void);

// information methods --------------------------------------------------------
    /// get atom index
    int                     GetIndex(void) const;

    /// get atom name
    const char*     GetName(void) const;

    /// get atom type
    const char*     GetType(void) const;

    /// get atom mass
    double                  GetMass(void) const;

    /// get atom position
    const CPoint&           GetPosition(void) const;

    /// get residue
    CPMFResidue*            GetResidue(void) const;

// section of private data ----------------------------------------------------
private:
    int             Index;      // atom index
    CPMFResidue*    Residue;    // residue
    CSmallString    Name;       // atom name
    CSmallString    Type;       // atom type
    double          Mass;       // atom mass
    CPoint          Position;   // atom position

    friend class CPMFTopology;
};

//---------------------------------------------------------------------------

#endif
