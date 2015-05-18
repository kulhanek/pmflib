#ifndef SnapshotH
#define SnapshotH
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

#include <PMFMainHeader.hpp>
#include <SimpleVector.hpp>

//------------------------------------------------------------------------------

class CXMLElement;

//------------------------------------------------------------------------------

class PMF_PACKAGE CSnapshot {
public:
    CSnapshot(void);
    ~CSnapshot(void);

// setup methods ---------------------------------------------------------------
    /// set number of atoms
    void  SetNumOfAtoms(unsigned int numofatoms);

// i/o methods -----------------------------------------------------------------
    /// load snapshot
    bool Load(CXMLElement* p_ele);

    /// save snapshot
    void Save(CXMLElement* p_ele);

// section of private data ----------------------------------------------------
public:
    // system snapshot
    CSimpleVector<double>   Crds;                       // coordinates
    CSimpleVector<double>   Vels;                       // velocities
    double                  BoxA,BoxB,BoxC;             // box dimmensions
    double                  BoxAlpha,BoxBeta,BoxGamma;  // box dimmensions
};

//------------------------------------------------------------------------------

extern CSnapshot      Snapshot;

//------------------------------------------------------------------------------

#endif
