#ifndef PMFLibCVsH
#define PMFLibCVsH
// =============================================================================
// CATS - Conversion and Analysis Tools
// -----------------------------------------------------------------------------
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <CATsMainHeader.hpp>
#include <PMFLibCV.hpp>
#include <SimpleList.hpp>
#include <AmberMaskAtoms.hpp>

class CSmallString;
class CPrmFile;
class CAmberRestart;

//------------------------------------------------------------------------------

class CPMFLibCVs : private CSimpleList<CPMFLibCV> {
public:
// constructor -----------------------------------------------------------------
    CPMFLibCVs(void);

    //! set topology
    bool SetTopology(CAmberTopology* p_top);

// methods ---------------------------------------------------------------------
    //! read CVs
    bool Read(const CSmallString& name);

    //! get CV value
    bool  GetValue(const CSmallString& name,CAmberRestart* p_crd,double& value);

// section of private data -----------------------------------------------------
private:
    int             fnatom; // total number of atoms
    CAmberMaskAtoms mask;

    //! create CV by type
    CPMFLibCV* CreateCV(const CSmallString& cv_type);

    //! get number of atoms in the mask
    int pmf_mask_natoms_in_mask(const CSmallString& smask);

    //! set mask
    bool pmf_mask_set_mask(const CSmallString& smask);

    //! is atom selected?
    bool pmf_mask_is_atom_selected(int i);

    friend class CPMFLibCV;
};

//------------------------------------------------------------------------------

#endif
