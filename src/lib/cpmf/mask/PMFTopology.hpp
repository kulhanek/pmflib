#ifndef PMFTopologyH
#define PMFTopologyH
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
#include <PMFResidue.hpp>
#include <PMFAtom.hpp>
#include <SimpleVector.hpp>

//---------------------------------------------------------------------------

/** \brief simple system topology for atom selection by masks
*/

class PMF_PACKAGE CPMFTopology {
public:
// constructor and destructor -------------------------------------------------
    CPMFTopology(void);
    ~CPMFTopology(void);

// setup topology -------------------------------------------------------------
    bool    Init(int natoms,int nres,bool has_box,const CPoint& box_centre);
    void    Clear(void);

    bool    SetResidue(int idx,const CSmallString& name,int first_atom);
    bool    SetAtom(int idx,const CSmallString& name,const CSmallString& type);
    bool    SetAtom(int idx,double mass,double x,double y,double z);
    bool    GetAtom(int idx,double& mass,double& x,double& y,double& z,
                    CSmallString& name,CSmallString& type,
                    int& resid,CSmallString& resname);

    bool    Finalize(void);

// setup topology -------------------------------------------------------------
    bool    Load(const CSmallString& name);

// informational methods ------------------------------------------------------
    int             GetNumberOfAtoms(void) const;
    CPMFAtom*       GetAtom(int index);

    int             GetNumberOfResidues(void) const;
    CPMFResidue*    GetResidue(int index);

    bool            HasBox(void) const;
    const CPoint&   GetBoxCenter(void) const;

    bool PrintTopology(void);

// section of private data ----------------------------------------------------
private:
    CSimpleVector<CPMFAtom>         Atoms;
    CSimpleVector<CPMFResidue>      Residues;
    bool                            BoxPresent;
    CPoint                          BoxCenter;

    // helper methods for loading Amber7 topology
    bool LoadBasicInfo(FILE* p_top,const char* p_format);
    bool LoadAtomNames(FILE* p_file,const char* p_format);
    bool LoadResidueNames(FILE* p_file,const char* p_format);
    bool LoadResidueIPRES(FILE* p_file,const char* p_format);
};

//---------------------------------------------------------------------------

#endif
