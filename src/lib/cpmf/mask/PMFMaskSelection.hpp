#ifndef PMFMaskSelectionH
#define PMFMaskSelectionH
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

#include <PMFMainHeader.hpp>
#include <SimpleVector.hpp>
#include "maskparser/MaskParser.hpp"

//---------------------------------------------------------------------------

class CPMFMask;
class CPMFAtom;

//---------------------------------------------------------------------------

/// selection of atoms

class PMF_PACKAGE CPMFMaskSelection {
public:
    CPMFMaskSelection(CPMFMask* p_owner);
    ~CPMFMaskSelection(void);

    // executive methods  ------------------------------------------------------
    bool ExpandAndReduceTree(struct SExpression* p_expr);

    // results -----------------------------------------------------------------
    int             GetNumberOfSelectedAtoms(void);
    CPMFAtom*       GetSelectedAtom(int index);

// section of private data ----------------------------------------------------
private:
    CPMFMask*                   Owner;
    CSimpleVector<CPMFAtom*>    Atoms;

    static bool ExpandAndReduceTree(CPMFMaskSelection* p_root,
                                    struct SExpression* p_expr);

    // individual selections
    bool Select(struct SSelection* p_sel);

    void SelectAtomByIndex(int index,int length);
    void SelectAtomByName(const char* p_name);
    void SelectAtomByType(const char* p_name);

    void SelectResidueByIndex(int index,int length);
    void SelectResidueByName(const char* p_name);

    bool SelectAtomByDistanceFromOrigin(SOperator dist_oper,double dist);
    bool SelectAtomByDistanceFromCentreOfBox(SOperator dist_oper,double dist);
    bool SelectAtomByDistanceFromList(CPMFMaskSelection* p_left,
                                      SOperator dist_oper,double dist);
    bool SelectAtomByDistanceFromCOM(CPMFMaskSelection* p_left,
                                     SOperator dist_oper,double dist);
    bool SelectAtomByDistanceFromPlane(CPMFMaskSelection* p_left,
                                       SOperator dist_oper,double dist);

    bool SelectResidueByDistanceFromOrigin(SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromCentreOfBox(SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromList(CPMFMaskSelection* p_left,
                                         SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromCOM(CPMFMaskSelection* p_left,
                                        SOperator dist_oper,double dist);
    bool SelectResidueByDistanceFromPlane(CPMFMaskSelection* p_left,
                                          SOperator dist_oper,double dist);

    int   strnlen(const char* p_s1,int len);
    bool  firstmatch(const char* p_s1,const char* p_s2,int len);
};

//---------------------------------------------------------------------------

#endif
