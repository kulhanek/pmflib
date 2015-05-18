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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include <PMFMask.hpp>
#include <PMFTopology.hpp>
#include <PMFMaskSelection.hpp>
#include "maskparser/MaskParser.hpp"

#include <FortranIO.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFMask::CPMFMask(void)
{
    Topology = NULL;
    Selection = NULL;
}

//------------------------------------------------------------------------------

CPMFMask::~CPMFMask(void)
{
    if(Selection != NULL) delete Selection;
    Selection = NULL;
    Topology = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFMask::AssignTopology(CPMFTopology* p_top)
{
    if(Selection != NULL) delete Selection;
    Topology = p_top;
    Selection = NULL;
    Mask = NULL;
    return(true);
}

//------------------------------------------------------------------------------

CPMFTopology* CPMFMask::GetTopology(void) const
{
    return(Topology);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFMask::SelectAllAtoms(void)
{
    return(SetMask("@*"));
}

//------------------------------------------------------------------------------

bool CPMFMask::SetMask(const CSmallString& mask)
{
    if(Topology == NULL) {
        ES_ERROR("no topology is assigned with mask");
        return(false);
    }

    if(mask == NULL) {
        ES_ERROR("mask is empty");
        return(false);
    }

// remove previous selection -----------------------------
    Mask = mask;
    if(Selection != NULL) delete Selection;
    Selection = NULL;

// init mask parser
    init_mask();

// parse mask
    if(parse_mask(Mask) != 0) {
        free_mask_tree();
        ES_ERROR("unable to parse mask");
        return(false);
    }

// get top mask expression
    struct SExpression* p_top_expr = get_expression_tree();
    if(p_top_expr == NULL) {
        ES_ERROR("top expression is NULL");
        return(false);
    }

// init selection tree
    try {
        Selection = new CPMFMaskSelection(this);
    } catch(...) {
        free_mask_tree();
        ES_ERROR("unable to allocate root selection");
        return(false);
    }

    if(Selection->ExpandAndReduceTree(p_top_expr) == false) {
        ES_ERROR("unable to expand and reduce expression tree");
        free_mask_tree();
        return(false);
    }

// free parser data
    free_mask_tree();

    return(true);
}

//---------------------------------------------------------------------------

const CSmallString& CPMFMask::GetMask(void) const
{
    return(Mask);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFMask::GetNumberOfTopologyAtoms(void)
{
    if(Topology == NULL) return(0);
    return(Topology->GetNumberOfAtoms());
}

//------------------------------------------------------------------------------

int CPMFMask::GetNumberOfSelectedAtoms(void)
{
    if(Selection == NULL) return(0);
    return(Selection->GetNumberOfSelectedAtoms());
}

//------------------------------------------------------------------------------

CPMFAtom* CPMFMask::GetSelectedAtom(int index)
{
    if(Selection == NULL) return(NULL);
    return(Selection->GetSelectedAtom(index));
}

//------------------------------------------------------------------------------

bool CPMFMask::IsAtomSelected(int index)
{
    if(Selection == NULL) return(false);
    return(Selection->GetSelectedAtom(index) != NULL);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


