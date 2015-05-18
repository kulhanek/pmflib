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
#include <SmallString.hpp>
#include <ErrorSystem.hpp>
#include <PMFTopology.hpp>
#include <PMFMask.hpp>

//------------------------------------------------------------------------------

CPMFTopology    PMFTopology;
CPMFMask        PMFMask;

//------------------------------------------------------------------------------

extern "C" {

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_topo_init_(int* rtstat,int* natoms,int* nres,int* has_box,
                          double* cbox_x,double* cbox_y,double* cbox_z)
{
    *rtstat = 0;
    if(PMFTopology.Init(*natoms,*nres,*has_box==1,CPoint(*cbox_x,*cbox_y,*cbox_z)) == false) {
        *rtstat = -1;
        return;
    }
    if(PMFMask.AssignTopology(&PMFTopology) == false) {
        *rtstat = -1;
        return;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_topo_load_(int* rtstat,char* name,unsigned int name_len)
{
    *rtstat = 0;

    CSmallString s_name;
    s_name.SetFromFortran(name,name_len);

    if( PMFTopology.Load(s_name) ){
        return;
    } else {
        *rtstat = 1;
        return;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_clear_(void)
{
    PMFMask.AssignTopology(NULL);
    PMFTopology.Clear();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_set_topo_residue_(int* rtstat,int* index,char* name,int* first_atom,
                                 unsigned int name_len)
{
    CSmallString s_name;
    s_name.SetFromFortran(name,name_len,true);

    *rtstat = 0;
    if(PMFTopology.SetResidue(*index-1,s_name,*first_atom-1) == false) {
        *rtstat = -1;
        return;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_set_topo_atom_(int* rtstat,int* index,char* name,char* type,
                              unsigned int name_len,unsigned int type_len)
{
    CSmallString s_name;
    s_name.SetFromFortran(name,name_len,true);

    CSmallString s_type;
    s_type.SetFromFortran(type,type_len,true);

    *rtstat = 0;
    if(PMFTopology.SetAtom(*index-1,s_name,s_type) == false) {
        *rtstat = -1;
        return;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_set_topo_atom_mcrd_(int* rtstat,int* index,
                                   double* mass, double* x,double* y,double* z)
{
    *rtstat = 0;
    if(PMFTopology.SetAtom(*index-1,*mass,*x,*y,*z) == false) {
        *rtstat = -1;
        return;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_topo_finalize_(int* rtstat)
{
    *rtstat = 0;
    PMFTopology.Finalize();
    PMFMask.AssignTopology(&PMFTopology);
}

//------------------------------------------------------------------------------

void cpmf_mask_get_topo_atom_(int* ret_st,int* idx,char* name,int* resid,
                              char* resname,double* x, double* y, double* z,
                              double* mass,char* type,
                              unsigned int name_len,
                              unsigned int resname_len,
                              unsigned int type_len)
{
    *ret_st = 0;
    if( (*idx < 1) || (*idx > PMFTopology.GetNumberOfAtoms() ) ){
        *ret_st = 1;
        return;
    }

    CSmallString sname,stype,sresname;
    PMFTopology.GetAtom(*idx-1,*mass,*x,*y,*z,sname,stype,*resid,sresname);

    memset(name,' ',name_len);
    memset(resname,' ',resname_len);
    memset(type,' ',type_len);
    strncpy(name,sname, name_len < strlen(sname) ? name_len : strlen(sname) );
    strncpy(resname,sresname, resname_len < strlen(sresname) ? resname_len : strlen(sresname) );
    strncpy(type,stype, type_len < strlen(stype) ? type_len : strlen(stype));
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_natoms_in_mask_(int* rtstat,char* mask,int* natoms,unsigned int mask_len)
{
    CSmallString s_mask;
    s_mask.SetFromFortran(mask,mask_len);

    *rtstat = 0;
    *natoms = 0;
    if(PMFMask.SetMask(s_mask) == false) {
        *rtstat = -1;
        return;
    }

    *natoms = PMFMask.GetNumberOfSelectedAtoms();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_set_mask_(int* rtstat,char* mask,unsigned int mask_len)
{
    CSmallString s_mask;
    s_mask.SetFromFortran(mask,mask_len);

    *rtstat = 0;
    if(PMFMask.SetMask(s_mask) == false) {
        *rtstat = -1;
        return;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_is_atom_selected_(int* rtstat,int* idx,int* sel)
{
    *rtstat = 0;
    *sel = 0;

    int index = *idx - 1;

    if((index < 0) || (index >= PMFTopology.GetNumberOfAtoms())) {
        *rtstat = 1;
        return;
    }

    if(PMFMask.IsAtomSelected(index) == true) {
        *sel = 1;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void cpmf_mask_print_errors_(void)
{
    ErrorSystem.PrintErrors(stderr);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

}
