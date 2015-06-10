// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2015 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include "PMFCATsDriver.hpp"

//------------------------------------------------------------------------------
extern "C" {
// FORTRAN INTERFACE ===========================================================
//subroutine pmf_cats_begin_init(mdin,anatom,anres, &
//                            antb,box_a,box_b,box_c,box_alpha,box_beta,box_gamma)
void pmf_cats_begin_init_(char* mdin,int* anatom,int* anres,int* antb,
                         double* box_a,double* box_b,double* box_c,
                         double* box_alpha,double* box_beta,double* box_gamma,
                         int len_mdin);

//subroutine pmf_cats_set_residue(idx,name,first_atom)
void pmf_cats_set_residue_(int* idx,char* name,int* first_atom,int len_name);

//subroutine pmf_cats_set_atom(idx,name,atype)
void pmf_cats_set_atom_(int* idx,char* name,char* atype,int len_name,int len_atype);

//subroutine pmf_cats_end_init(amass,ax)
void pmf_cats_end_init_(int* anatom,double* amass,double* ax);
// FORTRAN INTERFACE ===========================================================
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFCATsDriver::BeginInit(CSmallString mdin,int anatom,int anres,
             int antb,double box_a,double box_b,double box_c,
             double box_alpha,double box_beta,double box_gamma)
{
    pmf_cats_begin_init_(mdin.GetBuffer(),&anatom,&anres,&antb,&box_a,&box_b,&box_c,
                         &box_alpha,&box_beta,&box_gamma,mdin.GetLength());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFCATsDriver::SetResidue(int idx,CSmallString name,int first_atom)
{
    idx++; // c->fortran indexing
    first_atom++;
    pmf_cats_set_residue_(&idx,name.GetBuffer(),&first_atom,name.GetLength());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFCATsDriver::SetAtom(int idx,CSmallString name,CSmallString type)
{
    idx++; // c->fortran indexing
    pmf_cats_set_atom_(&idx,name.GetBuffer(),type.GetBuffer(),name.GetLength(),type.GetLength());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFCATsDriver::EndInit(int anatom,std::vector<double>& amass,std::vector<double>& xyz)
{
    pmf_cats_end_init_(&anatom,amass.data(),xyz.data());
}

//------------------------------------------------------------------------------

void CPMFCATsDriver::Finalize(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

