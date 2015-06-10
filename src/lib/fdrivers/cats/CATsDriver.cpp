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

#include "CATsDriver.hpp"

//------------------------------------------------------------------------------

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
void pmf_cats_end_init_(double* amass,double* ax);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCATsDriver::BeginInit(const CSmallString& mdin,int anatom,int anres,
             int antb,double box_a,double box_b,double box_c,
             double box_alpha,double box_beta,double box_gamma)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCATsDriver::SetResidue(int idx,const CSmallString& name,int first_atom)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCATsDriver::SetAtom(int idx,const CSmallString& name,const CSmallString& type)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCATsDriver::EndInit(std::vector<double>& amass,std::vector<double>& xyz)
{

}

//------------------------------------------------------------------------------

void CCATsDriver::Finalize(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


