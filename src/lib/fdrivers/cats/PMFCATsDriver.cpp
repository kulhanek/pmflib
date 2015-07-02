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

//subroutine pmf_cats_update_box(a,b,c,alpha,beta,gamma)
void pmf_cats_update_box_(double* a,double* b,double* c,double* alpha,double* beta,double* gamma);

//subroutine pmf_cats_update_x(natoms,x)
void pmf_cats_update_x_(int* natoms,double* x);

//subroutine pmf_cats_get_num_of_cvs(numofcvs)
void pmf_cats_get_num_of_cvs_(int* numofcvs);

//subroutine pmf_cats_get_value(value,name)
void pmf_cats_get_value_(double* value,char* name,int name_len);

//subroutine pmf_cats_get_value_by_indx(value,indx)
void pmf_cats_get_value_by_indx_(double* value,int* indx);

//subroutine pmf_cats_get_name(name,indx)
void pmf_cats_get_name_(char* name,int* indx,int name_len);

//subroutine pmf_cats_get_type_by_indx(ctype,indx)
void pmf_cats_get_type_by_indx_(char* ctype,int* indx,int ctype_len);

//subroutine pmf_cats_get_type(ctype,name)
void pmf_cats_get_type_(char* ctype,char* name,int ctype_len,int name_len);

//subroutine pmf_cats_finalize
void pmf_cats_finalize_(void);

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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFCATsDriver::SetCoordinates(int numofatoms,double* coords,double a,double b, double c, double alpha, double beta, double gamma)
{
    pmf_cats_update_box_(&a,&b,&c,&alpha,&beta,&gamma);
    pmf_cats_update_x_(&numofatoms,coords);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFCATsDriver::GetNumberOfCVs(void)
{
   int numofcvs = 0;
   pmf_cats_get_num_of_cvs_(&numofcvs);
   return(numofcvs);
}

//------------------------------------------------------------------------------

double CPMFCATsDriver::GetCVValue(CSmallString name)
{
    double value = 0.0;
    pmf_cats_get_value_(&value,name.GetBuffer(),name.GetLength());
    return(value);
}

//------------------------------------------------------------------------------

double CPMFCATsDriver::GetCVValue(int indx)
{
    indx++; // c->fortran indexing
    double value = 0.0;
    pmf_cats_get_value_by_indx_(&value,&indx);
    return(value);
}

//------------------------------------------------------------------------------

CSmallString CPMFCATsDriver::GetCVName(int indx)
{
    indx++; // c->fortran indexing
    CSmallString name;
    // FIXME
    name.SetLength(50); // PMF_MAX_CV_NAME
    pmf_cats_get_name_(name.GetBuffer(),&indx,name.GetLength());
    return(name);
}

//------------------------------------------------------------------------------

CSmallString CPMFCATsDriver::GetCVType(CSmallString name)
{
    CSmallString ctype;
    // FIXME
    ctype.SetLength(10); // PMF_MAX_TYPE
    pmf_cats_get_type_(ctype.GetBuffer(),name.GetBuffer(),ctype.GetLength(),name.GetLength());
    return(ctype);
}

//------------------------------------------------------------------------------

CSmallString CPMFCATsDriver::GetCVType(int indx)
{
    indx++; // c->fortran indexing
    CSmallString ctype;
    // FIXME
    ctype.SetLength(10); // PMF_MAX_TYPE
    pmf_cats_get_type_by_indx_(ctype.GetBuffer(),&indx,ctype.GetLength());
    return(ctype);
}

//------------------------------------------------------------------------------

void CPMFCATsDriver::Finalize(void)
{
    pmf_cats_finalize_();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


