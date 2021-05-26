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
#include <iostream>

// keep in sync with pmf_sizes + 1 (for "C" null terminator)

#define PMF_MAX_CV_NAME 51
#define PMF_MAX_TYPE    11

//------------------------------------------------------------------------------
extern "C" {
// FORTRAN INTERFACE ===========================================================

void pmf_cats_begin_init_(char* mdin,FTINT* anatom,FTINT* anres,FTINT* antb,
                         double* box_a,double* box_b,double* box_c,
                         double* box_alpha,double* box_beta,double* box_gamma,
                         UFTINT len_mdin);

//subroutine pmf_cats_set_residue(idx,name,first_atom)
void pmf_cats_set_residue_(FTINT* idx,char* name,FTINT* first_atom,UFTINT len_name);

//subroutine pmf_cats_set_atom(idx,name,atype)
void pmf_cats_set_atom_(FTINT* idx,char* name,char* atype,UFTINT len_name,UFTINT len_atype);

//subroutine pmf_cats_end_init(amass,ax)
void pmf_cats_end_init_(FTINT* anatom,double* amass,double* ax);

//subroutine pmf_cats_update_box(a,b,c,alpha,beta,gamma)
void pmf_cats_update_box_(double* a,double* b,double* c,double* alpha,double* beta,double* gamma);

//subroutine pmf_cats_update_x(natoms,x)
void pmf_cats_update_x_(FTINT* natoms,double* x);

//subroutine pmf_cats_get_num_of_cvs(numofcvs)
void pmf_cats_get_num_of_cvs_(FTINT* numofcvs);

//subroutine pmf_cats_get_value(value,name)
void pmf_cats_get_value_(double* value,char* name,UFTINT name_len);

//subroutine pmf_cats_get_value_by_indx(value,indx)
void pmf_cats_get_value_by_indx_(double* value,FTINT* indx);

//subroutine pmf_cats_get_name(name,indx)
void pmf_cats_get_name_(char* name,FTINT* indx,UFTINT name_len);

//subroutine pmf_cats_get_type_by_indx(ctype,indx)
void pmf_cats_get_type_by_indx_(char* ctype,FTINT* indx,UFTINT ctype_len);

//subroutine pmf_cats_get_type(ctype,name)
void pmf_cats_get_type_(char* ctype,char* name,UFTINT ctype_len,UFTINT name_len);

//subroutine pmf_cats_finalize
void pmf_cats_finalize_(FTINT* mode);

// FORTRAN INTERFACE ===========================================================
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE CPMFCATsDriver::BeginInit(CSmallString mdin,int anatom,int anres,
             int antb,double box_a,double box_b,double box_c,
             double box_alpha,double box_beta,double box_gamma)
{
    // flush C and C++ streams;
    fflush(stdout);
    fflush(stderr);
    std::cout.flush();
    std::cerr.flush();

    FTINT lnatom = anatom;
    FTINT lnres = anres;
    FTINT lntb = antb;

    pmf_cats_begin_init_(mdin.GetBuffer(),&lnatom,&lnres,&lntb,&box_a,&box_b,&box_c,
                         &box_alpha,&box_beta,&box_gamma,mdin.GetLength());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE CPMFCATsDriver::SetResidue(int idx,CSmallString name,int first_atom)
{
    FTINT lidx = idx;
    lidx++; // c->fortran indexing
    FTINT lfirst_atom = first_atom;
    lfirst_atom++;
    pmf_cats_set_residue_(&lidx,name.GetBuffer(),&lfirst_atom,name.GetLength());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE CPMFCATsDriver::SetAtom(int idx,CSmallString name,CSmallString type)
{
    FTINT lidx = idx;
    lidx++; // c->fortran indexing
    pmf_cats_set_atom_(&lidx,name.GetBuffer(),type.GetBuffer(),name.GetLength(),type.GetLength());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE CPMFCATsDriver::EndInit(int anatom,std::vector<double>& amass,std::vector<double>& xyz)
{
    FTINT lnatom = anatom;

    pmf_cats_end_init_(&lnatom,amass.data(),xyz.data());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void PMF_PACKAGE CPMFCATsDriver::SetCoordinates(int numofatoms,double* coords,double a,double b, double c, double alpha, double beta, double gamma)
{
    FTINT lnumofatoms = numofatoms;

    pmf_cats_update_box_(&a,&b,&c,&alpha,&beta,&gamma);
    pmf_cats_update_x_(&lnumofatoms,coords);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int PMF_PACKAGE CPMFCATsDriver::GetNumberOfCVs(void)
{
   FTINT numofcvs = 0;
   pmf_cats_get_num_of_cvs_(&numofcvs);
   return(numofcvs);
}

//------------------------------------------------------------------------------

double PMF_PACKAGE CPMFCATsDriver::GetCVValue(CSmallString name)
{
    double value = 0.0;
    pmf_cats_get_value_(&value,name.GetBuffer(),name.GetLength());
    return(value);
}

//------------------------------------------------------------------------------

double PMF_PACKAGE CPMFCATsDriver::GetCVValue(int indx)
{
    FTINT lidx = indx;
    lidx++; // c->fortran indexing
    double value = 0.0;
    pmf_cats_get_value_by_indx_(&value,&lidx);
    return(value);
}

//------------------------------------------------------------------------------

CSmallString PMF_PACKAGE CPMFCATsDriver::GetCVName(int indx)
{
    FTINT lidx = indx;
    lidx++; // c->fortran indexing

    CSmallString name;
    name.SetLength(PMF_MAX_CV_NAME);
    pmf_cats_get_name_(name.GetBuffer(),&lidx,name.GetLength());

    name[PMF_MAX_CV_NAME-1] = '\0';
    name.Trim();
    return(name);
}

//------------------------------------------------------------------------------

CSmallString PMF_PACKAGE CPMFCATsDriver::GetCVType(CSmallString name)
{
    CSmallString ctype;
    ctype.SetLength(PMF_MAX_TYPE);
    pmf_cats_get_type_(ctype.GetBuffer(),name.GetBuffer(),ctype.GetLength(),name.GetLength());

    ctype[PMF_MAX_TYPE-1] = '\0';
    ctype.Trim();
    return(ctype);
}

//------------------------------------------------------------------------------

CSmallString PMF_PACKAGE CPMFCATsDriver::GetCVType(int indx)
{
    FTINT lidx = indx;
    lidx++; // c->fortran indexing

    CSmallString ctype;
    ctype.SetLength(PMF_MAX_TYPE);
    pmf_cats_get_type_by_indx_(ctype.GetBuffer(),&lidx,ctype.GetLength());

    ctype[PMF_MAX_TYPE-1] = '\0';
    ctype.Trim();
    return(ctype);
}

//------------------------------------------------------------------------------

void PMF_PACKAGE CPMFCATsDriver::Finalize(int mode)
{
    FTINT lmode = mode;
    pmf_cats_finalize_(&lmode);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


