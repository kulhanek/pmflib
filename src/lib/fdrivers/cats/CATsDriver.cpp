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

#include <PMFLibCVs.hpp>
#include <PrmFile.hpp>
#include <AmberTopology.hpp>

#include <cv_dis.hpp>
#include <cv_odis.hpp>
#include <cv_plane.hpp>
#include <cv_puck6q.hpp>
#include <cv_puck6t.hpp>
#include <cv_puck6p.hpp>
#include <cv_dih.hpp>
#include <cv_mdisc.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFLibCVs::CPMFLibCVs(void)
{

}

//------------------------------------------------------------------------------

bool CPMFLibCVs::Read(const CSmallString& name)
{
    CPrmFile prm_file;

    fprintf(stdout,":::::::::::::::::::::::::::::::::::::::::: {CVS} :::::::::::::::::::::::::::::::\n");
    fprintf(stdout,"PMFLib CVs file name           : %s\n",(const char*) name);

    //-----------------------------
    if( prm_file.Read(name) == false ) {
        fprintf(stderr,">>> ERROR: unable to read PMFLib file (%s)!\n",(const char*)name);
        return(false);
    }

    if( prm_file.OpenGroup("CVS") == false ) {
        fprintf(stderr,">>> ERROR: unable to open {CVS} group!\n");
        return(false);
    }

    fprintf(stdout,"Total number of atoms          : %6d\n",fnatom);
    fprintf(stdout,"Number of collective variables : %6d\n",prm_file.CountGroup());

    //-----------------------------
    bool result = prm_file.FirstSection();

    int id = 1;
    while( result ) {
        CSmallString cv_type;
        cv_type = prm_file.GetSectionName();

        fprintf(stdout,"\n== Reading collective variable #%02d of type \"%s\"\n",id,(const char*)cv_type);

        CPMFLibCV* p_cv = CreateCV(cv_type);
        if( p_cv == NULL ) {
            return(false);
        }

        if( p_cv->Read(prm_file) == false ) {
            delete p_cv;
            return(false);
        }

        InsertToEnd(p_cv,0,true);

        result = prm_file.NextSection();
        id++;
    }
    fprintf(stdout,"\n");

    //-----------------------------

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFLibCVs::GetValue(const CSmallString& name,CAmberRestart* p_crd,double& value)
{
    CSimpleIterator<CPMFLibCV> I(this);
    CPMFLibCV* p_cv;

    value = 0.0;
    while( (p_cv = I.Current()) != NULL ) {
        if( p_cv->GetName() == name ) {
            value = p_cv->GetValue(p_crd);
            return(true);
        }
        I++;
    }

    CSmallString error;
    error << "CV with name '" << name << "' does not exist";
    ES_ERROR(error);

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFLibCVs::SetTopology(CAmberTopology* p_top)
{
    fnatom = p_top->AtomList.GetNumberOfAtoms();
    mask.AssignTopology(p_top);
    return(true);
}

//------------------------------------------------------------------------------

CPMFLibCV* CPMFLibCVs::CreateCV(const CSmallString& cv_type)
{
    if( cv_type == "DIS" ) {
        return(new CV_DIS(this));

    } else  if( cv_type == "ODIS" ) {
        return(new CV_ODIS(this));

    } else  if( cv_type == "PLANE" ) {
    return(new CV_PLANE(this));

    } else  if( cv_type == "PUCK6Q" ) {
        return(new CV_PUCK6Q(this));

    } else  if( cv_type == "PUCK6T" ) {
        return(new CV_PUCK6T(this));

    } else  if( cv_type == "PUCK6P" ) {
        return(new CV_PUCK6P(this));

    } else  if( cv_type == "DIH" ) {
        return(new CV_DIH(this));

    } else  if( cv_type == "MDISC" ) {
        return(new CV_MDISC(this));

    } else {
        fprintf(stderr,">>> ERROR: CV type (%s) is not implemented!\n",(const char*)cv_type);
        return(NULL);
    }
}

//------------------------------------------------------------------------------

int CPMFLibCVs::pmf_mask_natoms_in_mask(const CSmallString& smask)
{
    if( mask.SetMask(smask) == false ) return(0);
    return(mask.GetNumberOfSelectedAtoms());
}

//------------------------------------------------------------------------------

bool CPMFLibCVs::pmf_mask_set_mask(const CSmallString& smask)
{
    return( mask.SetMask(smask) );
}

//------------------------------------------------------------------------------

bool CPMFLibCVs::pmf_mask_is_atom_selected(int i)
{
    return( mask.IsAtomSelected(i) );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


