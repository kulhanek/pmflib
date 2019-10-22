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

#include <PMFTopology.hpp>
#include <ErrorSystem.hpp>
#include <string.h>
#include <errno.h>
#include <FortranIO.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFTopology::CPMFTopology(void)
{
    BoxPresent = false;
}

//------------------------------------------------------------------------------

CPMFTopology::~CPMFTopology(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFTopology::Init(int natoms,int nres,bool has_box,const CPoint& box_centre)
{
    try {
        Atoms.CreateVector(natoms);
        Residues.CreateVector(nres);
    } catch(...) {
        ES_ERROR("unable to allocate arrays for atoms and residues");
        return(false);
    }

    BoxPresent = has_box;
    BoxCenter = box_centre;
    return(true);
}

//------------------------------------------------------------------------------

void CPMFTopology::Clear(void)
{
    Atoms.FreeVector();
    Residues.FreeVector();
    BoxPresent = false;
    BoxCenter = CPoint();
}

//------------------------------------------------------------------------------

bool CPMFTopology::SetResidue(int idx,const CSmallString& name,int first_atom)
{
    if((idx < 0) || (idx >= (int)Residues.GetLength())) {
        ES_ERROR("residue out-of-range");
        return(false);
    }
    Residues[idx].Index = idx;
    Residues[idx].Name = name;
    Residues[idx].FirstAtomIndex = first_atom;
    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::SetAtom(int idx,const CSmallString& name,const CSmallString& type)
{
    if((idx < 0) || (idx >= (int)Atoms.GetLength())) {
        ES_ERROR("atom out-of-range");
        return(false);
    }
    Atoms[idx].Index = idx;
    Atoms[idx].Name = name;
    Atoms[idx].Type = type;
    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::SetAtom(int idx,double mass,double x,double y,double z)
{
    if((idx < 0) || (idx >= (int)Atoms.GetLength())) {
        ES_ERROR("atom out-of-range");
        return(false);
    }
    Atoms[idx].Mass = mass;
    Atoms[idx].Position = CPoint(x,y,z);
    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::GetAtom(int idx,double& mass,double& x,double& y,double& z,
                CSmallString& name,CSmallString& type,
                int& resid,CSmallString& resname)
{
    if((idx < 0) || (idx >= (int)Atoms.GetLength())) {
        ES_ERROR("atom out-of-range");
        return(false);
    }
    mass = Atoms[idx].Mass;
    x = Atoms[idx].Position.x;
    y = Atoms[idx].Position.y;
    z = Atoms[idx].Position.z;
    name = Atoms[idx].Name;
    type = Atoms[idx].Type;

    CPMFResidue* p_res = Atoms[idx].GetResidue();
    if(p_res != NULL) {
        resid = p_res->GetIndex()+1;
        resname = p_res->GetName();
    } else {
        resid = 1;
        resname = "XXX";
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::Finalize(void)
{
    for(int i=0; i < (int)Residues.GetLength(); i++) {
        int last_atm;
        if(i == (int)Residues.GetLength() - 1) {
            last_atm = (int)Atoms.GetLength();
        } else {
            last_atm = Residues[i+1].FirstAtomIndex;
        }
        Residues[i].NumOfAtoms = last_atm - Residues[i].FirstAtomIndex;
        for(int j= Residues[i].FirstAtomIndex; j < last_atm; j++) {
            Atoms[j].Residue = &Residues[i];
        }
    }
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFTopology::Load(const CSmallString& name)
{
    FILE* p_top = fopen(name,"rt");
    if( p_top == NULL ){
        CSmallString error;
        error << "unable to open topology file '" << name << "' (";
        error << strerror(errno) << ")";
        ES_ERROR(error);
        return(false);
    }

    CFortranIO fortranio(p_top,true);
    char        *p_sname;

    CSmallString fPOINTERS;
    CSmallString fATOM_NAME;
    CSmallString fRESIDUE_LABEL;
    CSmallString fRESIDUE_POINTER;

    bool result = true;

    while( (p_sname = fortranio.FindNewSection()) != NULL  ) {
        if( strcmp(p_sname,"%FLAG POINTERS") == 0 ) {
            fPOINTERS = fortranio.GetFormatOfSection("%FLAG POINTERS");
            if( fPOINTERS == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG POINTERS section");
                result = false;
                break;
            }
            if( LoadBasicInfo(p_top,fPOINTERS) == false ) return(false);
            continue;
        }
        //-----------------------------------
        if( strcmp(p_sname,"%FLAG ATOM_NAME") == 0 ) {
            fATOM_NAME = fortranio.GetFormatOfSection("%FLAG ATOM_NAME");
            if( fATOM_NAME == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG ATOM_NAME section");
                result = false;
                break;
            }
            if( LoadAtomNames(p_top,fATOM_NAME) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG RESIDUE_LABEL") == 0 ) {
            fRESIDUE_LABEL = fortranio.GetFormatOfSection("%FLAG RESIDUE_LABEL");
            if( fRESIDUE_LABEL == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG RESIDUE_LABEL section");
                result = false;
                break;
            }
            if( LoadResidueNames(p_top,fRESIDUE_LABEL) == false ) return(false);
            continue;
        }

        //-----------------------------------
        if( strcmp(p_sname,"%FLAG RESIDUE_POINTER") == 0 ) {
            fRESIDUE_POINTER = fortranio.GetFormatOfSection("%FLAG RESIDUE_POINTER");
            if( fRESIDUE_POINTER == NULL ) {
                ES_ERROR("unable to decode data format of %%FLAG RESIDUE_POINTER section");
                result = false;
                break;
            }
            if( LoadResidueIPRES(p_top,fRESIDUE_POINTER) == false ) return(false);
            continue;
        }
    }

    fclose(p_top);

    return(result);
}

//------------------------------------------------------------------------------

bool CPMFTopology::LoadBasicInfo(FILE* p_top,const char* p_format)
{
    CFortranIO fortranio(p_top);
    fortranio.SetFormat(p_format);

    int NATOM;     // total number of atoms
    int NTYPES;    // total number of distinct atom types
    int NBONH;     // number of bonds containing hydrogen
    int MBONA;     // number of bonds not containing hydrogen
    int NTHETH;    // number of angles containing hydrogen
    int MTHETA;    // number of angles not containing hydrogen
    int NPHIH;     // number of dihedrals containing hydrogen
    int MPHIA;     // number of dihedrals not containing hydrogen
    int NHPARM;    // currently not used  - member item is used
    int NPARM;     // currently not used  - member item is used
    int NEXT;      // number of excluded atoms
    int NRES;      // number of residues
    int NBONA;     // MBONA + number of constraint bonds
    int NTHETA;    // MTHETA + number of constraint angles
    int NPHIA;     // MPHIA + number of constraint dihedrals
    int NUMBND;    // number of unique bond types
    int NUMANG;    // number of unique angle types
    int NPTRA;     // number of unique dihedral types
    int NATYP;     // number of atom types in parameter file, see SOLTY below
    int NPHB;      // number of distinct 10-12 hydrogen bond pair types
    int IFPERT;    // set to 1 if perturbation info is to be read in
    int NBPER;     // number of bonds to be perturbed
    int NGPER;     // number of angles to be perturbed
    int NDPER;     // number of dihedrals to be perturbed
    int MBPER;     // number of bonds with atoms completely in perturbed group
    int MGPER;     // number of angles with atoms completely in perturbed group
    int MDPER;     // number of dihedrals with atoms completely in perturbed groups
    int IFBOX;     // set to 1 if standard periodic box, 2 when truncated octahedral
    int NMXRS;     // number of atoms in the largest residue
    int IFCAP;     // set to 1 if the CAP option from edit was specified

    if( fortranio.ReadInt(NATOM) == false ) {
        ES_ERROR("unable load NATOM item");
        return(false);
    }

    if( fortranio.ReadInt(NTYPES) == false ) {
        ES_ERROR("unable load NTYPES item");
        return(false);
    }

    if( fortranio.ReadInt(NBONH) == false ) {
        ES_ERROR("unable load NBONH item");
        return(false);
    }

    if( fortranio.ReadInt(MBONA) == false ) {
        ES_ERROR("unable load MBONA item");
        return(false);
    }

    if( fortranio.ReadInt(NTHETH) == false ) {
        ES_ERROR("unable load NTHETH item");
        return(false);
    }

    if( fortranio.ReadInt(MTHETA) == false ) {
        ES_ERROR("unable load MTHETA item");
        return(false);
    }

    if( fortranio.ReadInt(NPHIH) == false ) {
        ES_ERROR("unable load NPHIH item");
        return(false);
    }

    if( fortranio.ReadInt(MPHIA) == false ) {
        ES_ERROR("unable load MPHIA item");
        return(false);
    }

    if( fortranio.ReadInt(NHPARM) == false ) {
        ES_ERROR("unable load NHPARM item");
        return(false);
    }

    if( fortranio.ReadInt(NPARM) == false ) {
        ES_ERROR("unable load NPARM item");
        return(false);
    }

    if( fortranio.ReadInt(NEXT) == false ) {
        ES_ERROR("unable load NEXT item");
        return(false);
    }

    if( fortranio.ReadInt(NRES) == false ) {
        ES_ERROR("unable load atom NRES item");
        return(false);
    }

    if( fortranio.ReadInt(NBONA) == false ) {
        ES_ERROR("unable load NBONA item");
        return(false);
    }

    if( fortranio.ReadInt(NTHETA) == false ) {
        ES_ERROR("unable load atom NTHETA item");
        return(false);
    }

    if( fortranio.ReadInt(NPHIA) == false ) {
        ES_ERROR("unable load NPHIA item");
        return(false);
    }

    if( fortranio.ReadInt(NUMBND) == false ) {
        ES_ERROR("unable load NUMBND item");
        return(false);
    }

    if( fortranio.ReadInt(NUMANG) == false ) {
        ES_ERROR("unable load NUMANG item");
        return(false);
    }

    if( fortranio.ReadInt(NPTRA) == false ) {
        ES_ERROR("unable load NPTRA item");
        return(false);
    }

    if( fortranio.ReadInt(NATYP) == false ) {
        ES_ERROR("unable load NATYP item");
        return(false);
    }

    if( fortranio.ReadInt(NPHB) == false ) {
        ES_ERROR("unable load NPHB item");
        return(false);
    }

    if( fortranio.ReadInt(IFPERT) == false ) {
        ES_ERROR("unable load IFPERT item");
        return(false);
    }

    if( fortranio.ReadInt(NBPER) == false ) {
        ES_ERROR("unable load NBPER item");
        return(false);
    }

    if( fortranio.ReadInt(NGPER) == false ) {
        ES_ERROR("unable load atom NGPER item");
        return(false);
    }

    if( fortranio.ReadInt(NDPER) == false ) {
        ES_ERROR("unable load NDPER item");
        return(false);
    }

    if( fortranio.ReadInt(MBPER) == false ) {
        ES_ERROR("unable load MBPER item");
        return(false);
    }

    if( fortranio.ReadInt(MGPER) == false ) {
        ES_ERROR("unable load MGPER item");
        return(false);
    }

    if( fortranio.ReadInt(MDPER) == false ) {
        ES_ERROR("unable load MDPER item");
        return(false);
    }

    if( fortranio.ReadInt(IFBOX) == false ) {
        ES_ERROR("unable load IFBOX item");
        return(false);
    }

    if( fortranio.ReadInt(NMXRS) == false ) {
        ES_ERROR("unable load NMXRS item");
        return(false);
    }

    if( fortranio.ReadInt(IFCAP) == false ) {
        ES_ERROR("unable load IFCAP item");
        return(false);
    }

    // check number of atoms
    if( (int)Atoms.GetLength() !=  NATOM ){
        CSmallString error;
        error << "inconsistent number of atoms: " << NATOM;
        error << " provided, but the system has " << Atoms.GetLength();
        ES_ERROR(error);
        return(false);
    }

    // reallocate residues
    try{
        Residues.CreateVector(NRES);
    } catch(...) {
        ES_ERROR("unable to allocate residue array");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::LoadAtomNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    for(int i=0; i < (int)Atoms.GetLength(); i++) {
        if( fortranio.ReadString(Atoms[i].Name) == false ) {
            ES_ERROR("unable to load ATOM_NAME item");
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::LoadResidueNames(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    for(int i=0; i < (int)Residues.GetLength(); i++) {
        if( fortranio.ReadString(Residues[i].Name) == false ) {
            ES_ERROR("unable to load LABRES item");
            return(false);
        }
        Residues[i].Index = i;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPMFTopology::LoadResidueIPRES(FILE* p_file,const char* p_format)
{
    CFortranIO fortranio(p_file);
    fortranio.SetFormat(p_format);

    // read IPRES records
    for(int i=0; i < GetNumberOfResidues(); i++) {
        int IPRES = 0;
        if( fortranio.ReadInt(IPRES) == false ) {
            ES_ERROR("unable to load IPRES item");
            return(false);
        }
        Residues[i].FirstAtomIndex = IPRES - 1;
        if( i != 0 ) {
            Residues[i-1].NumOfAtoms = Residues[i].FirstAtomIndex - Residues[i-1].FirstAtomIndex;
            for(int j=0; j < Residues[i-1].NumOfAtoms; j++) {
                int idx = Residues[i-1].FirstAtomIndex+j;
                if( (idx < 0) || (idx >= GetNumberOfAtoms()) ){
                    ES_ERROR("illegal IPRES record in the topology");
                    return(false);
                }
                Atoms[idx].Residue = &Residues[i-1];
            }
        }
    }

    // last residue
    int last_res = Residues.GetLength()-1;
    Residues[last_res].NumOfAtoms = Atoms.GetLength() - Residues[last_res].FirstAtomIndex;
    for(int j=0; j < Residues[last_res].NumOfAtoms; j++) {
        int idx = Residues[last_res].FirstAtomIndex+j;
        if( (idx < 0) || (idx >= GetNumberOfAtoms()) ){
            ES_ERROR("illegal IPRES record in the topology");
            return(false);
        }
        Atoms[idx].Residue = &Residues[last_res];
    }

    // cross-check - atoms
    int natoms = 0;
    for(int i=0; i < GetNumberOfResidues(); i++) {
        natoms += Residues[i].GetNumberOfAtoms();
    }
    if( natoms != GetNumberOfAtoms() ){
        ES_ERROR("cross-check failed: number of atoms does not match");
        return(false);
    }

    // cross-check - residues
    for(int i=0; i < GetNumberOfAtoms(); i++) {
        if( Atoms[i].Residue == NULL ){
            ES_ERROR("cross-check failed: atom not in residue");
            return(false);
        }
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFTopology::GetNumberOfAtoms(void) const
{
    return(Atoms.GetLength());
}
//------------------------------------------------------------------------------

CPMFAtom* CPMFTopology::GetAtom(int index)
{
    return(&Atoms[index]);
}
//------------------------------------------------------------------------------

int CPMFTopology::GetNumberOfResidues(void) const
{
    return(Residues.GetLength());
}

//------------------------------------------------------------------------------

CPMFResidue* CPMFTopology::GetResidue(int index)
{
    return(&Residues[index]);
}
//------------------------------------------------------------------------------

bool CPMFTopology::HasBox(void) const
{
    return(BoxPresent);
}

//------------------------------------------------------------------------------

const CPoint&   CPMFTopology::GetBoxCenter(void) const
{
    return(BoxCenter);
}

//------------------------------------------------------------------------------

bool CPMFTopology::PrintTopology(void)
{
    printf("#    ID    Name  ResID  Res     X [A]      Y [A]      Z [A]    Mass  Type\n");
    printf("# -------- ---- ------- ---- ---------- ---------- ---------- ------ ----\n");

    for(int i=0; i < (int)Atoms.GetLength(); i++) {

        CPMFResidue* p_res = Atoms[i].GetResidue();
        if(p_res != NULL) {
            printf("  %8d %-4s %7d %-4s %10.3f %10.3f %10.3f %6.2f %-4s\n",i+1,
                   (const char*)Atoms[i].GetName(),
                   p_res->GetIndex()+1,(const char*)p_res->GetName(),
                   Atoms[i].Position.x,Atoms[i].Position.y,Atoms[i].Position.z,
                   Atoms[i].Mass,
                   (const char*)Atoms[i].GetType());
        } else {
            printf("  %8d %-4s              %10.3f %10.3f %10.3f %6.2f %-4s\n",i+1,
                   (const char*)Atoms[i].GetName(),
                   Atoms[i].Position.x,Atoms[i].Position.y,Atoms[i].Position.z,
                   Atoms[i].Mass,
                   (const char*)Atoms[i].GetType());
        }
    }

    fflush(stdout);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
