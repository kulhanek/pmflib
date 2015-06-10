#ifndef CATsDriverH
#define CATsDriverH
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

#include <PMFMainHeader.hpp>
#include <SmallString.hpp>
#include <vector>

//------------------------------------------------------------------------------

class PMF_PACKAGE CCATsDriver {
public:
// setup methods ---------------------------------------------------------------
    static void BeginInit(const CSmallString& mdin,int anatom,int anres,
                 int antb,double box_a,double box_b,double box_c,
                 double box_alpha,double box_beta,double box_gamma);
    static void SetResidue(int idx,const CSmallString& name,int first_atom);
    static void SetAtom(int idx,const CSmallString& name,const CSmallString& type);
    static void EndInit(std::vector<double>& amass,std::vector<double>& xyz);

    static void Finalize(void);
};

//------------------------------------------------------------------------------

#endif
