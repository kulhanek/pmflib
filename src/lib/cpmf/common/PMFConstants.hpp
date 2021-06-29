#ifndef PMFConstantsH
#define PMFConstantsH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

//------------------------------------------------------------------------------

//! gas constant 8.314 472(15) J mol-1 K-1
//real(PMFDP), parameter  :: PMF_Rgas     = 0.0019872065d0     ! kcal mol-1 K-1 = 8.314 472 / 4184

const double PMF_Rgas = 0.0019872065;

//------------------------------------------------------------------------------

#endif
