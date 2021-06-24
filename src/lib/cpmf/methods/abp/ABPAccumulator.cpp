// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//                       Martin Petrek, petrek@chemi.muni.cz
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

#include <errno.h>
#include <string.h>
#include <ABPAccumulator.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABPAccumulator::CABPAccumulator(void)
{
}

//------------------------------------------------------------------------------

CABPAccumulator::~CABPAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPAccumulator::GetNumberOfSamples(int ibin) const
{
    return( GetData("NSAMPLES",ibin) );
}

//------------------------------------------------------------------------------

double CABPAccumulator::GetEnergy(int ibin) const
{
    // FIXME
    return(0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



