#ifndef RSTAccumulatorH
#define RSTAccumulatorH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <PMFMainHeader.hpp>
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

class PMF_PACKAGE CRSTAccumulator : public CPMFAccumulator {
public:
// constructor and destructor --------------------------------------------------
    CRSTAccumulator(void);
    ~CRSTAccumulator(void);

// -----------------------------------------------------------------------------
    /// return number of samples for a given bin index
    int GetNumOfSamples(int ibin) const;
};

//------------------------------------------------------------------------------

#endif
