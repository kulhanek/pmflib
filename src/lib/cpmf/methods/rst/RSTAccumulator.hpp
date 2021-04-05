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
#include <ColVariable.hpp>
#include <SmallString.hpp>
#include <stdio.h>
#include <SimpleVector.hpp>

//------------------------------------------------------------------------------

class PMF_PACKAGE CRSTAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CRSTAccumulator(void);
    ~CRSTAccumulator(void);

// input and output methods ---------------------------------------------------
    /// load accumulator data from file with name
    void Load(const CSmallString& name);

    /// load accumulator data from file
    void Load(FILE* fin);

    /// save accumulator data from file with name
    void Save(const CSmallString& name);

    /// save accumulator data from file
    void Save(FILE* fout);

// dimension specification ----------------------------------------------------
    /// set number of coordinates, all previous data are destroyed
    void SetNumOfCVs(int numofcoords);

    /// set coordinate data
    void SetCV(int id,
                        const CSmallString& name,
                        const CSmallString& type,
                        double min_value,double max_value,int nbins);

    /// allocate all arrays specified via SetNumOfCoordinates and SetCV
    void FinalizeAllocation(void);

    /// deallocate all array
    void Deallocate(void);

// access data methods --------------------------------------------------------
    /// return number of coordinates
    int GetNumOfCVs(void) const;

    /// return total number of bins
    int GetNumOfBins(void) const;

    /// return number of bins with samples >= limit
    int GetNumOfBinsWithLimit(int limit) const;

    /// return coordinate definition
    const CColVariable* GetCV(int cv) const;

    /// return number of samples for a given bin position
    int GetNumberOfSamples(const CSimpleVector<int>& position) const;

    /// return global index
    int GetGlobalIndex(const CSimpleVector<int>& position) const;

    /// return number of samples for a given bin index
    int GetNumberOfSamples(int ibin) const;

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Clear(void);

// section of private data ----------------------------------------------------
private:
    int                         NCoords;    // number of coordinates
    CSimpleVector<CColVariable> Sizes;      // accumulator informations
    int                         TotNBins;   // number of total bins
    CSimpleVector<int>          NSamples;   // number of hits into bins
};

//------------------------------------------------------------------------------

#endif
