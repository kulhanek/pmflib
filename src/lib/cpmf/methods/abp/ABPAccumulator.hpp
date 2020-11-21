#ifndef ABPAccumulatorH
#define ABPAccumulatorH
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

#include <PMFMainHeader.hpp>
#include <ColVariable.hpp>
#include <SmallString.hpp>
#include <stdio.h>
#include <SimpleVector.hpp>

//------------------------------------------------------------------------------

/** \brief ABP accumulator
*/

class PMF_PACKAGE CABPAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CABPAccumulator(void);
    ~CABPAccumulator(void);

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
    /// set number of coordinates and optionaly groups, all previous data are destroyed
    void SetNumberOfCoords(int numofcoords);

    /// set coordinate data
    void SetCoordinate(int id,
                        const CSmallString& name,
                        const CSmallString& type,
                        double min_value,double max_value,int nbins,
                        double alpha,
                        double fconv, const CSmallString& unit);

    /// allocate all arrays specified via SetNumOfCoordinates and SetCoordinate
    void FinalizeAllocation(void);

    /// deallocate all array
    void Deallocate(void);

// access data methods --------------------------------------------------------
    /// return number of coordinates
    int GetNumberOfCoords(void) const;

    /// return coordinate definition
    const CColVariable* GetCoordinate(int cv) const;

    /// return total number of bins
    int GetNumberOfBins(void) const;

    /// return number of bins with samples >= limit
    int GetNumberOfBinsWithABPLimit(int limit) const;

    /// convert point index to point position
    void GetPoint(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

//---------------------------------------------------------------------------------

    /// return number of samples for a given bin index
    int GetNumberOfABPSamples(int ibin) const;

    /// get energy
    double GetEnergy(double temp,int ibin) const;

    /// return pointer to NSamples array
    int* GetNSamplesArray(void);

    /// return pointer to DPop array
    double* GetDPopArray(void);

    /// return pointer to Pop array
    double* GetPopArray(void);

    /// get value of pop
    double GetPop(int ibin) const;

    /// get temperature
    double GetTemperature(void) const;

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Clear(void);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele) const;

    /// check cv info between two accumulators
    bool CheckCVSInfo(const CABPAccumulator* p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// print cvs info
    void PrintCVSInfo(std::ostream& vout);

// data transfer methods ------------------------------------------------------
    /// load only abf data from XML file
    void ReadABPData(CXMLElement* p_ele);

    /// add abf data from XML file
    void AddABPData(CXMLElement* p_ele);

    /// write ABP data to XML file
    void WriteABPData(CXMLElement* p_rele);

// section of private data ----------------------------------------------------
private:
    int                         NCoords;    // number of coordinates
    CSmallString                Kernel;
    double                      Temperature;

    CSimpleVector<CColVariable> Sizes;      // accumulator information
    CSimpleVector<double>       Alphas;     // CV alphas
    int                         TotNBins;   // number of total bins

    CSimpleVector<int>          NSamples;   // number of hits into bins

    // all values are in internal units
    CSimpleVector<double>       DPop;       // accumulated derivatives
    CSimpleVector<double>       Pop;        // accumulated probability

// helper methods -------------------------------------------------------------
    /// return index to ABPForce array for particular item and bin
    int map(int item,int bin) const;

    /// return index to GRPForce array for particular group, item, and bin
    int map_g(int group,int item,int bin) const;
};

//------------------------------------------------------------------------------

#endif
