#ifndef ABFAccumulatorH
#define ABFAccumulatorH
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

/** \brief ABF accumulator
*/

class PMF_PACKAGE CABFAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CABFAccumulator(void);
    ~CABFAccumulator(void);

// input and output methods ---------------------------------------------------
    /// load accumulator data from file with name
    void Load(const CSmallString& name);

    /// load accumulator data from file
    void Load(FILE* fin);

    /// save accumulator data from file with name
    void Save(const CSmallString& name);

    /// save accumulator data from file
    void Save(FILE* fout);

    /// save mask
    void SaveMask(FILE* fout);

    /// save mask
    void SaveMaskGnuPlot(FILE* fout);

// dimension specification ----------------------------------------------------
    /// set number of coordinates and optionaly groups, all previous data are destroyed
    void SetNumberOfCoords(int numofcoords);

    /// set coordinate data
    void SetCoordinate(int id,
                        const CSmallString& name,
                        const CSmallString& type,
                        double min_value,double max_value,int nbins);

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
    int GetNumberOfBinsWithABFLimit(int limit) const;

    /// return global index
    int GetGlobalIndex(const CSimpleVector<int>& position) const;

    /// convert point index to point position
    void GetPoint(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

//---------------------------------------------------------------------------------

    /// return number of samples for a given bin position
    int GetNumberOfABFSamples(const CSimpleVector<int>& position) const;

    /// return number of samples for a given bin index
    int GetNumberOfABFSamples(int ibin) const;

    /// return mask weight for a given bin index
    double GetMaskWeight(int ibin) const;

    /// set mask weight for a given bin index
    void SetMaskWeight(int ibin,double weight);

    /// return ABF force sum
    const double& GetABFForceSum(int icoord,int ibin) const;

    /// return ABF force square sum
    const double& GetABFForceSquareSum(int icoord,int ibin) const;

    /// set ABF force sum
    void SetABFForceSum(int icoord,int ibin,const double& value);

    /// set ABF force sum square
    void SetABFForceSquareSum(int icoord,int ibin,const double& value);

    /// set ABF force sum square
    void SetNumberOfABFSamples(int ibin,const int& value);

    /// return pointer to NSamples array
    int* GetNSamplesArray(void);

    /// return pointer to ABF force sum array
    double* GetABFForceSumArray(void);

    /// return pointer to ABF force square sum array
    double* GetABFForceSquareSumArray(void);

    /// get value for integration
    double GetIntegratedValue(int icoord,int ibin,bool error) const;

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Clear(void);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele) const;

    /// check cv info between two accumulators
    bool CheckCVSInfo(const CABFAccumulator* p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// print cvs info
    void PrintCVSInfo(std::ostream& vout);

// data transfer methods ------------------------------------------------------
    /// load only abf data from XML file
    void ReadABFData(CXMLElement* p_ele);

    /// add abf data from XML file
    void AddABFData(CXMLElement* p_ele);

    /// write ABF data to XML file
    void WriteABFData(CXMLElement* p_rele);

// arithemtic operations ------------------------------------------------------
    /// add another accumulator
    void AddABFAccumulator(const CABFAccumulator* p_accu);

    /// substract another accumulator
    void SubABFAccumulator(const CABFAccumulator* p_accu);

// section of private data ----------------------------------------------------
private:
    int             NCoords;        // number of coordinates

    CSimpleVector<CColVariable> Sizes;          // accumulator informations
    int                         TotNBins;      // number of total bins

    CSimpleVector<int>      NSamples;       // number of hits into bins
    CSimpleVector<double>   Mask;

    // all values are in internal units       // mask weights
    CSimpleVector<double>   ABFForce;       // accumulated ABF force
    CSimpleVector<double>   ABFForce2;      // accumulated ABF force squares

// helper methods -------------------------------------------------------------
    /// return index to ABFForce array for particular item and bin
    int map(int item,int bin) const;

    /// return index to GRPForce array for particular group, item, and bin
    int map_g(int group,int item,int bin) const;

    void Load_v3(FILE* fin);
};

//------------------------------------------------------------------------------

#endif
