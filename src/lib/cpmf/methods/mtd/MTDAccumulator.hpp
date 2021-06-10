#ifndef MTDAccumulatorH
#define MTDAccumulatorH
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

/** \brief MTD accumulator
*/

class PMF_PACKAGE CMTDAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CMTDAccumulator(void);
    ~CMTDAccumulator(void);

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
    void SetNumOfCVs(int ncvs);

    /// set coordinate data
    void SetCV(int id,
                const CSmallString& name,
                const CSmallString& type,
                double min_value,double max_value,int nbins);

    /// set coordinate data
    void SetCV(int id,
                const CSmallString& name,
                const CSmallString& type,
                double min_value,double max_value,int nbins,
                double fconv, const CSmallString& unit);

    /// allocate all arrays specified via SetNumOfCVs and SetCV
    void FinalizeAllocation(void);

    /// deallocate all array
    void Deallocate(void);

// access data methods --------------------------------------------------------
    /// return number of coordinates
    int GetNumOfCVs(void) const;

    /// return coordinate definition
    const CColVariable* GetCV(int cv) const;

    /// return number of bins within limit
    int GetNumOfBins(int limit=0) const;

    /// return global index
    int GetGlobalIndex(const CSimpleVector<int>& position) const;

    /// convert point index to point position
    void GetPoint(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetPointRValues(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

    /// set temperature
    void SetTemperature(double temp);

    /// get temperature
    double GetTemperature(void);

    /// set temperature
    void SetEnergyUnit(double fconv,const CSmallString& unit);

//---------------------------------------------------------------------------------

    /// return number of samples for a given bin position
    int GetNumOfSamples(const CSimpleVector<int>& position) const;

    /// return number of samples for a given bin index
    int GetNumOfSamples(int ibin) const;

    /// return MTDPot
    double GetMTDPot(int ibin) const;

    /// return MTDForce
    double GetMTDForce(int icv,int ibin) const;

// -----------------

    /// set MTD force sum square
    void SetNumOfSamples(int ibin,const int& value);

    /// set MTDPot
    void SetMTDPot(int ibin,const double& value);

    /// set MTD force sum square
    void SetMTDForce(int icv,int ibin,const double& value);

// -----------------

    /// return pointer to NSamples array
    int* GetNSamplesArray(void);

    /// return pointer to MTDPot array
    double* GetMTDPotArray(void);

    /// return pointer to MTDForce sum array
    double* GetMTDForceArray(void);

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Clear(void);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele) const;

    /// check cv info between two accumulators
    bool CheckCVSInfo(const CMTDAccumulator* p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// print cvs info
    void PrintCVSInfo(std::ostream& vout);

    /// print cvs info
    void PrintCVSInfo(FILE* p_fout);

// data transfer methods ------------------------------------------------------
    /// load only abf data from XML file
    void ReadMTDData(CXMLElement* p_ele);

    /// add abf data from XML file
    void AddMTDData(CXMLElement* p_ele);

    /// write MTD data to XML file
    void WriteMTDData(CXMLElement* p_rele);

// arithmetic operations ------------------------------------------------------
    /// add another accumulator
    void AddMTDAccumulator(const CMTDAccumulator* p_accu);

    /// subtract another accumulator
    void SubMTDAccumulator(const CMTDAccumulator* p_accu);

// section of private data ----------------------------------------------------
private:
    int                         NCVs;           // number of CVs
    double                      Temperature;    // temperature
    double                      EnergyFConv;    // energy unit
    CSmallString                EnergyUnit;

    CSimpleVector<CColVariable> CVs;            // collective variables
    int                         TotNBins;       // number of total bins

    CSimpleVector<int>          NSamples;       // number of hits into bins

    // all values are in internal units
    CSimpleVector<double>       MTDPot;         // MTD potential
    CSimpleVector<double>       MTDForce;       // MTD force

// helper methods -------------------------------------------------------------
    /// return index to MICF array for particular item and bin
    int map(int item,int bin) const;

    void Load_v6(char* fline,FILE* fin);
};

//------------------------------------------------------------------------------

#endif
