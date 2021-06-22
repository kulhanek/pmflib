#ifndef PMFAccumulatorH
#define PMFAccumulatorH
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
#include <ColVariable.hpp>
#include <SmallString.hpp>
#include <stdio.h>
#include <SimpleVector.hpp>
#include <vector>
#include <map>

//------------------------------------------------------------------------------

/** \brief PMF accumulator data segment
*/

class PMF_PACKAGE CPMFAccuData {
public:
    CPMFAccuData(int nbins, int ncvs);

// I/O operation ---------------------------------------------------------------
    /// load data section
    void Load(FILE* p_fin,const CSmallString& keyline);

    /// save data section
    void Save(FILE* p_fout);

// setup methods ---------------------------------------------------------------
    // reset data to zero
    void Reset(void);

// access methods --------------------------------------------------------------
    // get data section name
    const CSmallString& GetName(void) const;

    /// get data
    double GetData(int ibin, int cv=0);

    /// set data
    void SetData(int ibin, double value);

    /// set data
    void SetData(int ibin, int cv, double value);

// section of private data -----------------------------------------------------
private:
    int                     NBins;
    int                     NCVs;
    CSmallString            Name;   // name of the section
    CSmallString            Op;     // data operation
    CSmallString            Type;   // data type: R - real, I - integer
    int                     Size;   // size of data
    CSimpleVector<double>   Data;   // all data are kept as real numbers
};

typedef std::shared_ptr<CPMFAccuData>    CPMFAccuDataPtr;

//------------------------------------------------------------------------------

/** \brief PMF accumulator
*/

class PMF_PACKAGE CPMFAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CPMFAccumulator(void);
    ~CPMFAccumulator(void);

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

    /// return number of coordinates
    int GetNumOfCVs(void) const;

    /// return coordinate definition
    const CColVariablePtr GetCV(int cv) const;

    /// destroy all data
    void Clear(void);

// access data methods --------------------------------------------------------

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

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Reset(void);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele) const;

    /// check cv info between two accumulators
    bool CheckCVSInfo(const CPMFAccumulator* p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// print accu info
    void PrintAccuInfo(std::ostream& vout);

    /// print cvs info
    void PrintCVSInfo(std::ostream& vout);

// section of private data ----------------------------------------------------
private:
    CSmallString                    Version;
    CSmallString                    Method;
    CSmallString                    Driver;

    double                          Temperature;        // temperature
    double                          TemperatureFConv;
    CSmallString                    TemperatureUnit;
    double                          EnergyFConv;        // energy unit
    CSmallString                    EnergyUnit;

    int                             NCVs;           // number of CVs
    std::vector<CColVariablePtr>    CVs;            // collective variables
    int                             TotNBins;       // number of total bins

    std::map<CSmallString,CPMFAccuDataPtr>  DataBlocks;

// helper methods -------------------------------------------------------------
    /// return index to MICF array for particular item and bin
    int map(int item,int bin) const;

    void Load_v6(char* fline,FILE* fin);

    /// is it a header section?
    bool IsHeaderSection(const CSmallString& keyline);

    /// get section name from keyline
    const CSmallString GetSectionName(const CSmallString& keyline) const;

    ///
    void ReadHeaderSection(FILE* fin,const CSmallString& keyline);

    ///
    void ReadDataSection(FILE* fin,const CSmallString& keyline);
};

//------------------------------------------------------------------------------

#endif
