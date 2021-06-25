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

    // get data length
    int GetLength(void) const;

    // get num of bins
    int GetNumOfBins(void) const;

    // get num of cvs
    int GetNumOfCVs(void) const;

    // get data section type
    const CSmallString& GetType(void) const;

    // get data section mode
    const CSmallString& GetMode(void) const;

    // get data section op
    const CSmallString& GetOp(void) const;

    /// get data
    double GetData(int ibin, int cv=0) const;

    /// set data
    void SetData(int ibin, double value);

    /// set data
    void SetData(int ibin, int cv, double value);

// section of private data -----------------------------------------------------
private:
    int                     NumOfBins;
    int                     NumOfCVs;
    CSmallString            Name;       // name of the section
    CSmallString            Op;         // data operation
    CSmallString            Type;       // data type: R - real, I - integer
    CSmallString            Mode;       // data mode: B - per bins, C - per CVs, M - mixed per bins and cvs
    int                     Size;       // size of data
    CSimpleVector<double>   Data;       // all data are kept as real numbers
};

typedef std::shared_ptr<CPMFAccuData>    CPMFAccuDataPtr;

//------------------------------------------------------------------------------

/** \brief PMF accumulator
*/

class PMF_PACKAGE CPMFAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CPMFAccumulator(void);
    virtual ~CPMFAccumulator(void);

// input and output methods ---------------------------------------------------
    /// load accumulator data from file with name
    void Load(const CSmallString& name);

    /// load accumulator data from file
    void Load(FILE* fin);

    /// load accumulator data from XML
    void Load(CXMLElement* p_ele);

    /// load accumulator data from trajectory file with name
    void LoadSnapshot(const CSmallString& name,int index);

    /// load accumulator data from trajectory file
    void LoadSnapshot(FILE* fin,int index);

    /// combine accumulator data from XML to existing data
    void Combine(CXMLElement* p_ele);

    /// save accumulator data from file with name
    void Save(const CSmallString& name);

    /// save accumulator data from file
    void Save(FILE* fout);

    /// save accumulator data to XML
    void Save(CXMLElement* p_ele);

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

    /// get method
    const CSmallString& GetMethod(void) const;

    /// get driver
    const CSmallString& GetDriver(void) const;

    /// return number of bins
    int GetNumOfBins(void) const;

    /// return global index
    int GetGlobalIndex(const CSimpleVector<int>& position) const;

    /// convert point index to point position
    void GetPoint(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetPointRValues(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

    /// get temperature in iu
    double GetTemperature(void);

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Reset(void);

// section data access -------------------------------------------------------
    /// get data
    double GetData(const CSmallString& name, int ibin, int cv=0) const;

    /// set data
    void SetData(const CSmallString& name, int ibin, double value);

    /// set data
    void SetData(const CSmallString& name, int ibin, int cv, double value);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele) const;

    /// check cv info between two accumulators
    bool CheckCVSInfo(const CPMFAccumulator* p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// print all info
    void PrintInfo(std::ostream& vout);

    /// print all info
    void PrintInfo(FILE* p_fout);

    /// print accu info
    void PrintAccuInfo(std::ostream& vout);

    /// print accu info
    void PrintAccuInfo(FILE* p_fout);

    /// print cvs info
    void PrintCVSInfo(std::ostream& vout);

    /// print cvs info
    void PrintCVSInfo(FILE* p_fout);

    /// print list of sections
    void ListSections(std::ostream& vout);

// section of private data ----------------------------------------------------
protected:
    CSmallString                    Version;
    CSmallString                    Method;
    CSmallString                    Driver;

    double                          Temperature;        // temperature
    double                          TemperatureFConv;
    CSmallString                    TemperatureUnit;
    double                          EnergyFConv;        // energy unit
    CSmallString                    EnergyUnit;

    int                             NumOfCVs;           // number of CVs
    std::vector<CColVariablePtr>    CVs;                // collective variables
    int                             NumOfBins;          // number of total bins

    std::map<CSmallString,CPMFAccuDataPtr>  DataBlocks;

// helper methods -------------------------------------------------------------
    /// return index to MICF array for particular item and bin
    int map(int item,int bin) const;

    void Load_v6(char* fline,FILE* fin);

    /// is it a header section?
    bool IsHeaderSection(const CSmallString& keyline);

    /// get section name from keyline
    const CSmallString GetSectionName(const CSmallString& keyline) const;

    /// read header
    void ReadHeaderSection(FILE* fin,const CSmallString& keyline);

    /// read section data
    void ReadDataSection(FILE* fin,const CSmallString& keyline);
};

typedef std::shared_ptr<CPMFAccumulator>    CPMFAccumulatorPtr;

//------------------------------------------------------------------------------

#endif
