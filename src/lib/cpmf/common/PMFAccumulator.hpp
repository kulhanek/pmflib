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
#include <list>
#include <map>
#include <PMFAccuData.hpp>

//------------------------------------------------------------------------------

class CPMFAccumulator;
typedef boost::shared_ptr<CPMFAccumulator>    CPMFAccumulatorPtr;

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

    /// load either single PMF accumulator or set of PMF accumulators from the trajectory and store only last "number' of them
    static std::list<CPMFAccumulatorPtr> LoadFinalSnapshots(const CSmallString& name, size_t number);

    /// load either single PMF accumulator or set of PMF accumulators from the trajectory and store only last "number' of them
    static std::list<CPMFAccumulatorPtr> LoadFinalSnapshots(FILE* fin, size_t number);

    /// load accumulator data from trajectory file with name
    void LoadSnapshot(const CSmallString& name,int index);

    /// load accumulator data from trajectory file
    void LoadSnapshot(FILE* fin,int index);

    /// combine accumulator data from other accumulator to existing data
    void Combine(CPMFAccumulatorPtr right);

    /// combine accumulator data from XML to existing data
    void Combine(CXMLElement* p_ele);

    /// save accumulator data from file with name
    void Save(const CSmallString& name);

    /// save accumulator data from file
    void Save(FILE* fout);

    /// save accumulator data to XML
    void Save(CXMLElement* p_ele);

// setup methods --------------------------------------------------------------
    /// set headers
    void SetHeaders(const CSmallString& method, const CSmallString& version, const CSmallString& driver,
                    double temp, const CSmallString& temp_unit, double temp_fconv,
                    const CSmallString& ene_unit, double ene_fconv);

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

    /// convert point index to point position - real values
    void GetPointRValues(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

    /// get temperature in iu
    double GetTemperature(void) const;

    /// get temperature in unit
    double GetRealTemperature(void) const;

    /// get temperature unit
    const CSmallString& GetTemperatureUnit(void) const;

    /// get energy conversion factor to given unit
    double GetEnergyFConv(void);

    /// get number of samples
    int GetNumOfSamples(int ibin) const;

    /// get total number of samples
    int GetTotalNumOfSamples(void) const;

    /// get real value of energy
    double GetEnergyRealValue(double value) const;

    // set NCorr
    void SetNCorr(double ncorr);

    // get NCorr
    double GetNCorr(void) const;

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Reset(void);

// section data access -------------------------------------------------------
    /// check if the section data exists
    bool HasSectionData(const CSmallString& name) const;

    /// delete section data
    void DeleteSectionData(const CSmallString& name);

    /// create section
    void CreateSectionData(const CSmallString& name,const CSmallString& op,const CSmallString& type,
                           const CSmallString& mode,int len=0);

    /// create section
    void CreateSectionData(const CSmallString& name,const CSmallString& op,const CSmallString& type,
                           const CSmallString& mode,const CSmallString& mxname);

    /// create section
    void CreateSectionData(const CSmallString& name,const CSmallString& op,const CSmallString& type,
                           const CSmallString& mode,const CSmallString& mxname,const CSmallString& myname);

    /// get access to data section
    CPMFAccuDataPtr GetSectionData(const CSmallString& name) const;

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
    bool CheckCVSInfo(CPMFAccumulatorPtr p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// update number of bins
    void UpdateNumOfBins(void);

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

    double                          NCorr;

    int                             NumOfCVs;           // number of CVs
    std::vector<CColVariablePtr>    CVs;                // collective variables
    int                             NumOfBins;          // number of total bins

    std::map<CSmallString,CPMFAccuDataPtr>  DataBlocks;

// helper methods -------------------------------------------------------------
    /// return index to MICF array for particular item and bin
    int map(int item,int bin) const;

    void Load_v6(char* fline,FILE* fin);

    /// is it a header section?
    static bool IsHeaderSection(const CSmallString& keyline);

    /// is it a trajectory header?
    static bool IsTrajectoryHeader(const CSmallString& keyline);

    /// is it a trajectory snapshot?
    static bool IsTrajectorySnapshot(const CSmallString& keyline);

    /// get section name from keyline
    const CSmallString GetSectionName(const CSmallString& keyline) const;

    /// read header
    void ReadHeaderSection(FILE* fin,const CSmallString& keyline);

    /// read section data
    void ReadDataSection(FILE* fin,const CSmallString& keyline);

    /// combine two sections
    CPMFAccuDataPtr Combine(CPMFAccumulatorPtr right,CPMFAccuDataPtr ldb,CPMFAccuDataPtr rdb);
};

//------------------------------------------------------------------------------

#endif
