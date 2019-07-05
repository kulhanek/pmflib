#ifndef MTDHistoryH
#define MTDHistoryH
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
#include <SimpleList.hpp>
#include <MTDBuffer.hpp>
#include <SmallString.hpp>
#include <ColVariable.hpp>
#include <SimpleVector.hpp>

//// calculate energy -----------------------------------------------------------
//    /// calculate FES from MTD history list
//    void CalculateFES(CMTDHistory& mtd_hist,unsigned int mtd_time=0);

//    /// calculate FES from MTD parameter list
//    void CalculateFES(unsigned int ncoords,CSimpleVector<double>& params);

//------------------------------------------------------------------------------

/*! \brief metadynamics gaussian history list
 */

class PMF_PACKAGE CMTDHistory {
public:
// constructor and destructor -------------------------------------------------
    CMTDHistory(void);
    ~CMTDHistory(void);

// input and output methods ---------------------------------------------------
    /// load history data from file with name
    void Load(const CSmallString& name);

    /// load history data from file
    void Load(FILE* fin);

    /// save history data from file with name
    void Save(const CSmallString& name);

    /// save history data from file
    void Save(FILE* fout);

// dimension specification ----------------------------------------------------
    /// set number of coordinates, all previous data are destroyed
    void SetNumberOfCoords(int numofcoords);

    /// set energy unit
    void SetEnergyUnit(const CSmallString& unit, double unit_fac);

    /// set coordinate data
    void SetCoordinate(unsigned int id,
                       const CSmallString& name,
                       const CSmallString& type,
                       double min_value,double max_value,unsigned int nbins);

    /// deallocate all array
    void Deallocate(void);

    /// reallocate as single buffer
    void ReallocateAsSingleBuffer(void);

// access data methods --------------------------------------------------------
    /// return number of coordinates
    unsigned int GetNumberOfCoords(void) const;

    /// get energy unit
    const CSmallString& GetEnergyUnit(void) const;

    /// get energy unit conversion factor
    double GetEnergyUnitFac(void) const;

    /// return number of hills
    unsigned int GetNumberOfHills(void) const;

    /// return number of hills in selected buffer list
    unsigned int GetNumberOfHills(const CSimpleList<CMTDBuffer>& list) const;

    /// return coordinate definition
    const CColVariable* GetCoordinate(unsigned int cv) const;

    /// calculate value for given point
    double CalculateValue(const CSimpleVector<double>& point,unsigned int mtdtime=0);

    /// get height
    const double& GetHeight(unsigned int hill_index);

    /// get CV value
    const double& GetValue(unsigned int hill_index,unsigned int cv);

    /// get width
    const double& GetWidth(unsigned int hill_index,unsigned int cv);

// buffer methods -------------------------------------------------------------
    /// allocate new buffer
    CMTDBuffer* GetNewBuffer(unsigned int size);

    /// get buffer
    CMTDBuffer* GetBuffer(unsigned int index);

// parameter vector: height,[value,width] -------------------------------------
    /// get history as parameter vector
    void GetPVector(CSimpleVector<double>& pvector);

    /// set history from parameter vector
    void SetPVector(const CSimpleVector<double>& pvector);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele);

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_iele);

// data transfer methods ------------------------------------------------------
    /// load all data from XML file
    void ReadMTDData(CXMLElement* p_ele);

    /// add data from XML file to the end of history
    void AddMTDData(CXMLElement* p_ele);

    /// add data from XML file to the and of history as single buffer
    CMTDBuffer* AddMTDDataAsSingleBuffer(CXMLElement* p_ele);

    /// write all data to XML file
    void WriteMTDData(CXMLElement* p_ele) const;

    /// write selected data to XML file
    void WriteMTDData(CXMLElement* p_ele,const CSimpleList<CMTDBuffer>& list) const;

// section of private data ----------------------------------------------------
protected:
    unsigned int                NCoords;        // number of coordinates
    CSmallString                EnergyUnit;     // energy unit
    double                      EnergyUnitFac;  // energy unit factor
    unsigned int                MaxBufferSize;  // max buffer size
    CSimpleVector<CColVariable> Sizes;          // accumulator informations
    CSimpleList<CMTDBuffer>     Buffers;        // all buffers in history
};

//------------------------------------------------------------------------------

#endif
