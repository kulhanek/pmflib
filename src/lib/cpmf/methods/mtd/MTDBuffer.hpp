#ifndef MTDBufferH
#define MTDBufferH
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

//------------------------------------------------------------------------------

class CXMLBinData;

//------------------------------------------------------------------------------

/*! \brief history list buffer
 */

class PMF_PACKAGE CMTDBuffer {
public:
// constructor and destructor -------------------------------------------------
    CMTDBuffer(void);
    ~CMTDBuffer(void);

// dimension specification ----------------------------------------------------
    /// allocate data
    void AllocateData(unsigned int nhills,unsigned int ncvs);

    /// destroy all previous data
    void DeallocateData(void);

// access data methods --------------------------------------------------------
    /// return number of coordinates
    unsigned int GetNumberOfCoords(void) const;

    /// return number of hills in this buffer
    unsigned int GetNumberOfHills(void) const;

    /// increment number of valid hills
    void IncNumberOfHills(void);

    /// return level ID
    unsigned int GetLevel(void) const;

    /// set level ID
    void SetLevel(unsigned int level);

    /// return start index
    unsigned int GetStart(void) const;

    /// set start index
    void SetStart(unsigned int start);

    /// return CV value
    const double& GetValue(unsigned int hill,unsigned int cvs) const;

    /// return width
    const double& GetWidth(unsigned int hill,unsigned int cvs) const;

    /// return height
    const double& GetHeight(unsigned int hill) const;

    /// set CV value
    void SetValue(unsigned int hill,unsigned int cvs,double value);

    /// set width
    void SetWidth(unsigned int hill,unsigned int cvs,double width);

    /// set height
    void SetHeight(unsigned int hill,double height);

// data transfer methods ------------------------------------------------------
    /// read data from XML element
    void ReadBufferData(CXMLBinData* p_ele);

    /// write data to XML element
    void WriteBufferData(CXMLBinData* p_ele) const;

// section of private data ----------------------------------------------------
private:
    unsigned int    level_id;               // level id - used in mtd-server/client
    unsigned int    start_id;
    unsigned int    length_of_buffer;       // length of buffer
    unsigned int    number_of_cvs;          // number of CVS
    unsigned int    num_of_values;          // actual number of snapshots
    double*         values;                 // CV values
    double*         widths;                 // gaussian widths for each CV
    double*         heights;                // gaussian heights
};

//------------------------------------------------------------------------------

#endif
