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

#include <stdlib.h>
#include <MTDBuffer.hpp>
#include <XMLBinData.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDBuffer::CMTDBuffer(void)
{
    level_id           = 0;       // origin id - used in mtd-server/client
    start_id           = 0;
    length_of_buffer   = 0;        // length of buffer
    number_of_cvs      = 0;
    num_of_values      = 0;        // actual number of snapshots
    values             = NULL;     // CV values
    widths             = NULL;     // gaussian widths for each CV
    heights            = NULL;    // gaussian heights
}

//------------------------------------------------------------------------------

CMTDBuffer::~CMTDBuffer(void)
{
    DeallocateData();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDBuffer::AllocateData(unsigned int nhills,unsigned int ncvs)
{
// first destroy all previous data
    if(values != NULL) DeallocateData();

// try to allocate new ones
    values = new double[nhills*ncvs];
    widths = new double[nhills*ncvs];
    heights = new double[nhills];
    length_of_buffer   = nhills;
    number_of_cvs      = ncvs;
}

//------------------------------------------------------------------------------

void CMTDBuffer::DeallocateData(void)
{
    if(values != NULL) delete[] values;
    if(widths != NULL) delete[] widths;
    if(heights != NULL)delete[] heights;

    level_id           = 0;       // origin id - used in mtd-server/client
    start_id           = 0;
    length_of_buffer   = 0;        // length of buffer
    number_of_cvs      = 0;
    num_of_values      = 0;        // actual number of snapshots
    values             = NULL;     // CV values
    widths             = NULL;     // gaussian widths for each CV
    heights            = NULL;    // gaussian heights
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDBuffer::GetNumOfCVs(void) const
{
    return(number_of_cvs);
}

//------------------------------------------------------------------------------

int CMTDBuffer::GetNumOfHills(void) const
{
    return(num_of_values);
}

//------------------------------------------------------------------------------

int CMTDBuffer::GetLevel(void) const
{
    return(level_id);
}

//------------------------------------------------------------------------------

void CMTDBuffer::SetLevel(int level)
{
    level_id = level;
}

//------------------------------------------------------------------------------

int CMTDBuffer::GetStart(void) const
{
    return(start_id);
}

//------------------------------------------------------------------------------

void CMTDBuffer::SetStart(int start)
{
    start_id = start;
}
//------------------------------------------------------------------------------

void CMTDBuffer::IncNumberOfHills(void)
{
    num_of_values++;
}

//------------------------------------------------------------------------------

const double& CMTDBuffer::GetValue(int hill,int cvs) const
{
    return(values[hill*number_of_cvs+cvs]);
}

//------------------------------------------------------------------------------

const double& CMTDBuffer::GetWidth(int hill,int cvs) const
{
    return(widths[hill*number_of_cvs+cvs]);
}

//------------------------------------------------------------------------------

const double& CMTDBuffer::GetHeight(int hill) const
{
    return(heights[hill]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDBuffer::SetValue(int hill, int cvs, double value)
{
    values[hill*number_of_cvs+cvs] = value;
}

//------------------------------------------------------------------------------

void CMTDBuffer::SetWidth(int hill, int cvs, double width)
{
    widths[hill*number_of_cvs+cvs] = width;
}

//------------------------------------------------------------------------------

void CMTDBuffer::SetHeight(int hill, double height)
{
    heights[hill] = height;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDBuffer::ReadBufferData(CXMLBinData* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    if(p_ele->GetAttribute("level",level_id) == false) {
        RUNTIME_ERROR("unable to get buffer level id");
    }

    if(p_ele->GetAttribute("start",start_id) == false) {
        RUNTIME_ERROR("unable to get buffer start id");
    }

    if(p_ele->GetAttribute("numofhills",num_of_values) == false) {
        RUNTIME_ERROR("unable to get number of values");
    }

    int history_size = num_of_values*(1 + 2*number_of_cvs)*sizeof(double);

    if((history_size == 0) || (history_size != (int)p_ele->GetLength())) {
        RUNTIME_ERROR("inconsistent history size");
    }

    double* src = (double*)p_ele->GetData();

    if(src == NULL) {
        RUNTIME_ERROR("data array is NULL");
    }

// copy data
    for(int i=0; i < num_of_values; i++) {
        SetHeight(i,*src++);
        for(int j=0; j < number_of_cvs; j++) {
            SetValue(i,j,*src++);
            SetWidth(i,j,*src++);
        }
    }
}

//------------------------------------------------------------------------------

void CMTDBuffer::WriteBufferData(CXMLBinData* p_ele) const
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    p_ele->SetAttribute("level",level_id);
    p_ele->SetAttribute("start",start_id);
    p_ele->SetAttribute("numofhills",num_of_values);

    int history_size = num_of_values*(1 + 2*number_of_cvs)*sizeof(double);
    if(history_size == 0) return;   // no data

    double* p_array = new double[history_size];
    p_ele->SetData(p_array,history_size,true); // p_history is owner of data

// copy all data
    double* dst = p_array;

// copy data
    for(int i=0; i < num_of_values; i++) {
        *dst++ = GetHeight(i);
        for(int j=0; j < number_of_cvs; j++) {
            *dst++ = GetValue(i,j);
            *dst++ = GetWidth(i,j);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
