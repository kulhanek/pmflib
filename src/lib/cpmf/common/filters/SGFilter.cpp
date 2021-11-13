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

#include <SGFilter.hpp>
#include <iostream>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// precalculated by MatLab
// https://www.mathworks.com/help/signal/ref/sgolay.html
// SG_2_3_5 - DOES NOT CONTAIN factorial(2) factor


double CSGFilter::SG_0_3_5[] = {
   -0.0857142857142858,
    0.3428571428571430,
    0.4857142857142860,
    0.3428571428571430,
   -0.0857142857142858
};

double CSGFilter::SG_1_3_5[] = {
    0.0833333333333332,
   -0.6666666666666670,
    0.0000000000000000,
    0.6666666666666670,
   -0.0833333333333336
};

double CSGFilter::SG_2_3_5[] = {
    0.1428571428571430,
   -0.0714285714285715,
   -0.1428571428571430,
   -0.0714285714285714,
    0.1428571428571430
};

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSGFilter::CSGFilter(void)
{
    FrameLen = 0;
    Order = 0;

    sg_0 = NULL;
    sg_1 = NULL;
    sg_2 = NULL;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CSGFilter::SetFilter(double timestep, int framelen, int order)
{
    SetTimeStep(timestep);

    if( (framelen == 5) && (order == 3) ){
        sg_0 = SG_0_3_5;
        sg_1 = SG_1_3_5;
        sg_2 = SG_2_3_5;
    } else {
        RUNTIME_ERROR("unsupported combination of SG order and framelen");
    }

    FrameLen = framelen;
}

//------------------------------------------------------------------------------

void CSGFilter::RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out)
{
    if( in->GetLength() != out->GetLength() ){
        RUNTIME_ERROR("sizes of in and out arrays differ");
    }
    int len = in->GetLength();
    if( len <= FrameLen ){
        RUNTIME_ERROR("size of in array is too small");
    }

    for(int t=0; t < FrameLen-1; t++){
        out->GetRawDataField()[t] = 0.0;
    }
    for(int t=FrameLen-1; t < len; t++){
        out->GetRawDataField()[t] = GetValue(in,t);
    }
}

//------------------------------------------------------------------------------

double CSGFilter::GetValue(const CVectorDataPtr& in,int t)
{
    double value = 0.0;
    t = t - FrameLen + 1;
    for(int idx=0; idx < FrameLen; idx++){
        // std::cout << in->GetRawDataField()[t+idx] << " " << sg_0[idx] << std::endl;
        value = value + in->GetRawDataField()[t+idx]*sg_0[idx];
    }
    return(value);
}

//------------------------------------------------------------------------------

double CSGFilter::Get1stDer(const CVectorDataPtr& in,int t)
{
    double value = 0.0;
    t = t - FrameLen + 1;
    for(int idx=0; idx < FrameLen; idx++){
        value = value + in->GetRawDataField()[t+idx]*sg_1[idx];
    }

    value = value * InvFDTX;

    return(value);
}

//------------------------------------------------------------------------------

double CSGFilter::Get2ndDer(const CVectorDataPtr& in,int t)
{
    double value = 0.0;
    t = t - FrameLen + 1;
    for(int idx=0; idx < FrameLen; idx++){
        value = value + in->GetRawDataField()[t+idx]*sg_2[idx];
    }

    // include factorial(2)
    value = 2.0 * value * InvFDTX * InvFDTX;

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
