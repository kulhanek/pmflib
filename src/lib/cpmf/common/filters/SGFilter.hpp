#ifndef SGFilterH
#define SGFilterH
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

#include <Filter.hpp>

//------------------------------------------------------------------------------

/// Savitzky-Golay filter
// https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter

class PMF_PACKAGE CSGFilter : public CFilter {
public:
// constructor and destructor -------------------------------------------------
    CSGFilter(void);

// executive methods ----------------------------------------------------------
    /// setup filter
    void SetFilter(double timestep, int framelen, int order);

    /// filter data
    virtual void RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out);

    /// get smoothed value
    double GetValue(const CVectorDataPtr& in,int t);

    /// get first derivative
    double Get1stDer(const CVectorDataPtr& in,int t);

    /// get second derivative
    double Get2ndDer(const CVectorDataPtr& in,int t);

// section of private data ----------------------------------------------------
protected:
    int     FrameLen;   // length of moving window
    int     Order;
    double* sg_0;
    double* sg_1;
    double* sg_2;

    static double SG_0_3_5[];
    static double SG_1_3_5[];
    static double SG_2_3_5[];
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CSGFilter>    CSGFilterPtr;

//------------------------------------------------------------------------------

#endif
