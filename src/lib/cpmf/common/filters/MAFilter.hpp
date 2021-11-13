#ifndef MAFilterH
#define MAFilterH
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

/// moving average filter

class PMF_PACKAGE CMAFilter : public CFilter {
public:
// constructor and destructor -------------------------------------------------
    CMAFilter(void);

// executive methods ----------------------------------------------------------
    /// setup filter
    void SetFilter(double timestep, double fcut);

    /// get frame length
    int GetFrameLen(void) const;

    /// filter data
    virtual void RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out);

// section of private data ----------------------------------------------------
protected:
    double  CutoffFreq;     // cut-off frequency in [cm-1]
    double  NormFreq;       // normalized frequency
    size_t  N;              // length of moving window
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CMAFilter>    CMAFilterPtr;

//------------------------------------------------------------------------------

#endif
