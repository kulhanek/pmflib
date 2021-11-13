#ifndef FilterH
#define FilterH
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
#include <SimpleVector.hpp>
#include <boost/shared_ptr.hpp>
#include <PMFAccuData.hpp>

//------------------------------------------------------------------------------

/// base class for data filters

class PMF_PACKAGE CFilter {
public:
// constructor and destructor -------------------------------------------------
    CFilter(void);
    virtual ~CFilter(void);

// executive methods ----------------------------------------------------------
    /// set timestep
    void SetTimeStep(double timestep);

    /// filter data
    void RunFilter(CVectorDataPtr& inout);

    /// filter data
    virtual void RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out);

// section of private data ----------------------------------------------------
protected:
    double  TimeStep;       // in [fs]
    double  SamplingFreq;   // in [cm-1]
    double  InvFDTX;        // reciprocal timestep for derivatives
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CFilter>    CFilterPtr;

//------------------------------------------------------------------------------

#endif
