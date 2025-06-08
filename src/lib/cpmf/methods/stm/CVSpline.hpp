#ifndef CVSplineH
#define CVSplineH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <memory>

//------------------------------------------------------------------------------

class PMF_PACKAGE CCVSpline {
public:
    CCVSpline(void);
    virtual ~CCVSpline(void);

// setup method ----------------------------------------------------------------
    /// clear all data
    virtual void Clear(void) = 0;

    /// allocate memory for data
    virtual void Allocate(int numofknots) = 0;

    /// register data point
    virtual void AddPoint(int knotid,double alpha,double cv) = 0;

    /// finalize spline
    virtual void Finalize(void) = 0;

// information methods ---------------------------------------------------------
    /// get CV value for given alpha
    virtual double GetCV(double alpha) = 0;

    /// get CV first derivatives for given alpha
    virtual double GetCVFirstDer(double alpha) = 0;
};

//------------------------------------------------------------------------------

typedef std::shared_ptr<CCVSpline>  CCVSplinePtr;

//------------------------------------------------------------------------------

#endif
