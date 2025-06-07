#ifndef CVSplineNaturalCubicH
#define CVSplineNaturalCubicH
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
#include <CVSpline.hpp>
#include <vector>

//------------------------------------------------------------------------------

class PMF_PACKAGE CCVSplineNaturalCubic : public CCVSpline {
public:
    CCVSplineNaturalCubic(void);
    virtual ~CCVSplineNaturalCubic(void);

// setup method ----------------------------------------------------------------
    /// clear all data
    virtual void Clear(void);

    /// allocate memory for data
    virtual void Allocate(int numofknots);

    /// register data point
    virtual bool AddPoint(int knotid,double alpha,double cv);

    /// finalize spline
    virtual bool Finalize(void);

// information methods ---------------------------------------------------------
    /// get CV value for given alpha
    virtual double GetCV(double alpha);

    /// get CV first derivatives for given alpha
    virtual double GetCVFirstDer(double alpha);

// section of private data -----------------------------------------------------
private:
    int                     NumOfKnots;
    std::vector<double>     X;
    std::vector<double>     Y;
    std::vector<double>     a, b, c, d;
};

//------------------------------------------------------------------------------

typedef std::shared_ptr<CCVSplineNaturalCubic> CCVSplineNaturalCubicPtr;

//------------------------------------------------------------------------------

#endif
