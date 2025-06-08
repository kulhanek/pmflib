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
#include <SimpleVector.hpp>

// Smoothing with cubic splines
// D. S. G. Pollock
// Queen Mary and Westfield College, The University of London

// Cubic Spline Interpolation, algorithm (18)

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
    virtual void AddPoint(int knotid,double alpha,double cv);

    /// finalize spline
    virtual void Finalize(void);

// information methods ---------------------------------------------------------
    /// get CV value for given alpha
    virtual double GetCV(double alpha);

    /// get CV first derivatives for given alpha
    virtual double GetCVFirstDer(double alpha);

// section of private data -----------------------------------------------------
private:
    int                     n;  // number of knots - 1
    CSimpleVector<double>   x;  // alphas       indexing: 0,1,...,n
    CSimpleVector<double>   y;  // CV values    indexing: 0,1,...,n

    // spline
    CSimpleVector<double>   sa, sb, sc, sd;   //indexing: 0,1,...,n
};

//------------------------------------------------------------------------------

typedef std::shared_ptr<CCVSplineNaturalCubic> CCVSplineNaturalCubicPtr;

//------------------------------------------------------------------------------

#endif
