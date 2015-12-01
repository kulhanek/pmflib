#ifndef CVSplineH
#define CVSplineH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <XMLElement.hpp>
#include <SimpleVector.hpp>
#include <FortranMatrix.hpp>

//------------------------------------------------------------------------------

class PMF_PACKAGE CCVSpline {
public:
    CCVSpline(void);
    ~CCVSpline(void);

// setup method ----------------------------------------------------------------
    /// clear all data
    void Clear(void);

    /// allocate memory for data
    void Allocate(int numofknots);

    /// register data point
    void AddPoint(int knotid,double alpha,double cv);

    /// finalize spline
    void Finalize(void);

// information methods ---------------------------------------------------------
    /// get CV value for given aplha
    double GetCV(double alpha);

    /// get CV first derivatives for given alpha
    double GetCVFirstDer(double alpha);

// section of private data -----------------------------------------------------
private:
    int     NumOfKnots;
    double* t;
    double* y;
    double* ypp;
};

//------------------------------------------------------------------------------

#endif
